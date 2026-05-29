################################################################################
## This file contains:
## - Experimental features. Not yet exported!
################################################################################

## Interface to newton::Tag
Tag <- function(x) {
    if (inherits(x, "advector"))
        LowRankTag(x)
    else
        x
}

TMB_TransformADFunObject <- get("TransformADFunObject", envir = asNamespace("TMB"), inherits = FALSE)
TMB_info <- get("info", envir = asNamespace("TMB"), inherits = FALSE)
## Post vectorization heuristic (Move to TMB if useful)
post_vectorize <- function(ADFun, verbose=FALSE) {
    if (verbose) info1 <- TMB_info(ADFun)
    clear_inv_pos(ADFun$ptr)
    TMB_TransformADFunObject(ADFun, method = "accumulation_tree_split")
    TMB_TransformADFunObject(ADFun, method = "reorder_sub_expressions")
    TMB_TransformADFunObject(ADFun, method = "optimize")
    TMB_TransformADFunObject(ADFun, method = "fuse_and_replay")
    if (verbose) info2 <- TMB_info(ADFun)
    if (verbose) {
        cat("Before:\n")
        print(unlist(info1[-1]))
        cat("After:\n")
        print(unlist(info2[-1]))        
    }
    invisible(ADFun)
}

## Decompose F(x)=L(T(x)) where L is *linear* and *maximal*
## Or apply an 'inner' scale transformation S(F(x)) using linearity of L.
term_split <- function(F, scale=NULL) {
  ptr <- .pointer(environment(F)$mod)
  nodes <- get_term_nodes(ptr)
  L <- .copy(environment(F)$mod)
  T <- .copy(environment(F)$mod)
  vars <- op2var(ptr, nodes)
  setinvIndex(.pointer(L), vars)
  inactivate(.pointer(L), nodes)
  L <- .expose(L)
  L$simplify("eliminate")
  setdepIndex(.pointer(T), vars)
  T <- .expose(T)
  T$simplify("eliminate")
  ## Linearize:
  ## L(x) = L(x0) + L'(x0) * (x - x0)
  zero <- rep(0, length(L$par()))
  L0 <- L(zero)
  J <- L$jacfun(sparse=TRUE)(zero)
  if (!is.null(scale)) {
    J <- AD(J)
    J@x[] <- scale[J@i+1L] * J@x
    ans <- MakeTape(function(x) scale * L0 + AD(J) %*% T(x), T$par())
    return(ans)
  }
  L <- MakeTape(function(x) L0 + AD(J) %*% x, L$par())
  list(T=T, L=L)
}

################################################################################
## GPU codegen
cuda <- new.env() ## Hackable
cuda$include <-
  "#include <cmath>
   template<class T>T sign(const T &x) { return (x > 0) - (x < 0); }
"
cuda$access <-
  "struct access {
      double* x;
      int offset, stride;
      __device__ double& operator[](int i) {
        return x[offset + i * stride];
      }
      __device__ access(double *x) : x(x) {
        stride = gridDim.x * blockDim.x;
        offset = threadIdx.x + blockIdx.x * blockDim.x;
      }
    };
"
cuda$control <-
  '__global__ void forward_global(double* v) {
  forward(v);
}
__global__
void reverse_global(double* v, double* d) {
  reverse(v, d);
}
// Memory alloc
struct device_alloc_t {
  double* value;
  double* deriv;
  int n;
} dev = {0, 0, 0};
// Allocate (n>0) or free (n==0) memory on device
extern "C"
void dev_alloc(int *n) {
  if (dev.value != NULL) cudaFree(dev.value);
  if (dev.deriv != NULL) cudaFree(dev.deriv);
  dev = {0, 0, 0};
  if (n[0] > 0) {
    cudaError_t result = cudaMalloc(&dev.value, n[0] * sizeof(double));
    if (result == cudaSuccess) {
      dev.n = n[0];
    }
  }
}
extern "C"
void clear_deriv() {
  if (dev.deriv == NULL) {
    cudaError_t result = cudaMalloc(&dev.deriv, dev.n * sizeof(double));
  }
  cudaMemset(dev.deriv, 0., dev.n * sizeof(double));
}
// Get / Set array subsets
extern "C"
void get_value(double *x, int *offset, int *n) {
  cudaMemcpy(x, dev.value + offset[0], n[0] * sizeof(double), cudaMemcpyDeviceToHost);
}
extern "C"
void set_value(double *x, int *offset, int *n) {
  cudaMemcpy(dev.value + offset[0], x, n[0] * sizeof(double), cudaMemcpyHostToDevice);
}
extern "C"
void get_deriv(double *x, int *offset, int *n) {
  cudaMemcpy(x, dev.deriv + offset[0], n[0] * sizeof(double), cudaMemcpyDeviceToHost);
}
extern "C"
void set_deriv(double *x, int *offset, int *n) {
  cudaMemcpy(dev.deriv + offset[0], x, n[0] * sizeof(double), cudaMemcpyHostToDevice);
}
extern "C"
void forward_kernel(int* n) {
  forward_global<<<n[0],n[1]>>>(dev.value);
}
extern "C"
void reverse_kernel(int* n) {
  reverse_global<<<n[0],n[1]>>>(dev.value, dev.deriv);
}
'

codegen <- function(F, file=tempfile(), gpu=TRUE, remap=FALSE, ...) {
  names(file) <- basename(file)
  if (gpu) {
    file[] <- paste0(file, ".cu")
    sink(file)
    on.exit(sink(NULL))
    cat(cuda$include)
    if (!remap) {
      cat(cuda$access)
    } else {
      v <- remap_values(.pointer(environment(F)$mod))
      cat(paste("__device__ static int remap[",length(v),"] = {",paste(v, collapse=","),"};\n"))
      cat(sub(" i ", " remap[i] ", cuda$access))
    }
    src_transform(.pointer(environment(F)$mod),
                  config=list(gpu=gpu, ...))
    cat(cuda$control)
  }
  file
}

## Tape -> GPU
gpu_atomic <- function(F, verbose=TRUE, remap=FALSE) {
  inv <- getinvIndex(.pointer(environment(F)$mod))
  dep <- getdepIndex(.pointer(environment(F)$mod))
  if (!all(diff(inv)==1))
    stop("All tape inputs must be consecutive on tape")
  if (!all(diff(dep)==1))
    stop("All tape outputs must be consecutive on tape")
  src <- codegen(F, remap=remap)
  DLL <- names(src)
  dll <- sub(".cu$", ".so", src)
  cmd <- paste("nvcc",
               "--ptxas-options=-v"[verbose],
               "--compiler-options '-fPIC'",
               "-o",
               dll,
               " --shared",
               src)
  system(cmd)
  dyn.load(dll)
  nrep_prev <- 0L
  isInt <- function(x) round(x) == x
  f <- function(x, blksize=32L) {
    nrep <- length(x) / length(inv)
    if (!isInt(nrep)) stop("Invalid length")
    nrep <- as.integer(nrep)
    if (nrep != nrep_prev) {
      value <- rep(F$data.frame()$Value, each=nrep)
      (.C)("dev_alloc", length(value), PACKAGE=DLL)
      (.C)("set_value", value, 0L, length(value), PACKAGE=DLL)
      nrep_prev <<- nrep
    }
    (.C)("set_value", x, inv[1] * nrep, length(inv) * nrep, PACKAGE=DLL)
    if ( !isInt (nrep / blksize) ) blksize <- 1
    parms <- as.integer(c(blksize, nrep/blksize))
    (.C)("forward_kernel", parms, PACKAGE=DLL)
    y <- numeric(length(dep)*nrep)
    ans <- (.C)("get_value", y, dep[1] * nrep, length(dep) * nrep, PACKAGE=DLL)
    ans[[1]]
  }
}
