#include "RTMB.h"

#define EPSILON DBL_EPSILON

/* R_zeroin2() is faster for "expensive" f(), in those typical cases where
 *             f(ax) and f(bx) are available anyway : */

template <class FUNC>
double R_zeroin2(
    double ax,
    double bx,
    double fa,
    double fb,
    FUNC &f,
    void *info, // Extra passed to f
    double *Tol,
    int *Maxit)
{
    double a,b,c, fc;			/* Abscissae, descr. see above,  f(c) */
    double tol;
    int maxit;

    a = ax;  b = bx;
    c = a;   fc = fa;
    maxit = *Maxit + 1; tol = * Tol;

    /* First test if we have found a root at an endpoint */
    if(fa == 0.0) {
	*Tol = 0.0;
	*Maxit = 0;
	return a;
    }
    if(fb ==  0.0) {
	*Tol = 0.0;
	*Maxit = 0;
	return b;
    }

    while(maxit--)		/* Main iteration loop	*/
    {
	double prev_step = b-a;		/* Distance from the last but one
					   to the last approximation	*/
	double tol_act;			/* Actual tolerance		*/
	double p;			/* Interpolation step is calcu- */
	double q;			/* lated in the form p/q; divi-
					 * sion operations is delayed
					 * until the last moment	*/
	double new_step;		/* Step at this iteration	*/

	if( fabs(fc) < fabs(fb) )
	{				/* Swap data for b to be the	*/
	    a = b;  b = c;  c = a;	/* best approximation		*/
	    fa=fb;  fb=fc;  fc=fa;
	}
	tol_act = 2*EPSILON*fabs(b) + tol/2;
	new_step = (c-b)/2;

	if( fabs(new_step) <= tol_act || fb == (double)0 )
	{
	    *Maxit -= maxit;
	    *Tol = fabs(c-b);
	    return b;			/* Acceptable approx. is found	*/
	}

	/* Decide if the interpolation can be tried	*/
	if( fabs(prev_step) >= tol_act	/* If prev_step was large enough*/
	    && fabs(fa) > fabs(fb) ) {	/* and was in true direction,
					 * Interpolation may be tried	*/
	    double t1,cb,t2;
	    cb = c-b;
	    if( a==c ) {		/* If we have only two distinct	*/
					/* points linear interpolation	*/
		t1 = fb/fa;		/* can only be applied		*/
		p = cb*t1;
		q = 1.0 - t1;
	    }
	    else {			/* Quadric inverse interpolation*/

		q = fa/fc;  t1 = fb/fc;	 t2 = fb/fa;
		p = t2 * ( cb*q*(q-t1) - (b-a)*(t1-1.0) );
		q = (q-1.0) * (t1-1.0) * (t2-1.0);
	    }
	    if( p>(double)0 )		/* p was calculated with the */
		q = -q;			/* opposite sign; make p positive */
	    else			/* and assign possible minus to	*/
		p = -p;			/* q				*/

	    if( p < (0.75*cb*q-fabs(tol_act*q)/2) /* If b+p/q falls in [b,c]*/
		&& p < fabs(prev_step*q/2) )	/* and isn't too large	*/
		new_step = p/q;			/* it is accepted
						 * If p/q is too large then the
						 * bisection procedure can
						 * reduce [b,c] range to more
						 * extent */
	}

	if( fabs(new_step) < tol_act) {	/* Adjust the step to be not less*/
	    if( new_step > (double)0 )	/* than tolerance		*/
		new_step = tol_act;
	    else
		new_step = -tol_act;
	}
	a = b;	fa = fb;			/* Save the previous approx. */
	b += new_step;	fb = f(b, info);	/* Do step to a new approxim. */
	if( (fb > 0 && fc > 0) || (fb < 0 && fc < 0) ) {
	    /* Adjust c for it to have a sign opposite to that of b */
	    c = a;  fc = fa;
	}

    }
    /* failed! */
    *Tol = fabs(c-b);
    *Maxit = -1;
    return b;
}


struct RootOp : TMBad::global::DynamicOperator< -1 , -1 > {
  static const bool have_input_size_output_size = true;
  static const bool add_forward_replay_copy = true;
  TMBAD_SHARED_PTR<TMBad::ADFun<> > f;
  double a, b;
  double tol;
  int maxiter;
  RootOp(const TMBad::ADFun<> &F, double a, double b, Rcpp::List cfg) :
    f(std::make_shared<TMBad::ADFun<> >(F)),
    a(a), b(b),
    tol(REAL(cfg["tol"])[0]),
    maxiter(INTEGER(cfg["maxiter"])[0]) { }
  struct feval {
    TMBad::ADFun<> &f;
    double operator() (double x, void *info) {
      f.glob.value_inv(0) = x;
      f.glob.forward(); // FIXME: From start
      return f.glob.value_dep(0);
    }
  };
  TMBad::Index input_size() const { return (*f).Domain() - 1 ; } // parms
  TMBad::Index output_size() const { return 1 ; } // root, but maybe more
  template<class Type>
  void forward(TMBad::ForwardArgs<Type> &args) {
    TMBAD_ASSERT(false);
  }
  void forward(TMBad::ForwardArgs<double> &args) {
    int n = input_size();
    std::vector<double> t = args.x_segment(0, n);
    std::vector<double> at = {a};
    at.insert(at.end(), t.begin(), t.end());
    double fa = (*f)(at)[0];
    at[0] = b;
    double fb = (*f)(at)[0];
    void* info = NULL;
    feval F = {*f};
    args.y(0) = R_zeroin2(a,b,fa,fb,F,info,&tol,&maxiter);
  }
  template<class Type> void reverse(TMBad::ReverseArgs<Type> &args) {
    size_t n = input_size();
    std::vector<Type> t = args.x_segment(0, n);
    std::vector<Type> xt = args.y_segment(0, 1);
    xt.insert(xt.end(), t.begin(), t.end());
    std::vector<Type> one = {1};
    std::vector<Type> J = (*f).Jacobian(xt, one); // n+1
    Type s = -(args.dy(0) / J[0]);
    for (size_t i = 0; i<n; i++)
      args.dx(i) += s * J[i+1];
  }
  void reverse(TMBad::ReverseArgs<TMBad::Writer> &args) {
    TMBAD_ASSERT(false);
  }
  // 'f' can be used to uniqely identify this operator
  static const bool have_custom_identifier = true;
  void* custom_identifier() { return &(*f); }
  const char* op_name() { return "RootOp"; }
};


// [[Rcpp::export]]
ADrep ad_uniroot(Rcpp::XPtr<TMBad::ADFun<> > adf, double a, double b, ADrep parms, Rcpp::List cfg) {
  std::vector<ad> x(parms.adptr(), parms.adptr() + parms.size());
  std::vector<ad> y = TMBad::global::Complete<RootOp>(*adf, a, b, cfg) (x);
  return ADrep(y.data(), y.data() + y.size());
}
