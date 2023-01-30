// Autogenerated - do not edit
#include <Rcpp.h>
#include "TMB.h"
typedef TMBad::ad_aug ad;
ad* adptr(const Rcpp::ComplexVector &x);
Rcpp::ComplexVector& as_advector(Rcpp::ComplexVector &x);
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dexp ( Rcpp::ComplexVector x, Rcpp::ComplexVector rate, bool give_log )
{
int n1=x.size();
int n2=rate.size();
int nmax = std::max({n1, n2});
int nmin = std::min({n1, n2});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(rate);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dexp(X1[i % n1], X2[i % n2], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dweibull ( Rcpp::ComplexVector x, Rcpp::ComplexVector shape, Rcpp::ComplexVector scale, bool give_log )
{
int n1=x.size();
int n2=shape.size();
int n3=scale.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(shape); const ad* X3 = adptr(scale);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dweibull(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dbinom ( Rcpp::ComplexVector x, Rcpp::ComplexVector size, Rcpp::ComplexVector prob, bool give_log )
{
int n1=x.size();
int n2=size.size();
int n3=prob.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(size); const ad* X3 = adptr(prob);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dbinom(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dbinom_robust ( Rcpp::ComplexVector x, Rcpp::ComplexVector size, Rcpp::ComplexVector logit_p, bool give_log )
{
int n1=x.size();
int n2=size.size();
int n3=logit_p.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(size); const ad* X3 = adptr(logit_p);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dbinom_robust(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dbeta ( Rcpp::ComplexVector x, Rcpp::ComplexVector shape1, Rcpp::ComplexVector shape2, bool give_log )
{
int n1=x.size();
int n2=shape1.size();
int n3=shape2.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(shape1); const ad* X3 = adptr(shape2);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dbeta(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_df ( Rcpp::ComplexVector x, Rcpp::ComplexVector df1, Rcpp::ComplexVector df2, bool give_log )
{
int n1=x.size();
int n2=df1.size();
int n3=df2.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(df1); const ad* X3 = adptr(df2);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = df(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dlogis ( Rcpp::ComplexVector x, Rcpp::ComplexVector location, Rcpp::ComplexVector scale, bool give_log )
{
int n1=x.size();
int n2=location.size();
int n3=scale.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(location); const ad* X3 = adptr(scale);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dlogis(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dsn ( Rcpp::ComplexVector x, Rcpp::ComplexVector alpha, bool give_log )
{
int n1=x.size();
int n2=alpha.size();
int nmax = std::max({n1, n2});
int nmin = std::min({n1, n2});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(alpha);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dsn(X1[i % n1], X2[i % n2], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dt ( Rcpp::ComplexVector x, Rcpp::ComplexVector df, bool give_log )
{
int n1=x.size();
int n2=df.size();
int nmax = std::max({n1, n2});
int nmin = std::min({n1, n2});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(df);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dt(X1[i % n1], X2[i % n2], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dSHASHo ( Rcpp::ComplexVector x, Rcpp::ComplexVector mu, Rcpp::ComplexVector sigma, Rcpp::ComplexVector nu, Rcpp::ComplexVector tau, bool give_log )
{
int n1=x.size();
int n2=mu.size();
int n3=sigma.size();
int n4=nu.size();
int n5=tau.size();
int nmax = std::max({n1, n2, n3, n4, n5});
int nmin = std::min({n1, n2, n3, n4, n5});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(mu); const ad* X3 = adptr(sigma); const ad* X4 = adptr(nu); const ad* X5 = adptr(tau);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dSHASHo(X1[i % n1], X2[i % n2], X3[i % n3], X4[i % n4], X5[i % n5], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dtweedie ( Rcpp::ComplexVector x, Rcpp::ComplexVector mu, Rcpp::ComplexVector phi, Rcpp::ComplexVector p, bool give_log )
{
int n1=x.size();
int n2=mu.size();
int n3=phi.size();
int n4=p.size();
int nmax = std::max({n1, n2, n3, n4});
int nmin = std::min({n1, n2, n3, n4});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(mu); const ad* X3 = adptr(phi); const ad* X4 = adptr(p);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dtweedie(X1[i % n1], X2[i % n2], X3[i % n3], X4[i % n4], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dnorm ( Rcpp::ComplexVector x, Rcpp::ComplexVector mean, Rcpp::ComplexVector sd, bool give_log )
{
int n1=x.size();
int n2=mean.size();
int n3=sd.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(mean); const ad* X3 = adptr(sd);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dnorm(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dnbinom ( Rcpp::ComplexVector x, Rcpp::ComplexVector size, Rcpp::ComplexVector prob, bool give_log )
{
int n1=x.size();
int n2=size.size();
int n3=prob.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(size); const ad* X3 = adptr(prob);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dnbinom(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dnbinom2 ( Rcpp::ComplexVector x, Rcpp::ComplexVector mu, Rcpp::ComplexVector var, bool give_log )
{
int n1=x.size();
int n2=mu.size();
int n3=var.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(mu); const ad* X3 = adptr(var);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dnbinom2(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dnbinom_robust ( Rcpp::ComplexVector x, Rcpp::ComplexVector log_mu, Rcpp::ComplexVector log_var_minus_mu, bool give_log )
{
int n1=x.size();
int n2=log_mu.size();
int n3=log_var_minus_mu.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(log_mu); const ad* X3 = adptr(log_var_minus_mu);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dnbinom_robust(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dpois ( Rcpp::ComplexVector x, Rcpp::ComplexVector lambda, bool give_log )
{
int n1=x.size();
int n2=lambda.size();
int nmax = std::max({n1, n2});
int nmin = std::min({n1, n2});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(lambda);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dpois(X1[i % n1], X2[i % n2], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dgamma ( Rcpp::ComplexVector x, Rcpp::ComplexVector shape, Rcpp::ComplexVector scale, bool give_log )
{
int n1=x.size();
int n2=shape.size();
int n3=scale.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(shape); const ad* X3 = adptr(scale);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dgamma(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dlgamma ( Rcpp::ComplexVector x, Rcpp::ComplexVector shape, Rcpp::ComplexVector scale, bool give_log )
{
int n1=x.size();
int n2=shape.size();
int n3=scale.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(shape); const ad* X3 = adptr(scale);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dlgamma(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_dzipois ( Rcpp::ComplexVector x, Rcpp::ComplexVector lambda, Rcpp::ComplexVector zip, bool give_log )
{
int n1=x.size();
int n2=lambda.size();
int n3=zip.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(x); const ad* X2 = adptr(lambda); const ad* X3 = adptr(zip);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = dzipois(X1[i % n1], X2[i % n2], X3[i % n3], give_log);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_pnorm_approx ( Rcpp::ComplexVector q )
{
int n1=q.size();
int nmax = std::max({n1});
int nmin = std::min({n1});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(q);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = pnorm_approx(X1[i % n1]);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_pnorm ( Rcpp::ComplexVector q, Rcpp::ComplexVector mean , Rcpp::ComplexVector sd  )
{
int n1=q.size();
int n2=mean .size();
int n3=sd .size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(q); const ad* X2 = adptr(mean ); const ad* X3 = adptr(sd );
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = pnorm(X1[i % n1], X2[i % n2], X3[i % n3]);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_pgamma ( Rcpp::ComplexVector q, Rcpp::ComplexVector shape, Rcpp::ComplexVector scale  )
{
int n1=q.size();
int n2=shape.size();
int n3=scale .size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(q); const ad* X2 = adptr(shape); const ad* X3 = adptr(scale );
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = pgamma(X1[i % n1], X2[i % n2], X3[i % n3]);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_ppois ( Rcpp::ComplexVector q, Rcpp::ComplexVector lambda )
{
int n1=q.size();
int n2=lambda.size();
int nmax = std::max({n1, n2});
int nmin = std::min({n1, n2});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(q); const ad* X2 = adptr(lambda);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = ppois(X1[i % n1], X2[i % n2]);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_pexp ( Rcpp::ComplexVector q, Rcpp::ComplexVector rate )
{
int n1=q.size();
int n2=rate.size();
int nmax = std::max({n1, n2});
int nmin = std::min({n1, n2});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(q); const ad* X2 = adptr(rate);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = pexp(X1[i % n1], X2[i % n2]);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_pweibull ( Rcpp::ComplexVector q, Rcpp::ComplexVector shape, Rcpp::ComplexVector scale )
{
int n1=q.size();
int n2=shape.size();
int n3=scale.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(q); const ad* X2 = adptr(shape); const ad* X3 = adptr(scale);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = pweibull(X1[i % n1], X2[i % n2], X3[i % n3]);
return as_advector(ans);
}
// [[Rcpp::export]]
Rcpp::ComplexVector distr_pbeta ( Rcpp::ComplexVector q, Rcpp::ComplexVector shape1, Rcpp::ComplexVector shape2 )
{
int n1=q.size();
int n2=shape1.size();
int n3=shape2.size();
int nmax = std::max({n1, n2, n3});
int nmin = std::min({n1, n2, n3});
int n = (nmin == 0 ? 0 : nmax);
Rcpp::ComplexVector ans(n);
const ad* X1 = adptr(q); const ad* X2 = adptr(shape1); const ad* X3 = adptr(shape2);
ad* Y = adptr(ans);
for (int i=0; i<n; i++) Y[i] = pbeta(X1[i % n1], X2[i % n2], X3[i % n3]);
return as_advector(ans);
}
