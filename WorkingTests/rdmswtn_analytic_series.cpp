// Research prototype: analytic threshold integration, no production dispatch.
#include "../src/gaussian.h"
#include "../src/gl_quad.h"
using namespace Rcpp;

static double mills(double x) {
  if(x<26) return std::sqrt(M_PI/2)*std::exp(x*x/2)*std::erfc(x/std::sqrt(2.0));
  double term=1,sum=1;
  for(int k=1;k<=12;++k) {term*=-(2*k-1)/(x*x);sum+=term;}
  return sum/x;
}

// Tail-only initialization fallback. The threshold integration remains analytic;
// this evaluates one bounded posterior-normal expectation at the right endpoint.
static double reflection_seed(double t,double mu,double r,double s,double sv) {
  const double q=s*s+sv*sv*t,S=std::sqrt(t*q),c=s*std::sqrt(t);
  const double m=(sv*sv*r+mu*s*s)/q,sd=s*sv/std::sqrt(q);
  const double lo=std::max(-m/sd,-12.0),hi=12.0;
  if(lo>=hi) return 0;
  const GLRule &gl=gl_get_rule(64); double sum=0;
  for(int j=0;j<64;++j) {
    const double z=(hi+lo)/2+(hi-lo)/2*gl.x[j];
    sum+=gl.w[j]*R::dnorm(z,0,1,false)*mills((r+(m+sd*z)*t)/c);
  }
  return c*R::dnorm((r-mu*t)/S,0,1,false)/S*(hi-lo)/2*sum;
}

static void initial(double t, double mu, double r, double s, double sv,
                    bool pos, double &T, double &H) {
  const double q=s*s+sv*sv*t, S=std::sqrt(t*q);
  const double rho=sv*std::sqrt(t/q), mp=mu+2*r*sv*sv/(s*s);
  const double tilt=2*r*mu/(s*s)+2*r*r*sv*sv/(s*s*s*s);
  if (pos) {
    // Direct rectangles avoid subtraction of nearly equal probabilities.
    T=norm_cdf_2d(mu/sv,(mu*t-r)/S,rho);
    const double h2=(-mp*t-r)/S;
    const double full=std::exp(tilt+R::pnorm(h2,0,1,true,true));
    // For nonnegative mu, take the smaller negative-drift rectangle out
    // of the full Gaussian term. This avoids BVN anticorrelation cancellation.
    if(mu>=0 && (tilt>20 || rho>.9)) H=reflection_seed(t,mu,r,s,sv);
    else H=mu>=0 ? full-std::exp(tilt)*norm_cdf_2d(-mp/sv,h2,rho)
                 : std::exp(tilt)*norm_cdf_2d(mp/sv,h2,-rho);
  } else {
    T=R::pnorm((mu*t-r)/S,0,1,true,false);
    H=std::exp(tilt+R::pnorm((-mp*t-r)/S,0,1,true,true));
  }
}

static double segment(double t,double mu,double c,double h,double s,double sv,
                      bool pos,int order,double *state=nullptr) {
  const double s2=s*s, v=sv*sv, q=s2+v*t, S=std::sqrt(t*q);
  const double noise=s*std::sqrt(t), w=(v*c+mu*s2)/(s*sv*std::sqrt(q));
  const double wp=sv/(s*std::sqrt(q));
  double T,H;
  if(state) {T=state[0]; H=state[1];}
  else initial(t,mu,c,s,sv,pos,T,H);
  double Tend=T,Hend=H;
  double J=R::dnorm((c-mu*t)/S,0,1,false)/S;
  double K=pos ? J*R::dnorm(w,0,1,false) : 0;
  if(pos) J*=R::pnorm(w,0,1,true,false);
  double Q=R::pnorm(-c/noise,0,1,true,false);
  double E=R::dnorm(c/noise,0,1,false);
  const double a0=h*(2*mu/s2+4*v*c/(s2*s2)), a1=4*v*h*h/(s2*s2);
  const double k=1+2*v*t/s2;
  const double d=pos ? 2*sv/s2*R::dnorm(mu/sv,0,1,false) : 0;
  const double l0=-h*(c-mu*t)/(S*S), l1=-h*h/(S*S);
  const double m0=l0-h*w*wp, m1=l1-h*h*wp*wp;
  const double e0=-h*c/(noise*noise), e1=-h*h/(noise*noise);
  double Hp=0,Jp=0,Kp=0,Ep=0;
  double sum=T+H;
  for(int n=0;n<order;++n) {
    const double nt=-h*J/(n+1);
    const double nh=(a0*H+a1*Hp-k*h*J+d*h*Q)/(n+1);
    const double nj=(l0*J+l1*Jp+h*wp*K)/(n+1);
    const double nk=(m0*K+m1*Kp)/(n+1);
    const double nq=-h/noise*E/(n+1);
    const double ne=(e0*E+e1*Ep)/(n+1);
    Hp=H; Jp=J; Kp=K; Ep=E;
    T=nt; H=nh; J=nj; K=nk; Q=nq; E=ne;
    if(state || (n+1)%2==0) sum+=(T+H)/(n+2);
    Tend+=T; Hend+=H;
  }
  if(state) {state[0]=Tend; state[1]=Hend;}
  return sum/(pos ? R::pnorm(mu/sv,0,1,true,false) : 1);
}

// [[Rcpp::export]]
NumericVector backward_cdf(NumericMatrix pars,int order=24) {
  NumericVector out(pars.nrow());
  for(int i=0;i<pars.nrow();++i) {
    const double t=pars(i,0),mu=pars(i,1),b=pars(i,2),A=pars(i,3),s=pars(i,4),sv=pars(i,5);
    const bool pos=pars(i,6)>0;
    double state[2]; initial(t,mu,b,s,sv,pos,state[0],state[1]);
    // Bound dimensionless coefficients and Gaussian variation per step.
    const double growth=A*(2*std::abs(mu)/(s*s)+4*sv*sv*b/std::pow(s,4));
    const int count=std::max(1,(int)std::ceil(std::max(growth/2,2*A/(s*std::sqrt(t)))));
    const double step=A/count; double ans=0;
    for(int j=0;j<count;++j)
      ans+=segment(t,mu,b-j*step,-step,s,sv,pos,order,state);
    out[i]=ans/count;
  }
  return out;
}

// [[Rcpp::export]]
NumericVector series_cdf(NumericMatrix pars,int order=32,int pieces=1) {
  NumericVector out(pars.nrow());
  for(int i=0;i<pars.nrow();++i) {
    double t=pars(i,0),mu=pars(i,1),b=pars(i,2),A=pars(i,3),s=pars(i,4),sv=pars(i,5);
    int count=pieces;
    if(count==0) {
      const double growth=A/2*(2*std::abs(mu)/(s*s)+4*sv*sv*b/std::pow(s,4));
      count=std::max(1,(int)std::ceil(std::max(growth/2,A/(s*std::sqrt(t)))));
    }
    const double h=A/(2*count); double ans=0;
    for(int j=0;j<count;++j) ans+=segment(t,mu,b-A+(2*j+1)*h,h,s,sv,pars(i,6)>0,order);
    out[i]=ans/count;
  }
  return out;
}

// Same direct-rectangle initializer and loop boundary for fair kernel timing.
// [[Rcpp::export]]
NumericVector quadrature_cdf(NumericMatrix pars,int nodes=20,bool original=false) {
  NumericVector out(pars.nrow()); const GLRule &gl=gl_get_rule(nodes);
  for(int i=0;i<pars.nrow();++i) {
    double sum=0;
    for(int j=0;j<nodes;++j) {
      double T,H;
      const double r=pars(i,2)-pars(i,3)/2+pars(i,3)/2*gl.x[j];
      if(original && pars(i,6)>0) {
        // Exact ordinary-range algebra and BVN dispatch from
        // positive_trunc_swtn_cdf_k0; excludes its zero-result fallback.
        const double t=pars(i,0),mu=pars(i,1),s=pars(i,4),sv=pars(i,5);
        const double S=std::sqrt(t*(s*s+sv*sv*t)),rho=sv*std::sqrt(t/(s*s+sv*sv*t));
        const double h1=(mu*t-r)/S,mp=mu+2*r*sv*sv/(s*s),h2=(-mp*t-r)/S;
        T=std::max(0.0,R::pnorm(h1,0,1,true,false)-norm_cdf_2d_hybrid(-mu/sv,h1,-rho));
        H=std::exp(2*r*mu/(s*s)+2*r*r*sv*sv/std::pow(s,4))*
          std::max(0.0,R::pnorm(h2,0,1,true,false)-norm_cdf_2d_hybrid(-mp/sv,h2,rho));
      } else initial(pars(i,0),pars(i,1),r,pars(i,4),pars(i,5),pars(i,6)>0,T,H);
      sum+=gl.w[j]*(T+H);
    }
    out[i]=sum/2/(pars(i,6)>0 ? R::pnorm(pars(i,1)/pars(i,5),0,1,true,false) : 1);
  }
  return out;
}
