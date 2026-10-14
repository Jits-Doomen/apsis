#include "apsis/numerics.hpp"
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace apsis::numerics {
namespace {
void finite(double x, const char* what) {
    if (!std::isfinite(x)) throw std::runtime_error(std::string(what) + " produced a non-finite value");
}
State addScaled(const State& a, const State& b, double s) {
    if (a.size()!=b.size()) throw std::invalid_argument("ODE state dimension mismatch");
    State r(a.size()); for (std::size_t i=0;i<a.size();++i) r[i]=a[i]+s*b[i]; return r;
}
State combine(const State& y, double h, const State& a, double ca,
              const State& b, double cb, const State& c, double cc, const State& d, double cd) {
    if (y.size()!=a.size()||y.size()!=b.size()||y.size()!=c.size()||y.size()!=d.size())
        throw std::invalid_argument("ODE derivative dimension mismatch");
    State r(y.size());
    for(std::size_t i=0;i<y.size();++i) r[i]=y[i]+h*(ca*a[i]+cb*b[i]+cc*c[i]+cd*d[i]);
    return r;
}
double simpson(double a,double b,double fa,double fm,double fb) { return (b-a)*(fa+4*fm+fb)/6.0; }
double adaptive(const std::function<double(double)>& f,double a,double b,double fa,double fm,double fb,
                double whole,double tol,std::size_t depth) {
    const double m=(a+b)/2, lm=(a+m)/2, rm=(m+b)/2;
    const double fl=f(lm), fr=f(rm); finite(fl,"integrand"); finite(fr,"integrand");
    const double left=simpson(a,m,fa,fl,fm), right=simpson(m,b,fm,fr,fb);
    const double delta=left+right-whole;
    if(depth==0) {
        if(std::abs(delta)>15*tol) throw std::runtime_error("adaptive Simpson failed to converge");
        return left+right+delta/15.0;
    }
    if(std::abs(delta)<=15*tol) return left+right+delta/15.0;
    return adaptive(f,a,m,fa,fl,fm,left,tol/2,depth-1)+
           adaptive(f,m,b,fm,fr,fb,right,tol/2,depth-1);
}
}
RootResult bisection(const std::function<double(double)>& f,double a,double b,double tol,std::size_t max_iter) {
    if(!(tol>0)||max_iter==0||!std::isfinite(a)||!std::isfinite(b)||a==b) throw std::invalid_argument("invalid bisection settings");
    double fa=f(a), fb=f(b); finite(fa,"function"); finite(fb,"function");
    if(fa==0) return {a,0,0}; if(fb==0) return {b,0,0};
    if(std::signbit(fa)==std::signbit(fb)) throw std::invalid_argument("bisection requires a sign-changing bracket");
    for(std::size_t i=1;i<=max_iter;++i) {
        const double m=a+(b-a)/2, fm=f(m); finite(fm,"function");
        if(fm==0 || std::abs(b-a)/2<=tol) return {m,i,std::abs(fm)};
        if(std::signbit(fa)!=std::signbit(fm)) { b=m; fb=fm; } else { a=m; fa=fm; }
    }
    throw std::runtime_error("bisection did not converge within max_iter");
}
double derivative(const std::function<double(double)>& f,double x,double h) {
    if(!(h>0)||!std::isfinite(x)) throw std::invalid_argument("derivative requires finite x and positive h");
    const double r=(f(x+h)-f(x-h))/(2*h); finite(r,"derivative"); return r;
}
RootResult newton(const std::function<double(double)>& f,double x,double tol,std::size_t max_iter,double dh) {
    if(!(tol>0)||max_iter==0) throw std::invalid_argument("invalid Newton settings");
    for(std::size_t i=0;i<max_iter;++i) {
        const double fx=f(x); finite(fx,"function");
        if(std::abs(fx)<=tol) return {x,i,std::abs(fx)};
        const double d=derivative(f,x,dh);
        if(std::abs(d)<1e-15) throw std::runtime_error("Newton derivative is near zero");
        const double next=x-fx/d; finite(next,"Newton iteration");
        if(std::abs(next-x)<=tol) { const double r=std::abs(f(next)); return {next,i+1,r}; }
        x=next;
    }
    throw std::runtime_error("Newton method did not converge within max_iter");
}
double integrateAdaptiveSimpson(const std::function<double(double)>& f,double a,double b,double tol,std::size_t depth) {
    if(!(tol>0)||depth==0||!std::isfinite(a)||!std::isfinite(b)) throw std::invalid_argument("invalid integration settings");
    if(a==b) return 0; if(b<a) return -integrateAdaptiveSimpson(f,b,a,tol,depth);
    const double m=(a+b)/2, fa=f(a), fm=f(m), fb=f(b);
    finite(fa,"integrand"); finite(fm,"integrand"); finite(fb,"integrand");
    return adaptive(f,a,b,fa,fm,fb,simpson(a,b,fa,fm,fb),tol,depth);
}
OdeResult integrateRK4(const OdeFunction& f,double t,const State& y,double t1,double h) {
    if(y.empty()||!std::isfinite(t)||!std::isfinite(t1)||!(h>0)) throw std::invalid_argument("invalid RK4 initial conditions");
    if(t1<t) throw std::invalid_argument("RK4 currently integrates forward in time only");
    State state=y; std::size_t steps=0;
    while(t<t1) {
        const double dt=std::min(h,t1-t);
        const State k1=f(t,state);
        const State k2=f(t+dt/2,addScaled(state,k1,dt/2));
        const State k3=f(t+dt/2,addScaled(state,k2,dt/2));
        const State k4=f(t+dt,addScaled(state,k3,dt));
        state=combine(state,dt,k1,1.0/6,k2,1.0/3,k3,1.0/3,k4,1.0/6);
        for(double v:state) finite(v,"ODE integration");
        t+=dt; ++steps;
        if(steps>10000000) throw std::runtime_error("RK4 step budget exceeded");
    }
    return {t,state,steps,0};
}

OdeResult integrateRK45(const OdeFunction& f,double t,const State& y0,double t1,double h,
                        double atol,double rtol,double min_h,double max_h) {
    if(y0.empty()||!std::isfinite(t)||!std::isfinite(t1)||t1<t||h<=0||atol<=0||rtol<0||min_h<=0||max_h<=0)
        throw std::invalid_argument("invalid RK45 settings or initial conditions");
    State y=y0; std::size_t accepted=0,rejected=0,attempts=0;
    auto eval=[&](double time,const State& state){State k=f(time,state);if(k.size()!=state.size())throw std::runtime_error("ODE derivative dimension mismatch");for(double v:k)finite(v,"ODE derivative");return k;};
    while(t<t1) {
        if(++attempts>10000000)throw std::runtime_error("RK45 attempt budget exceeded");
        h=std::min({h,max_h,t1-t});
        if(h<min_h && t1-t>min_h)throw std::runtime_error("RK45 step underflow");
        const State k1=eval(t,y);
        const State y2=combine(y,h,k1,1.0/5,State(y.size()),0,State(y.size()),0,State(y.size()),0);
        const State k2=eval(t+h/5,y2);
        State q=y;for(std::size_t i=0;i<q.size();++i)q[i]=y[i]+h*(3.0/40*k1[i]+9.0/40*k2[i]);
        const State k3=eval(t+3*h/10,q);
        for(std::size_t i=0;i<q.size();++i)q[i]=y[i]+h*(44.0/45*k1[i]-56.0/15*k2[i]+32.0/9*k3[i]);
        const State k4=eval(t+4*h/5,q);
        for(std::size_t i=0;i<q.size();++i)q[i]=y[i]+h*(19372.0/6561*k1[i]-25360.0/2187*k2[i]+64448.0/6561*k3[i]-212.0/729*k4[i]);
        const State k5=eval(t+8*h/9,q);
        for(std::size_t i=0;i<q.size();++i)q[i]=y[i]+h*(9017.0/3168*k1[i]-355.0/33*k2[i]+46732.0/5247*k3[i]+49.0/176*k4[i]-5103.0/18656*k5[i]);
        const State k6=eval(t+h,q);
        State y5=y,y4=y;
        for(std::size_t i=0;i<y.size();++i) {
            y5[i]=y[i]+h*(35.0/384*k1[i]+500.0/1113*k3[i]+125.0/192*k4[i]-2187.0/6784*k5[i]+11.0/84*k6[i]);
        }
        const State k7=eval(t+h,y5);
        double err=0;
        for(std::size_t i=0;i<y.size();++i) {
            y4[i]=y[i]+h*(5179.0/57600*k1[i]+7571.0/16695*k3[i]+393.0/640*k4[i]-92097.0/339200*k5[i]+187.0/2100*k6[i]+1.0/40*k7[i]);
            const double scale=atol+rtol*std::max(std::abs(y[i]),std::abs(y5[i]));
            err=std::max(err,std::abs(y5[i]-y4[i])/scale);
        }
        if(err<=1.0) { t+=h;y=std::move(y5);++accepted; }
        else ++rejected;
        const double factor=err==0?5.0:std::clamp(0.9*std::pow(err,-0.2),0.1,5.0);
        h=std::max(min_h,h*factor);
    }
    return {t,y,accepted,rejected};
}
}
