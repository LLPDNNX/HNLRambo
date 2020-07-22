#include <iostream>
#include <vector>
#include <math.h>
#include <cmath>
#include <stdlib.h>

#include "HNLRambo/Rambo/interface/Rambo.h"

Rambo::Rambo(size_t seed):
    rndEngine_(seed),
    dist_(0.0,1.0)
{
}


vector<Rambo4Vector> Rambo::generate(
    double scale, 
    const std::vector<double>& xm, 
    double& wt
)
{
/**********************************************************************
 *                       rambo                                         *
 *    ra(ndom)  m(omenta)  b(eautifully)  o(rganized)                  *
 *                                                                     *
 *    a democratic multi-particle phase space generator                *
 *    authors:  s.d. ellis,  r. kleiss,  w.j. stirling                 *
 *    this is version 1.0 -  written by r. kleiss                      *
 *    -- adjusted by hans kuijf, weights are logarithmic (20-08-90)    *
 *                                                                     *
 *    n  = number of particles                                         *
 *    scale = total centre-of-mass energy                                 *
 *    xm = particle masses ( dim=nexternal-nincoming )                 *
 *    p  = particle momenta ( dim=(4,nexternal-nincoming) )            *
 *    wt = weight of the event                                         *
 ***********************************************************************/
  const size_t n = xm.size();
  std::vector<Rambo4Vector> q(n), p(n);
  std::vector<double> z(n,0), r(4,0), b(3,0), p2(n,0), xm2(n,0), e(n,0), v(n,0);
  
  
  constexpr int itmax = 6;
  constexpr double acc = 1e-14;
  constexpr double twopi=8.*atan(1.);
  constexpr double po2log=log(twopi/4.);

// initialization step: factorials for the phase space weight

    z[1]=po2log;
    for(size_t k=2;k < n;k++)
      z[k]=z[k-1]+po2log-2.*log(double(k-1));
    for(size_t k=2;k < n;k++)
      z[k]=(z[k]-log(double(k)));
  
// check on the number of particles
  if(n<1 || n>101){
    throw std::runtime_error("Too few or many particles: "+std::to_string(n));
  }
// check whether total energy is sufficient; count nonzero masses
  double xmt=0.;
  int nm=0;
  for(size_t i=0; i<n; i++){
    if(xm[i]!=0.) nm=nm+1;
    xmt=xmt+std::fabs(xm[i]);
  }
  if (xmt>scale){
    throw std::runtime_error("Too low energy: "+std::to_string(scale)+" needed "+std::to_string(xmt));

  }
// the parameter values are now accepted

// generate n massless momenta in infinite phase space
  for(size_t i=0; i<n;i++){
    double r1=rnd();
    double c=2.*r1-1.;
    double s=sqrt(1.-c*c);
    double f=twopi*rnd();
    r1=rnd();
    double r2=rnd();
    q[i][0]=-log(r1*r2);
    q[i][3]=q[i][0]*c;
    q[i][2]=q[i][0]*s*cos(f);
    q[i][1]=q[i][0]*s*sin(f);
  }
// calculate the parameters of the conformal transformation
  for (size_t i=0;i < 4;i++)
    r[i]=0.;
  for(size_t i=0;i < n;i++){
    for (size_t k=0; k<4; k++)
      r[k]=r[k]+q[i][k];
  }
  double rmas=sqrt(pow(r[0],2)-pow(r[3],2)-pow(r[2],2)-pow(r[1],2));
  for(size_t k=1;k < 4; k++)
    b[k-1]=-r[k]/rmas;
  double g=r[0]/rmas;
  double a=1./(1.+g);
  double x=scale/rmas;

// transform the q's conformally into the p's
  for(size_t i=0; i< n;i++){
    double bq=b[0]*q[i][1]+b[1]*q[i][2]+b[2]*q[i][3];
    for (int k=1;k<4;k++)
      p[i][k]=x*(q[i][k]+b[k-1]*(q[i][0]+a*bq));
    p[i][0]=x*(g*q[i][0]+bq);
  }

// calculate weight and possible warnings
  wt=po2log;
  if(n!=2) wt=(2.*n-4.)*log(scale)+z[n-1];
  if(wt<-180.){
    std::cout << "Too small wt, risk for underflow: " << wt << std::endl;
  }
  if(wt> 174.){
    std::cout << "Too large wt, risk for overflow: " << wt << std::endl;
  }

// return for weighted massless momenta
  if(nm==0){
// return log of weight
    return p;
  }

// massive particles: rescale the momenta by a factor x
  double xmax=sqrt(1.-pow(xmt/scale, 2));
  for(size_t i=0;i < n; i++){
    xm2[i]=pow(xm[i],2);
    p2[i]=pow(p[i][0],2);
  }
  int iter=0;
  x=xmax;
  double accu=scale*acc;
  while(true){
    double f0=-scale;
    double g0=0.;
    double x2=x*x;
    for(size_t i=0; i < n; i++){
      e[i]=sqrt(xm2[i]+x2*p2[i]);
      f0=f0+e[i];
      g0=g0+p2[i]/e[i];
    }
    if(abs(f0)<=accu) break;
    iter=iter+1;
    if(iter>itmax){
      std::cout << "Too many iterations without desired accuracy: " << itmax << std::endl;
      break;
    }
    x=x-f0/(x*g0);
  }
  for(size_t i=0;i < n;i++){
    v[i]=x*p[i][0];
    for(size_t k=1;k < 4; k++)
      p[i][k]=x*p[i][k];
    p[i][0]=e[i];
  }

// calculate the mass-effect weight factor
  double wt2=1.;
  double wt3=0.;
  for(size_t i=0;i < n; i++){
    wt2=wt2*v[i]/e[i];
    wt3=wt3+pow(v[i],2)/e[i];
  }
  double wtm=(2.*n-3.)*log(x)+log(wt2/wt3*scale);

// return for  weighted massive momenta
  wt=wt+wtm;
  if(wt<-180.){
    std::cout << "Too small wt, risk for underflow: " << wt << std::endl;
    
  }
  if(wt> 174.){
    std::cout << "Too large wt, risk for overflow: " << wt << std::endl;
  }
// return log of weight
  return p;
}


