// $Id: Funcs.cc,v 1.1 2009/06/11 20:31:19 loizides Exp $

#include "MitCommon/MathTools/interface/Funcs.h"
#include <TF1.h>

using namespace mithep;

//--------------------------------------------------------------------------------------------------
Double_t Funcs::BreitGaus(Double_t x, Double_t m, Double_t mwidth, Double_t msig, 
                            Double_t fintf, Double_t xmin, Double_t xmax)
{
  // Breit-Wigner convoluted with Gaussian.
  // Fit parameters:
  //  m: Most probable location) Breit mean
  //  mwidth: Width (scale) Breit-Wigner
  //  msig: Width (sigma) of convoluted Gaussian function
  //  fintf: interference fraction between Z and gamma
  //  xmin: minimum x (for normalization)
  //  xmax: maximum x (for normalization)

  // Numeric constants
  const Double_t invsq2pi = 0.3989422804014; //(2pi)^(-1/2)

  // Control constants
  const Double_t np = 300; // number of convolution steps
  const Double_t sc = 3.0; // convolution extends to +-sc Gaussian sigmas

  // Range of convolution integral
  const Double_t xlow = x - sc * msig;
  const Double_t xupp = x + sc * msig;
  const Double_t step = (xupp-xlow) / np;

  // Convolution integral of Breit and Gaussian by sum
  Double_t sum0 = 0.0;
  if (fintf<1) {
    for(Double_t i=1.0; i<=np/2; i++) {
      Double_t xx = xlow + (i-0.5) * step;
      Double_t fland = BreitWignerZ(xx,m,mwidth);
      sum0 += fland * TMath::Gaus(x,xx,msig);
      xx = xupp - (i-0.5) * step;
      fland = BreitWignerZ(xx,m,mwidth);
      sum0 += fland * TMath::Gaus(x,xx,msig);
    }
    sum0 = (step * sum0 * invsq2pi / msig);
  }

  Double_t sum1 = 0.0;
  if (fintf>0) {
    for(Double_t i=1.0; i<=np/2; i++) {
      Double_t xx = xlow + (i-0.5) * step;
      Double_t fland = BreitWignerGamma(xx,m,mwidth);
      sum1 += fland * TMath::Gaus(x,xx,msig);
      xx = xupp - (i-0.5) * step;
      fland = BreitWignerGamma(xx,m,mwidth);
      sum1 += fland * TMath::Gaus(x,xx,msig);
    }
    sum1 = (step * sum1 * invsq2pi / msig / (xmax-xmin));
  }
  return ((1.-fintf)*sum0+fintf*sum1);
}

//--------------------------------------------------------------------------------------------------
Double_t Funcs::BreitGaus(Double_t *x, Double_t *par)
{
  // Breit-Wigner convoluted with Gaussian.
  // Fit parameters:
  //  norm: Normalization factor
  //  m: Most probable location) Breit mean
  //  mwidth: Width (scale) Breit-Wigner
  //  msig: Width (sigma) of convoluted Gaussian function
  //  fintf: interference fraction between Z and gamma
  //  xmin: minimum x (for normalization)
  //  xmax: maximum x (for normalization)

  Double_t xx     = x[0];
  Double_t norm   = par[0];
  Double_t m      = par[1];
  Double_t mwidth = par[2];
  Double_t msig   = par[3];
  Double_t fintf  = par[4];
  Double_t xmin   = par[5];
  Double_t xmax   = par[6];
  Double_t ret = norm * BreitGaus(xx, m, mwidth, msig, fintf, xmin, xmax);
  return ret;
}

//--------------------------------------------------------------------------------------------------
Double_t Funcs::BreitWignerZ(Double_t x, Double_t mean, Double_t gamma)
{
  // Breit-Wigner for Z peak.

  Double_t bw =  (gamma*mean*mean)/
                 ((x*x-mean*mean)*(x*x-mean*mean) + x*x*x*x*(gamma*gamma)/(mean*mean));
  return bw / (TMath::Pi()/2.0);
}

//--------------------------------------------------------------------------------------------------
Double_t Funcs::BreitWignerZ(Double_t *x, Double_t *par)
{
  // Breit-Wigner for Z peak.

  Double_t xx    = x[0];
  Double_t norm  = par[0];
  Double_t mean  = par[1];
  Double_t gamma = par[2];

  Double_t ret = norm * BreitWignerZ(xx,mean,gamma);
  return ret;
}

//--------------------------------------------------------------------------------------------------
Double_t Funcs::BreitWignerGamma(Double_t x, Double_t mean, Double_t gamma)
{
  // Breit-Wigner for gamma interference.

  Double_t bw =  ((x*x-mean*mean)*(x*x-mean*mean))/
                 ((x*x-mean*mean)*(x*x-mean*mean) + x*x*x*x*(gamma*gamma)/(mean*mean));
  return bw;
}

//--------------------------------------------------------------------------------------------------
Double_t Funcs::BreitWignerGamma(Double_t *x, Double_t *par)
{
  // Breit-Wigner for gamma interference.

  Double_t xx    = x[0];
  Double_t norm  = par[0];
  Double_t mean  = par[1];
  Double_t gamma = par[2];

  Double_t ret = norm * BreitWignerGamma(xx,mean,gamma);
  return ret;
}

//--------------------------------------------------------------------------------------------------
Double_t Funcs::ExpRange(Double_t mass, Double_t lambda, Double_t xmin, Double_t xmax)
{
  // Probability of an exponential in a given interval.

  if (mass < xmin || mass > xmax)
    return 0.0;

  if (lambda == 0.0)
    return 1.0/(xmax-xmin);

  Double_t xmed  = 0.5*(xmin+xmax);
  Double_t integ = lambda * TMath::Exp(-lambda*xmed) /
    (TMath::Exp(-lambda*xmin)-TMath::Exp(-lambda*xmax));

  return integ*TMath::Exp(-lambda*(mass-xmed));
}

//--------------------------------------------------------------------------------------------------
Double_t Funcs::ExpRange(Double_t *x, Double_t *par)
{
  // Probability of an exponential in a given interval.

  Double_t xx     = x[0];
  Double_t norm   = par[0];
  Double_t lambda = par[1];
  Double_t xmin   = par[2];
  Double_t xmax   = par[3];

  Double_t ret = norm * ExpRange(xx,lambda,xmin,xmax);
  return ret;
}

//--------------------------------------------------------------------------------------------------
Double_t Funcs::ZLineShapePlusBkg(Double_t *x, Double_t *par)
{ 
  // Z line shape (BreitGaus) plus background (ExpRange).

  // Parameters:
  //  0: overall norm
  //  1: relative fraction of signal and background
  //  2: relative fraction of exponential and polynomial
  //  3: most probable location) Breit mean
  //  4: width (scale) Breit-Wigner
  //  5: width (sigma) of convoluted Gaussian function
  //  6: interference fraction between Z and gamma
  //  7: lambda of exponential
  //  8: minimum x (for normalization)
  //  9: maximum x (for normalization)

  Double_t xmin=par[8];
  Double_t xmax=par[9];
  Double_t s = (1.0 - par[1]) * BreitGaus(x[0], par[3], par[4], par[5], par[6], xmin, xmax);
  Double_t b = par[1] * ((1.0 - par[2]) * ExpRange(x[0], par[7], xmin, xmax) + par[2]/(xmax-xmin));

  return (s+b)*par[0];
}

//--------------------------------------------------------------------------------------------------
TF1 *Funcs::CreateZLineShapePlusBkg(Double_t norm, Double_t xmin, Double_t xmax, const char *n)
{
  // Create TF1 for Z line shape and background.

  TF1 *f = new TF1(n,mithep::Funcs::ZLineShapePlusBkg,xmin,xmax,10);
  f->SetParName(0,"norm");
  f->SetParName(1,"f_sb");
  f->SetParName(2,"f_eb");
  f->SetParName(3,"mass");
  f->SetParName(4,"width");
  f->SetParName(5,"gsig");
  f->SetParName(6,"f_int");
  f->SetParName(7,"lambda");
  f->SetParName(8,"min");
  f->SetParName(9,"max");

  f->SetParameters(norm,0.5,0.5,91.,1.,1.,0.,0.1,xmin,xmax);
  f->SetParLimits(0,0,1e6);
  f->SetParLimits(1,0,1);
  f->SetParLimits(2,0,1);
  f->SetParLimits(3,xmin,xmax);
  f->SetParLimits(4,0.1,5);
  f->SetParLimits(5,0.1,3);
  f->SetParLimits(6,0,1);
  f->SetParLimits(7,0,10);
  f->FixParameter(8,xmin);
  f->FixParameter(9,xmax);

  return f;
}
