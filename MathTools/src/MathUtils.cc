#include "MitCommon/MathTools/interface/MathUtils.h"
#include <TError.h>
#include <TH1D.h>
#include <TGraphAsymmErrors.h>
#include "PhysicsTools/RooStatsCms/interface/ClopperPearsonBinomialInterval.h"
#include "PhysicsTools/RooStatsCms/interface/FeldmanCousinsBinomialInterval.h"

ClassImp(mithep::MathUtils)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::AddInQuadrature(Double_t a, Double_t b)
{
  // Add quantities in quadrature.

  return(TMath::Sqrt(a*a + b*b));
}

//--------------------------------------------------------------------------------------------------
void MathUtils::CalcRatio(Double_t n1, Double_t n2, Double_t &r, Double_t &rlow, Double_t &rup, Int_t type = 0)
{
  // Calculate ratio and lower/upper errors from given values using various methods, dependent on type:
  //  0: Bayes
  //  1: Feldman-Cousins
  //  2: Clopper-Pearson
  
  if (n1>n2) {
    Error("CalcRatio", "First value should be smaller than second: %f > %f", n1, n2);
    r = n1/n2;
    rlow = 0;
    rup  = 0;
    return;
  }
 
  if (n2 >= 1) {
    if (type == 1) {
      //compute errors using Feldman-Cousins
      r = double(n1) / double(n2);
      FeldmanCousinsBinomialInterval fc;
      const double alpha = (1-0.682);
      fc.init(alpha);
      fc.calculate(n1, n2);
      rlow = r - fc.lower();
      rup = fc.upper() - r;    
      
    } else if (type == 2) {
      //compute errors using Clopper-Pearson
      r = double(n1) / double(n2);
      ClopperPearsonBinomialInterval cp;
      const double alpha = (1-0.682);
      cp.init(alpha);
      cp.calculate(n1, n2);
      rlow = r - cp.lower();
      rup = cp.upper() - r;
      
    } else {
      //compute using Bayes
      TH1D h1("dummy1","",1,1,2);
      h1.SetBinContent(1,n1);
      
      TH1D h2("dummy2","",1,1,2);
      h2.SetBinContent(1,n2);
      
      TGraphAsymmErrors g;
      g.BayesDivide(&h1,&h2);
      r = g.GetY()[0];
      rup = g.GetErrorYhigh(0);
      rlow = g.GetErrorYlow(0);
    } 
  } else {
    r = 0;
    rlow = 0;
    rup  = 0;
  }
}	


//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaPhi(Double_t phi1, Double_t phi2)
{
  // Compute DeltaPhi between two given angles. Results is in [-pi/2,pi/2].

  Double_t dphi = phi1-phi2;
  while (dphi > TMath::Pi())
    dphi -= TMath::TwoPi();
  while (dphi <= -TMath::Pi())
    dphi += TMath::TwoPi();
  
  
  return(dphi);
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaR(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2)
{
  // Compute DeltaR between two given points in the eta/phi plane.

  Double_t dR = TMath::Sqrt(DeltaR2(phi1,eta1,phi2,eta2));
  return(dR);
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaR2(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2)
{
  // Compute DeltaR between two given points in the eta/phi plane.

  Double_t dphi = DeltaPhi(phi1, phi2);
  Double_t deta = eta1-eta2;
  Double_t dR = dphi*dphi + deta*deta;

  return(dR);
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::Eta2Theta(Double_t eta) 
{ 
  // Compute theta from given eta value.

  return 2.*TMath::ATan(exp(-eta)); 
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::Theta2Eta(Double_t theta) 
{ 
  // Compute eta from given theta value.

  return -TMath::Log(TMath::Tan(theta/2.)); 
}
