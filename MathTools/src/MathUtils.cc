// $Id: MathUtils.cc,v 1.3 2008/09/10 03:27:35 loizides Exp $

#include "MitCommon/MathTools/interface/MathUtils.h"

using namespace mithep;

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::AddInQuadrature(Double_t a, Double_t b)
{
  return(TMath::Sqrt(a*a + b*b));
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaPhi(Double_t phi1, Double_t phi2)
{
  Double_t dphi = TMath::Abs(phi1-phi2);
  while (dphi>TMath::Pi())
    dphi = TMath::Abs(dphi - TMath::TwoPi());
  return(dphi);
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaPhi(const FourVector &v1, const FourVector &v2)
{
  return DeltaPhi(v1.Phi(),v2.Phi());
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaR(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2)
{
  Double_t dphi = DeltaPhi(phi1, phi2);
  Double_t deta = eta1-eta2;
  Double_t dR = TMath::Sqrt(dphi*dphi + deta*deta);
  return(dR);
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaR(const FourVector &v1, const FourVector &v2)
{
  return MathUtils::DeltaR(v1.Phi(),v1.Eta(),v2.Phi(),v2.Eta());
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::Eta2Theta(Double_t eta) 
{ 
  return 2.*TMath::ATan(exp(-eta)); 
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::Theta2Eta(Double_t theta) 
{ 
  return -TMath::Log(TMath::Tan(theta/2.)); 
}
