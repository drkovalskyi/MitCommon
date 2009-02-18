// $Id: MathUtils.cc,v 1.4 2008/11/24 14:23:01 loizides Exp $

#include "MitCommon/MathTools/interface/MathUtils.h"

using namespace mithep;

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::AddInQuadrature(Double_t a, Double_t b)
{
  // Add quantities in quadrature.

  return(TMath::Sqrt(a*a + b*b));
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaPhi(Double_t phi1, Double_t phi2)
{
  // Compute DeltaPhi between two given angles. Results is in [-pi/2,pi/2].

  Double_t dphi = TMath::Abs(phi1-phi2);
  while (dphi>TMath::Pi())
    dphi = TMath::Abs(dphi - TMath::TwoPi());
  return(dphi);
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaPhi(const FourVector &v1, const FourVector &v2)
{
  // Compute DeltaPhi between two given angles. Results is in [-pi/2,pi/2].

  return DeltaPhi(v1.Phi(),v2.Phi());
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaPhi(const FourVectorM &v1, const FourVectorM &v2)
{
  // Compute DeltaPhi between two given angles. Results is in [-pi/2,pi/2].

  return DeltaPhi(v1.Phi(),v2.Phi());
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaR(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2)
{
  // Compute DeltaR between two given points in the eta/phi plane.

  Double_t dphi = DeltaPhi(phi1, phi2);
  Double_t deta = eta1-eta2;
  Double_t dR = TMath::Sqrt(dphi*dphi + deta*deta);
  return(dR);
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaR(const FourVector &v1, const FourVector &v2)
{
  // Compute DeltaR between two given points in the eta/phi plane.

  return MathUtils::DeltaR(v1.Phi(),v1.Eta(),v2.Phi(),v2.Eta());
}

//--------------------------------------------------------------------------------------------------
Double_t MathUtils::DeltaR(const FourVectorM &v1, const FourVectorM &v2)
{
  // Compute DeltaR between two given points in the eta/phi plane.

  return MathUtils::DeltaR(v1.Phi(),v1.Eta(),v2.Phi(),v2.Eta());
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
