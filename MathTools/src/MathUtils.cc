// $Id:$

#include "MitCommon/MathTools/interface/MathUtils.h"

using namespace mitcommon;

//--------------------------------------------------------------------------------------------------
double MathUtils::AddInQuadrature(double a, double b)
{
  return(TMath::Sqrt(a*a + b*b));
}

//--------------------------------------------------------------------------------------------------
double MathUtils::DeltaPhi(double phi1, double phi2)
{
  double ans = fabs(phi1-phi2);
  while(ans>M_PI)
    ans = fabs(ans - 2*M_PI);
  return(ans);
}

//--------------------------------------------------------------------------------------------------
double MathUtils::DeltaR(double phi1, double eta1, double phi2, double eta2)
{
  double dphi = fabs(phi1-phi2);
  while(dphi>M_PI) dphi = fabs(2*M_PI-dphi);
  double deta = fabs(eta1-eta2);
  double dR = TMath::Sqrt(dphi*dphi + deta*deta);
  return(dR);
}

//--------------------------------------------------------------------------------------------------
double MathUtils::DeltaR(const FourVector &v1, const FourVector &v2)
{
  return MathUtils::DeltaR(v1.Phi(),v1.Eta(),v2.Phi(),v2.Eta());
}

//--------------------------------------------------------------------------------------------------
double MathUtils::Eta2Theta(double eta) 
{ 
  return 2.*TMath::ATan(exp(-eta)); 
}

//--------------------------------------------------------------------------------------------------
double MathUtils::Theta2Eta(double theta) 
{ 
  return -TMath::Log(TMath::Tan(theta/2.)); 
}
