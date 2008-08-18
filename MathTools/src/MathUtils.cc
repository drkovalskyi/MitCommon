// $Id $

#include "MitCommon/MathTools/interface/MathUtils.h"

using namespace mitMath;

//ClassImp(mitMath::MathUtils)

double mitMath::deltaPhi(double phi1, double phi2)
{
  double ans = fabs(phi1-phi2);
  while(ans>M_PI)
    ans = fabs(ans - 2*M_PI);
  return(ans);
}

double mitMath::phiEtaDeltaR(double phi1, double eta1, double phi2, double eta2)
{
  double dphi = fabs(phi1-phi2);
  while(dphi>M_PI) dphi = fabs(2*M_PI-dphi);
  double deta = fabs(eta1-eta2);
  double dR = TMath::Sqrt(dphi*dphi + deta*deta);
  return(dR);
}

double mitMath::deltaR(mithep::FourVector v1, mithep::FourVector v2)
{
  return mitMath::phiEtaDeltaR(v1.Phi(),v1.Eta(),v2.Phi(),v2.Eta());
}

double mitMath::addInQuadrature(double a, double b)
{
  return(TMath::Sqrt(a*a + b*b));
}

std::string mitMath::ftoa(double x)
{
  char a[100];
  sprintf(a,"%g",x);
  return(a);
}

double mitMath::Eta2Theta(double eta) { return 2.*atan(exp(-eta)) ; }
double mitMath::Theta2Eta(double theta) { return -log(tan(theta/2.)) ; }
