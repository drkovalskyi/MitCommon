#include "MitCommon/MathTools/interface/Helix.h"
#include <TSystem.h>

ClassImp(mithep::Helix)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
Helix::Helix(const TVector3 &MomentumGev, const TVector3 &PositionCm,
	     double q, double BFieldTesla)
{
  if (BFieldTesla != 0.0 && q != 0.0) {
    double CurvatureConstant = 0.0029979;
    double Helicity          = -1.0*fabs(BFieldTesla)*fabs(q)/(BFieldTesla*q);
    double Radius            = fabs(MomentumGev.Perp()/(CurvatureConstant*BFieldTesla*q));
    fCurvature               = 1.0/(2.0*MomentumGev.Perp()/(CurvatureConstant*BFieldTesla*q));
    
    double W = HUGE_VAL;
    if (Radius != 0.0)
      W = Helicity/Radius;

    Angle  phi1     = MomentumGev.Phi();
    double x        = PositionCm.x(),
           y        = PositionCm.y(),
           z        = PositionCm.z(); 
    double sinPhi1  = sin(phi1),
           cosPhi1  = cos(phi1);
    double gamma    = atan((x*cosPhi1 + y*sinPhi1)/(x*sinPhi1-y*cosPhi1 -1/W));
    fPhi0            = phi1+gamma;
    fD0              = ((1/W + y*cosPhi1 - x*sinPhi1) /cos(gamma) - 1/W);
    fCotTheta        = MomentumGev.z()/MomentumGev.Perp();
    fZ0              = z + gamma*fCotTheta/W;
  }
  // For the special case that we have no Helix, never happens for us .. right?
  else {
    std::cout << "MAJOR ERROR in HELIX!!!! Should not happen. STOP!!!!\n";
    gSystem->Exit(123);
  }
}

//--------------------------------------------------------------------------------------------------
bool Helix::operator == (const Helix &right) const
{
  return
    fCotTheta  == right.fCotTheta  &&
    fCurvature == right.fCurvature &&
    fZ0        == right.fZ0        &&
    fD0        == right.fD0        &&
    fPhi0      == right.fPhi0;
}

//--------------------------------------------------------------------------------------------------
bool Helix::operator != (const Helix & right) const
{
  return !((*this)==right);
}

//--------------------------------------------------------------------------------------------------
Helix::~Helix()
{
  delete fVParameters;
}

//--------------------------------------------------------------------------------------------------
TVector3 Helix::Position(double s) const
{
  fCacheSinesAndCosines(s);
  if (s == 0.0 || fCurvature == 0.0) {
    return TVector3(-fD0*fSinPhi0+s*fCosPhi0*fSinTheta,
		     fD0*fCosPhi0+s*fSinPhi0*fSinTheta,
		     fZ0+s*fCosTheta);
  }
  else {
    return TVector3((fCosPhi0*fSs-fSinPhi0*(2.0*fCurvature*fD0+1.0-fCc))/(2.0*fCurvature),
		    (fSinPhi0*fSs+fCosPhi0*(2.0*fCurvature*fD0+1.0-fCc))/(2.0*fCurvature),   
		    fS*fCosTheta + fZ0);
  }
}

//--------------------------------------------------------------------------------------------------
TVector3 Helix::Direction(double s) const
{
  fCacheSinesAndCosines(s);
  if (s == 0.0) {
    return TVector3(fCosPhi0*fSinTheta,fSinPhi0*fSinTheta,fCosTheta);
  }
  double xtan = fSinTheta*(fCosPhi0*fCc -fSinPhi0*fSs); 
  double ytan = fSinTheta*(fCosPhi0*fSs +fSinPhi0*fCc); 
  double ztan = fCosTheta;

  return TVector3(xtan,ytan,ztan);
}

//--------------------------------------------------------------------------------------------------
double Helix::PathLengthAtRhoEquals(double rho) const
{
  return (SinTheta()?(L2DAtR(rho)/SinTheta()):0.0);
}

//--------------------------------------------------------------------------------------------------
double Helix::InverseRadius() const
{
  return fCurvature*2.0;
}

//--------------------------------------------------------------------------------------------------
double Helix::Radius() const
{
  return fabs(1.0/InverseRadius());
}

//--------------------------------------------------------------------------------------------------
SignedAngle Helix::TurningAngle(double s) const
{
  return s/Radius();
}

//--------------------------------------------------------------------------------------------------
double Helix::Curvature() const
{
  return fCurvature;
}

//--------------------------------------------------------------------------------------------------
double Helix::Helicity() const
{
  return fCurvature>0 ? 1.0 : -1.0 ;
}

//--------------------------------------------------------------------------------------------------
double Helix::CotTheta() const
{
  return fCotTheta;
}

//--------------------------------------------------------------------------------------------------
Angle Helix::Phi0() const
{
  return fPhi0;
}

//--------------------------------------------------------------------------------------------------
double Helix::D0() const
{
  return fD0;
}

//--------------------------------------------------------------------------------------------------
double Helix::SignLz() const
{
  return (fD0>0) ? -1.0 : 1.0;
}

//--------------------------------------------------------------------------------------------------
double Helix::Z0() const
{
  return fZ0;
}

//--------------------------------------------------------------------------------------------------
TVector3 Helix::SecondDerivative(double s) const
{
  double phi1    = fPhi0+s*2.0*fCurvature*fSinTheta;
  double sp1     = sin(phi1);
  double xsecond = -fSinTheta*sp1*2.0*fCurvature*fSinTheta;
  double ysecond = fSinTheta*sqrt(1.0-sp1*sp1)*2.0*fCurvature*fSinTheta;

  return TVector3(xsecond,ysecond,0.0);
}

//--------------------------------------------------------------------------------------------------
double Helix::SinPhi0() const
{
  fRefreshCache();

  return fSinPhi0;
}

//--------------------------------------------------------------------------------------------------
double Helix::CosPhi0() const
{
  fRefreshCache();

  return fCosPhi0;
}

//--------------------------------------------------------------------------------------------------
double Helix::SinTheta() const
{
  fRefreshCache();

  return fSinTheta;
}

//--------------------------------------------------------------------------------------------------
double Helix::CosTheta() const
{
  fRefreshCache();

  return fCosTheta;
}

//--------------------------------------------------------------------------------------------------
Angle Helix::PhiAtR(double rho) const
{
  double c = Curvature();
  double d = D0();
  double a = c/(1.0+2.0*c*d);
  double b = d - a*d*d;
  double arcsin = a*rho+b/rho;
  if (arcsin>1.0 || arcsin<-1.0) {
    return (arcsin > 0.) ? M_PI : -M_PI;
  }
  Angle phi = Phi0() + asin(arcsin);

  return phi;
}

//--------------------------------------------------------------------------------------------------
double Helix::ZAtR(double rho) const
{
  return fZ0 + CotTheta()*L2DAtR(rho);
}

//--------------------------------------------------------------------------------------------------
double Helix::L2DAtR(double rho) const
{
  double L2D;

  double c = Curvature();
  double d = D0();

  if (c!=0.0) {
    double rad = (rho*rho-d*d)/(1.0+2.0*c*d);
    double rprime;

    if (rad<0.0)
      rprime = 0.0;
    else
      rprime = sqrt( rad );

    if (c*rprime > 1.0 || c*rprime < -1.0)
      L2D = c*rprime > 0. ? M_PI/c : -M_PI/c;
    else
      L2D = asin(c*rprime)/c;

  }
  else {
    L2D = rho;
  }

  return L2D;
}

//--------------------------------------------------------------------------------------------------
double Helix::CosAlphaAtR(double rho) const
{
  double c   = Curvature();
  double d   = D0();
  double a   = c/(1.0+2.0*c*d);
  double phi = PhiAtR(rho);

  double x   = rho*cos(phi);
  double y   = rho*sin(phi);

  double dir_x0 = (1.0+2.0*c*d) * (1.0-2.0*a*y);
  double dir_y0 = (1.0+2.0*c*d) * 2.0*a*x;

  double phi0 = Phi0();
  double dir_x = dir_x0*cos(phi0) - dir_y0*sin(phi0);
  double dir_y = dir_x0*sin(phi0) + dir_y0*cos(phi0);

  double cosalpha = (x*dir_x + y*dir_y)/rho;

  return cosalpha;
}
