//--------------------------------------------------------------------------------------------------
// $Id: Helix.h,v 1.4 2008/09/14 15:28:19 loizides Exp $
//
// Class Helix
//
// Implementation nof a general helix class with a set of tools for working with objects like 
// tracks and finding intersections (vertices).
//
// Author: C.Paus, stolen and adjusted from CDF implementation of Kurt Rinnert
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_MATHTOOLS_HELIX_H
#define MITCOMMON_MATHTOOLS_HELIX_H

#include <TVector.h>
#include <TVector3.h>
#include "MitCommon/MathTools/interface/Angle.h"

namespace mithep {
  class Helix {

    public:

      //  Constructor
      inline Helix();

      // Construct a helix in various ways
      Helix(const TVector3 &mom, const TVector3 &pos, double q, double Field);
      inline Helix(const Helix &right);
      inline Helix(double cotTheta, double curvature, double z0, double d0, Angle phi0);
      static Helix create(const TVector & v);

      // Destructor
      virtual ~Helix();

      // Operators
      inline const          Helix & operator=(const Helix &right);
      bool                  operator == (const Helix & right) const;
      bool                  operator != (const Helix & right) const;

      inline void           SetParameters(const TVector &p);
      inline const TVector &Parameters() const;
      inline void           SetCotTheta(double cotTheta);
      inline void           SetCurvature(double curvature);
      inline void           SetZ0(double z0);
      inline void           SetD0(double d0);
      inline void           SetPhi0(Angle phi0);

      virtual TVector3      Position(double s = 0.0) const;
      virtual TVector3      Direction(double s = 0.0) const;
      //  the second derivative of the helix vs (three-dimensional) path length
      virtual TVector3      SecondDerivative(double s = 0.0) const;
      // Get both position and direction at once.
      //virtual void Location(Trajectory::Location &loc, double s = 0.0) const;
      // Get pathlength at fixed rho=sqrt(x^2 + y^2)
      virtual double        PathLengthAtRhoEquals(double rho) const;

      // Get certain parameters as a function of two-dimensional R.
      Angle                 PhiAtR(double r) const;
      double                ZAtR(double r) const;
      double                L2DAtR(double r) const;
      double                CosAlphaAtR(double r) const;
      double                InverseRadius() const;
      double                Radius() const;
      SignedAngle           TurningAngle(double s) const; // turning angle as function of path length
      double                Curvature() const;
      double                Helicity() const;             // helicity, positive = counterclockwise
      double                CotTheta() const;
      Angle                 Phi0() const;
      double                D0() const;
      double                Z0() const;
      double                SignLz() const;
      double                SinPhi0() const;
      double                CosPhi0() const;
      double                SinTheta() const;
      double                CosTheta() const;
      static unsigned int   ParameterSpaceSize() { return 5; }

      //============================================================================================
      // KCDF: analytical computation of helix/plane intersection.
      //
      // What we really compute is the intersection of a line and a circle (projected helix) 
      // in the x-y plane.
      //
      //     >>>>>>>>>>  W A R N I N G    W A R N I N G    W A R N I N G  <<<<<<<<<<
      //     >                                                                     <
      //     > We assume the plane to be parallel or perpendicular                 <
      //     > to the z axis (i.e. the B-field),                                   <
      //     > since otherwise there is no analytical solution. (One would end up  <
      //     > with an equation of type cos(alpha) == alpha.)                      <
      //     > Although we know this assumption doesn´t hold exactly, we think it  <
      //     > is a reasonable first approximation.                                <
      //     > In cases when our assumption has to be dropped, one can use the     <
      //     > intersection point computed here as a *good* starting point of a    <
      //     > numerical method, or one uses another analytical method to compute  <
      //     > the intersection of the tangent line at the point and the plane.    <
      //     > We plan to use one of these approaches in the near future, but      <
      //     > this is NOT YET IMPLEMENTED!                                        <
      //     > For the time being, we invoke the old numerical                     <
      //     > Trajectory::newIntersectionWith in such circumstances.              <
      //     >                                                                     <
      //     >>>>>>>>>>  W A R N I N G    W A R N I N G    W A R N I N G  <<<<<<<<<<
      //
      // Kurt Rinnert,  08/31/1998
      //============================================================================================
      //Location* newIntersectionWith(const HepPlane3D& plane) const;

    private:
      // This is the Helix
      double           fCotTheta;
      double           fCurvature;
      double           fZ0;
      double           fD0;
      Angle            fPhi0;

      // This is the cache
      mutable TVector *fVParameters;
      mutable bool     fIsStale;
      mutable double   fSinPhi0;
      mutable double   fCosPhi0;
      mutable double   fSinTheta;
      mutable double   fCosTheta;
      mutable double   fS;
      mutable double   fAa;
      mutable double   fSs;
      mutable double   fCc;
      mutable bool     fCenterIsValid; //needed by newIntersectionWith KR
      mutable double   fMx;
      mutable double   fMy;

      // Needed whenever fSinPhi0, fCosPh0, fSinTheta, or fCosTheta is used.
      inline void      fRefreshCache() const;
  
      // Needed whenever fSs or fCc are used.
      inline void      fCacheSinesAndCosines(double s) const;
  
  };
}

// Update fSinTheta,fCosTheta,fSinPhi0, and fCosPhi0
inline
void mithep::Helix::fRefreshCache() const {
  if (fIsStale) {
    fIsStale=false;
    double theta;
    if ( fCotTheta==0.0 ) {
      theta = M_PI/2.0;          
    }
    else {
      theta=atan(1.0/fCotTheta);
      if (theta<0.0) theta+=M_PI;
    }
    if (theta==0.0) {
      fSinTheta=0.0;
      fCosTheta=1.0;
    }
    else {
      fCosTheta=cos(theta);
      fSinTheta=sqrt(1-fCosTheta*fCosTheta);
    }
    if (fPhi0==0.0) {
      fSinPhi0=0.0;
      fCosPhi0=1.0;
    }
    else {
      fCosPhi0 = cos(fPhi0);
      fSinPhi0 = sqrt(1.0-fCosPhi0*fCosPhi0);
      if (fPhi0>M_PI) fSinPhi0 = -fSinPhi0;
    }
  }
}

// Update fS, fAa, fSs, and fCc if the arclength has changed.
inline
void mithep::Helix::fCacheSinesAndCosines(double s) const {
  fRefreshCache();
  if (fS!=s){
    fS=s;
    fAa=2.0*fS*fCurvature*fSinTheta;
    if (fAa==0.0) {
      fSs=0.0;
      fCc=1.0;
    }
    else {
      fSs=sin(fAa);
      fCc=cos(fAa);
    }
  }
}


inline
mithep::Helix::Helix() :
  fCotTheta(0.0),
  fCurvature(0.0),
  fZ0(0.0),
  fD0(0.0),
  fPhi0(0.0),
  fVParameters(0),
  fIsStale(1),
  fSinPhi0(2),
  fCosPhi0(2),
  fSinTheta(2),
  fCosTheta(2),
  fS(-999.999),
  fAa(2),
  fSs(2),
  fCc(2),
  fCenterIsValid(false),
  fMx(0.0),
  fMy(0.0)
{
}

inline
mithep::Helix::Helix(const Helix &right) :
  fCotTheta(right.fCotTheta),
  fCurvature(right.fCurvature),
  fZ0(right.fZ0),
  fD0(right.fD0),
  fPhi0(right.fPhi0),
  fVParameters(0),
  fIsStale(right.fIsStale),
  fSinPhi0(right.fSinPhi0),
  fCosPhi0(right.fCosPhi0),
  fSinTheta(right.fSinTheta),
  fCosTheta(right.fCosTheta),
  fS(right.fS),
  fAa(right.fAa),
  fSs(right.fSs),
  fCc(right.fCc),
  fCenterIsValid(right.fCenterIsValid),
  fMx(right.fMx),
  fMy(right.fMy)
{
  if (right.fVParameters)
    fVParameters=new TVector(*(right.fVParameters));
}

inline
mithep::Helix::Helix(double cotTheta,double curvature,double z0,double d0, Angle  phi0) :
  fCotTheta(cotTheta),
  fCurvature(curvature),
  fZ0(z0),
  fD0(d0),
  fPhi0(phi0),
  fVParameters(0),
  fIsStale(1),
  fSinPhi0(2),
  fCosPhi0(2),
  fSinTheta(2),
  fCosTheta(2),
  fS(-999.999),
  fAa(2),
  fSs(2),
  fCc(2),
  fCenterIsValid(false),
  fMx(0.0),
  fMy(0.0)
{
}

// Assign helix and cache
inline
const mithep::Helix & mithep::Helix::operator=(const mithep::Helix &right)
{
  if (this != &right) {
    fCotTheta=right.fCotTheta;
    fCurvature=right.fCurvature;
    fZ0=right.fZ0;
    fD0=right.fD0;
    fPhi0=right.fPhi0;
    fIsStale= right.fIsStale;
    fSinTheta=right.fSinTheta;
    fCosTheta=right.fCosTheta;
    fSinPhi0=right.fSinPhi0;
    fCosPhi0=right.fCosPhi0;
    fS=right.fS;
    fAa=right.fAa;
    fSs=right.fSs;
    fCc=right.fCc;
    if (fVParameters)
      delete fVParameters;
    fVParameters=0;
    if (right.fVParameters)
      fVParameters = new TVector(*(right.fVParameters));
    fCenterIsValid = right.fCenterIsValid;
    fMx = right.fMx;
    fMy = right.fMy;
  }
  return *this;
}

inline
void mithep::Helix::SetCotTheta(double cotTheta)
{
  fCotTheta=cotTheta;
  fIsStale=true;
}

inline
void mithep::Helix::SetZ0(double z0)
{
  fZ0 = z0;
  fIsStale=true;
}

inline
void mithep::Helix::SetCurvature(double curvature)
{
  fCurvature=curvature;
  fIsStale=true;
}

inline
void mithep::Helix::SetD0(double d0)
{
  fD0=d0;
  fIsStale=true;
}

inline
void mithep::Helix::SetPhi0(Angle phi0)
{
  fPhi0=phi0;
  fIsStale=true;
}

inline
void mithep::Helix::SetParameters(const TVector &p)
{
  // Check we're getting a sensible vector
  if (p.GetNrows() < 5) {
    return;
  }
  SetCotTheta(p[0]);
  SetCurvature(p[1]);
  SetZ0(p[2]);
  SetD0(p[3]);
  SetPhi0(p[4]);
  fIsStale = true;
}

inline
const TVector & mithep::Helix::Parameters() const
{
  if (!fVParameters)
    fVParameters = new TVector(5);
  (*fVParameters)(1)=CotTheta(); 
  (*fVParameters)(2)=Curvature();
  (*fVParameters)(3)=Z0();
  (*fVParameters)(4)=D0();
  (*fVParameters)(5)=Phi0();
  return *fVParameters;;
}
 
inline
mithep::Helix mithep::Helix::create(const TVector & v)  {
  return mithep::Helix(v(1),v(2),v(3),v(4),v(5));
}
#endif
