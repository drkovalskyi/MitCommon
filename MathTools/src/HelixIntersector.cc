// $Id: HelixIntersector.cc,v 1.3 2009/03/20 13:33:19 loizides Exp $

#include "MitCommon/MathTools/interface/HelixIntersector.h"
#include <iostream>

ClassImp(mithep::HelixIntersector)

using namespace std;
using namespace mithep;

//--------------------------------------------------------------------------------------------------
HelixIntersector::HelixIntersector(const TVectorD *tr1, const TVector3 *momentum1,
				   const TVectorD *tr2, const TVector3 *momentum2) :
  fISec0(tr1,momentum1,tr2,momentum2),
  fISec1(tr1,momentum1,tr2,momentum2)
{
  // Constructor.

  // short cuts
  TrackParams &trk0 = fISec0.fTrk0;
  TrackParams &trk1 = fISec0.fTrk1;
  fISecs[0]         = &fISec0;
  fISecs[1]         = &fISec1;

  // Calculate intersections and test for existence

  // line connecting the two centers
  TVector3 connector =  trk1.Center() - trk0.Center();
  
  // intersection positions relative to connector
  double para  = (connector.Mag2()+trk0.Radius()*trk0.Radius()
		  -trk1.Radius()*trk1.Radius())/(2.0*connector.Mag());
  double perpsqr = trk0.Radius()*trk0.Radius() - para*para;
  
  // intersection exist if perp is real -> perp^2 is positive
  fHasIntersections = (perpsqr >= 0.0);
  
  // fill Intersection objects appropriately
  if (fHasIntersections) {

    // vector perpendicular to connector (sign doesn't matter)
    TVector3 perpvec(connector.y(),-connector.x(),0.0);
    perpvec.SetMag(sqrt(perpsqr));
    
    // make both combinations
    for (int i = 0; i < 2; i++) {
      // the intersection
      fISecs[i]->fLocation =
	trk0.Center() + para * connector.Unit() + double(2*i-1) * perpvec;
      
      // calculate tracks/intersection specific quantities for real intersections both tracks have
      // the same point
      fISecs[i]->SetLocation(fISecs[i]->fLocation,fISecs[i]->fLocation);
    }

    if (fabs(fISecs[0]->DeltaZ()) > fabs(fISecs[1]->DeltaZ())) {
      Intersection* tmp = fISecs[1];
      fISecs[1] = fISecs[0];
      fISecs[0] = tmp;
    }

  }
  else {
    // no intersection so we pick points of closest approach
    TVector3 locTrk0 = trk0.Center() + trk0.Radius() * connector.Unit();
    TVector3 locTrk1 = trk1.Center() - trk1.Radius() * connector.Unit();

    // mean is the "nominal" location
    fISecs[0]->fLocation = (locTrk0 + locTrk1) * 0.5;

    // for closest approach points actually on the tracks are used calculate tracks/intersection
    // specific quantities
    fISecs[0]->SetLocation(locTrk0,locTrk1);
  }  
}

//--------------------------------------------------------------------------------------------------
HelixIntersector::~HelixIntersector()
{
  // Destructor.
}

//--------------------------------------------------------------------------------------------------
HelixIntersector::TrackParams::TrackParams(const TVectorD *trk, const TVector3 *momentum) :
  Helix  ((*trk)(0),(*trk)(1),(*trk)(2),(*trk)(3),Angle((*trk)(4))),
  fTrkMom(*momentum)
{
  // Calculate center of r-phi plane circles.

  double rho = Radius() +  Helicity() * D0();
  double phi = Helicity() * M_PI/2.0 + Phi0();
  double z   = 0.0;
  if (0 && rho < 0) {
    printf("Cylindrical coordinates supplied with negative Rho\n");
    printf("Radius = %5f, Helicity = %5f, D0 = %5f, Rho=%5f\n", Radius(), Helicity(), D0(), rho);
  }
  fCenter.SetXYZ(rho*cos(phi),rho*sin(phi),z);

  // Original code ....
  //fCenter.SetRhoPhiZ(Radius() +  Helicity() * D0(),  
  //		       Helicity() * M_PI/2.0 + Phi0(),0.0);
}

//--------------------------------------------------------------------------------------------------
void HelixIntersector::TrackParams::SetLocation(TVector3& location)
{
  // Rotation about circle center.

  double deltaphi = M_PI - fCenter.Angle(location-fCenter);
  
  // check for backward locations
  if (fTrkMom.Dot(location)<0.0)
    deltaphi = - deltaphi;
    

  // calculate arc length to and z of intersection
  fArcLen = Radius()*deltaphi/SinTheta();
  fZAtISec = Z0() + Radius()*deltaphi*CotTheta();

  // calculate momentum of track at intersection
  fMomentum = Direction(fArcLen) * fTrkMom.Mag();
}

//--------------------------------------------------------------------------------------------------
HelixIntersector::Intersection::Intersection(const TVectorD *tr1, const TVector3 *momentum1,
					     const TVectorD *tr2, const TVector3 *momentum2) :
  fDeltaZ(0.0),
  fTrk0  (tr1,momentum1),
  fTrk1  (tr2,momentum2)
{
  // Constructor.

  fTrks[0] = &fTrk0;
  fTrks[1] = &fTrk1;
}

//--------------------------------------------------------------------------------------------------
void HelixIntersector::Intersection::SetLocation(TVector3& loc1, TVector3& loc2)
{
  // Calculate tracks specific quantities.

  fTrk0.SetLocation(loc1);
  fTrk1.SetLocation(loc2);
  
  // combined variables
  fDeltaZ = fTrk0.fZAtISec - fTrk1.fZAtISec;
  fLocation.SetZ((fTrk0.fZAtISec + fTrk1.fZAtISec)/2.0);
  fMomentum = fTrk0.fMomentum + fTrk1.fMomentum;
}
