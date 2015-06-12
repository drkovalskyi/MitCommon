//--------------------------------------------------------------------------------------------------
// Class HelixIntersector
//
// Finds the intersection of two tracks. 
// If they do not intersect it finds the point of closest approach.
//
// Author List: C.Paus (stolen from CDF implementation of E. Lipeles,
//                      therefore not all our coding conventions fulfilled) 
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_MATHTOOLS_HELIXINTERSECTOR_H
#define MITCOMMON_MATHTOOLS_HELIXINTERSECTOR_H

#include <iostream>
#include "MitCommon/MathTools/interface/Helix.h"

namespace mithep 
{
  class HelixIntersector {

    public:
      //--------------------------------------------------------------------------------------------
      // Class for track properties at intersection
      //--------------------------------------------------------------------------------------------
      class Intersection;
      class TrackParams : public Helix
      {
          friend class HelixIntersector;
          friend class HelixIntersector::Intersection;

        public:
          // Sets the track and calculates center
          TrackParams(const TVectorD *params, const TVector3 *momentum);

          const  TVector3 &Momentum              () const { return fMomentum; }
          const  TVector3 &Center                () const { return fCenter;   }
          double           ArcLengthToInterection() const { return fArcLen;   }
          double           ZAtIntersection       () const { return fZAtISec;  }

        private:
          // Make it illegal to copy or assign
          TrackParams(const TrackParams &);
          TrackParams & operator = (const TrackParams &);
    
          // Calculates the location dependent quantities at the specified location
          void             SetLocation           (TVector3 &location);

          // Data
          TVector3         fMomentum;
          TVector3         fCenter;
          double           fArcLen;
          double           fZAtISec;
          TVector3         fTrkMom;
      };  

      //--------------------------------------------------------------------------------------------
      // Class for Intersection results
      //--------------------------------------------------------------------------------------------
      class Intersection
      {
          friend class HelixIntersector;

        public:
          // Properties of the vertex
          const TVector3    &Location() const { return fLocation; }
          double             DeltaZ  () const { return fDeltaZ; }
    
          // Combined momenta of the two tracks
          const TVector3    &Momentum() const { return fMomentum; }
    
          // Properties of the vertex daughters
          // i = 0 or 1 for the two daughters
          const TrackParams &TrackParamsI(int i) const { return *(fTrks[i]); }
    
        private:
          Intersection(const TVectorD *tr1, const TVector3 *momentum1,
                       const TVectorD *tr2, const TVector3 *momentum2);
          virtual ~Intersection() {}

          void               SetLocation(TVector3 &loc1, TVector3 &loc2);

          double             fDeltaZ;
          TVector3           fLocation;
          TVector3           fMomentum;
          TrackParams       *fTrks[2];
          TrackParams        fTrk0;
          TrackParams        fTrk1;
      };

      //--------------------------------------------------------------------------------------------
      // HelixIntersector methods
      //--------------------------------------------------------------------------------------------
      // Constructors, does intersection work
      HelixIntersector(const TVectorD *tr1, const TVector3 *momentum1,
                       const TVectorD *tr2, const TVector3 *momentum2);

      // Destructor
      virtual ~HelixIntersector();

      // Intersections or points of closest approach
      //   i = 0 or 1 (1 only if there are real intersections)
      //   ordered by deltaZ (so in general you only want i=0);
      const Intersection& IntersectionI(int i) const { return *(fISecs[i]); }

      // Has intersections
      bool                HasIntersections() const { return fHasIntersections; }

    private:
      // Make it illegal to copy or assign
      HelixIntersector(const HelixIntersector &);
      HelixIntersector   &operator = (const HelixIntersector &);
  
      bool                fHasIntersections;
      Intersection       *fISecs[2];
      Intersection        fISec0;
      Intersection        fISec1;
  };
}
#endif
