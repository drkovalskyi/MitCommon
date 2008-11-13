//--------------------------------------------------------------------------------------------------
// $Id: MultiVertexFitter.h,v 1.3 2008/11/13 16:34:28 paus Exp $
//
// MultiVertexFitter class header file
//
// Author: Craig Blocker, Brandeis CDF Group
//         Christoph Paus, MIT
//           + porting to CMS
//           + major cleanup
//           + removed track dependency
//           - removed CLHEP dependence, now using root for vector/matrix operations
//           - re-enbale the track extrapolation
//
// Description:  Does generic wrapping of the multiple vertex fitter (wrapper of CTVMFT).
//
// Revision History:
//   May 27, 1999  Craig Blocker     Creation
//   Jun 22, 2001  Craig Blocker     Fix a few bugs in checking on vertices (such as pointing to
//                                   primary vertex wasn't allowed)
//   May 07 2002   Mat Martin        Added support for "zero track" vertices: Pointing multi track
//                 Jonas Rademacker  vertices at another vertex which has no tracks directly	  
//                                   associated with it.  This is acheived through the routine:	  
//                                   vertexPoint_0track. See the implementation in		  
//                                   MultiVertexFitter.cc                                         
//   Sep 07, 2002 Joao Guimaraes     Added accessors for ijk errors and track-id of tracks that
//                                   produce some of those errors. These methods are:
//                                        getIJKErr(int&,int&,int&) const
//                                        getIJKErr0() const
//                                        getIJKErr1() const
//                                        getIJKErr2() const
//                                        getErrTrackId() const
//   Sep 15, 2002 Joao Guimaraes     Increased maximum number of tracks to 50 (_maxtrk).
//   Nov 11, 2002 Craig Blocker      Put protection in Fortran for potential bombs [divide by zero,
//                                   sqrt(neg. number]. Added particle specific addTrack routines
//                                   (addElectron, addPion, etc.) Added methods to directly access
//                                   the ctvmft error matrix and pointers.
//
//   Mar 18, 2003 Dmitry Litvintsev  added ability to restore the state of MultiVertexFitter
//                                   object from FORTRAN common blocks function is
//                                   'restoreFromCommons()' use it with care, fortran common blocks
//                                   have to be properly initialised. Contact me if unsure how to
//                                   use it. This routine is not needed for standard use of
//                                   MultiVertexFitter.
//   Oct 15, 2003 Sal Rappoccio      Added accessor to add tracks to ctvmft via a track id and the
//                                   track parameters and covariance
//   Dec 02, 2003 Sal Rappoccio      Added an interface to use extrapolated track errors in a
//                                   "transparent" way. The new interface is used as:
//
//                                     MultiVertexFitter fit;
//                                     fit.init();
//                                     fit.setPrimaryVertex(p1);    // Note: the order matters
//                                     fit.setReferencePoint(p2);
//                                     fit.addTrack(trk1, M_PION, 1);
//                                     fit.addTrack(trk2, M_PION, 1);
//                                     fit.fit();
//                                     Hep3Vector &pos = fit.getVertex(1);
//
//                                   "pos" will now hold the new vertex constructed from tracks
//                                   with errors referenced at the new reference point "p2". This
//                                   track extrapolation code is OFF by default. To turn it ON use:
//
//                                     fit.setReferencePoint(p2);
//
//                                   AddTrack was modified to implement this feature.
//   Jul 18, 2008  Christoph Paus    Works for CMS and any other type of Helix fitting.
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_VERTEXFIT_MULTIVERTEXFITTER_H
#define MITCOMMON_VERTEXFIT_MULTIVERTEXFITTER_H

#include <string>
#include <iostream>

#include <Rtypes.h>
#include <TVectorD.h>
#include <TMatrixDSym.h>
#include <CLHEP/Vector/LorentzVector.h>
#include <CLHEP/Vector/ThreeVector.h>
#include <CLHEP/Matrix/SymMatrix.h>
#include <CLHEP/Matrix/Vector.h>

#include "MitCommon/DataFormats/interface/Types.h"
#include "MitCommon/Ctvmft/interface/dimensions.hh"

//-------------------------------------------------------------------------------------------------
// Fortran routines to get address of the start of the ctvmq and ctvmfr common blocks
//-------------------------------------------------------------------------------------------------
extern "C" {
   int ctvmq_address_ (void);
   int ctvmfr_address_(void);
   int fiddle_address_(void);
   int trkprm_address_(void);
}

namespace mithep {
  class MultiVertexFitter {

    public:

      //--------------------------------------------------------------------------------------------
      // Enumerations
      //--------------------------------------------------------------------------------------------
      enum vertexNumber { PRIMARY_VERTEX,VERTEX_1 };
      enum vertexIndex  { X_INDEX=0, Y_INDEX, Z_INDEX, P1_INDEX, P2_INDEX };
      enum trackIndex   { CURVATURE_INDEX=0, PHI_INDEX, COTTH_INDEX };

      //--------------------------------------------------------------------------------------------
      // *structors
      //--------------------------------------------------------------------------------------------
      MultiVertexFitter();
      ~MultiVertexFitter() {}

      //--------------------------------------------------------------------------------------------
      // Fundamental funtions
      //--------------------------------------------------------------------------------------------
      void init                 (double bfield = 3.8);
      void setChisqMax          (const float chisqmx);
      
      //--------------------------------------------------------------------------------------------
      // CMS parameter ordering for the vector/matrix, which is assumed here:
      //
      //   qoverp, cotTheta, phi0, d0, z0;
      //   mapping to CDF is therefore { 1, 0*, 4, 3, 2 } 
      //
      // Note that the radius of curvature (in cm) is then:
      //
      //   Rc = cos(theta) / [ 0.0029979.... * (q/p) * B ],
      //
      // where B is the magnetic field in Tesla and tht is the angle between the field and the
      // direction. With p * cos(theta) = pT it follows:
      //
      //   Rc = pT / [ 0.0029979.... * q * B ],
      //   fullCurvature = 1 / Rc = 0.0029979 * q * B / pT = - 0.0029979 * B / pT.
      //
      //--------------------------------------------------------------------------------------------

      //--------------------------------------------------------------------------------------------
      // Parameter ordering for the vector/matrix, assumed here:
      //
      //   cotTheta, curvature, z0, d0, phi0;
      //
      // And also: curvature = -0.5 * CurvConst * BFieldT / pT (strictly speaking: half curvature)
      //                        CurvConst = 0.0029979, BField in Tesla
      //
      // Incidentally this is the ordering used in the CDF experiment.
      //--------------------------------------------------------------------------------------------
      bool addTrack             (const HepVector &pars, const HepSymMatrix &cov, int trackid,
                                 float mass, vertexNumber jv);

      bool addTrack             (const TVectorD &pars, const TMatrixDSym &cov, int trackid,
                                 double mass, vertexNumber jv);


      bool vertexPoint_2d       (vertexNumber jv1, vertexNumber jv2);
      bool vertexPoint_3d       (vertexNumber jv1, vertexNumber jv2);
      bool vertexPoint_1track   (vertexNumber jv1, vertexNumber jv2);
      bool vertexPoint_0track   (vertexNumber jv1, vertexNumber jv2);
      bool conversion_2d        (vertexNumber jv);
      bool conversion_3d        (vertexNumber jv);
      bool massConstrain        (int ntrk, const int trkIds[], float mass);

      void setPrimaryVertex     (float xv, float yv, float zv);
      void setPrimaryVertex     (Hep3Vector pv);
      bool setPrimaryVertexError(const HepSymMatrix &xverr);
      void setPrimaryVertexError(const float xverr[3][3]);

      bool fit                  ();

      void print                (std::ostream& os) const;
      void print                () const;
      void printErr             (std::ostream& os) const;
      void printErr             () const;

      void restoreFromCommons   ();

      // -------------------------------------------------------------------------------------------
      // Turn ON this option only when fitting the primary with no additional constraints. 
      // The routineis protected and will not function if those options are selected as well as the
      // beamline constraint.  In case this option is chosen the fit will calculate a result also 
      // with just one track in input.
      //
      // In order to enable this option the user must do the following: Fit the primary vertex as
      // VERTEX_1 (without any other vertices).  Enable beamlineConstraint by setting the pointing
      // constraint variable vtxpnt[0][0] to -100 (no pointing constraints when <0) Note that index
      // 0,0 corresponds to index 1,1 in fortran Provide the beamline parameters:
      //
      //  --> assume xv = xb + xzbslope*z, 
      //             yv = yb + yzbslope*z, where (xv, yv) is the transverse
      //      beam position at z, (xb, yb) is the transverse beam position at z = 0 
      //      and (xzbslope, yzbslope) are the slopes of the beamline
      //  --> set primary vertex at z=0 (xb, yb, 0)
      //  --> set the diagonal elements of the primary vertex error matrix to the
      //      sigma**2 of beam spread in x, y and z:
      //      xverr[1][1] =  (~30 - 50 um)**2
      //      xverr[2][2] =  (~30 - 50 um)**2
      //      xverr[3][3] =  (~30 cm) **2
      //      All other elements should be 0
      // For help email: Joao Guimaraes  guima@huhepl.harvard.edu, Franco Bedeschi bed@fnal.gov 
      // See more information in the ctvmft.f
      // -------------------------------------------------------------------------------------------
      bool beamlineConstraint(float xb, float yb, HepSymMatrix berr,float xzbslope,float yzbslope);
      bool beamlineConstraint(Hep3Vector pv, HepSymMatrix berr, float xzbslope, float yzbslope);

      //--------------------------------------------------------------------------------------------
      // Accessors
      //--------------------------------------------------------------------------------------------

      void             setBField  (double bField) { _bField = bField; }
      double           bField     () const { return _bField; }
      void             setExcuse  ();
      void             setNoExcuse();

      std::string      expert     () const;                // name/email of MultiVertexFitter expert
      int              status     () const;                // return status of fit
      // overall fit quality paramters
      int              ndof       () const;                // number of degrees of freedom
      float            prob       () const;                // return probability of chi-square
      float            chisq      () const;                // return chi-square of fit
      float            chisq      (const int trkId) const;
      float            chisq_rphi () const;
      float            chisq_rphi (const int trkId) const;
      float            chisq_z    () const;
      float            chisq_z    (const int trkId) const;
      float            chisq_rphiz() const;
      float            chisq_rphiz(const int trkId) const;

      // return fit track four momentum
      //HepLorentzVector getTrackP4 (const int trkId) const;
      FourVector       getTrackP4 (const int trkId) const;

      //// return fit track parameters
      //Helix            getHelix   (const int trkId) const;
      // return fit mass and get error
      float            getMass    (int ntrk, const int trkIds[], float& dmass) const;

      // return decay length
      float            getDecayLength(vertexNumber nv, vertexNumber mv, const Hep3Vector& dir,
                                      float& dlerr) const;
      float            getDecayLength(vertexNumber nv, vertexNumber mv, const ThreeVector& dir,
                                      float& dlerr) const;                                      
      float            getZDecayLength(vertexNumber nv, vertexNumber mv,
                                      const Hep3Vector& dir, float& dlerr) const;  
      float            getZDecayLength(vertexNumber nv, vertexNumber mv,
                                      const ThreeVector& dir, float& dlerr) const;                                        
      float            getImpactPar(vertexNumber prdV, vertexNumber dcyV,
                                      const Hep3Vector &v, float &dxyerr) const;     
      float            getImpactPar(vertexNumber prdV, vertexNumber dcyV,
                                      const ThreeVector &v, float &dxyerr) const;                      
      float            get_dr(vertexNumber nv, vertexNumber mv, float& drerr) const;
      float            get_dz(vertexNumber nv, vertexNumber mv, float& dzerr) const;
      // return location of vertex
      Hep3Vector       getVertexHep(vertexNumber nv) const;
      ThreeVector      getVertex(vertexNumber nv) const;
      // return error matrix element.
      ThreeSymMatrix   getErrorMatrix(vertexNumber nv) const;
      double           getErrorMatrixHep(int j, int k) const;
      HepSymMatrix     getErrorMatrixHep(vertexNumber nv) const;
      HepSymMatrix     getErrorMatrixHep(const int trkId) const;
      void             getPosMomErr  (HepMatrix& errors) const;
      int              vOff          (vertexNumber jv) const;
      int              tOff          (const int trkId) const;
      int              pOff          (int lp) const;
      int              cOff          (int lc) const;
      int              mOff          () const;
      double           VMat          (int i, int j) const;
      float            getPtError    (const int trkId) const;
      MultiVertexFitter::vertexNumber
      allocateVertexNumber();
      void             resetAllocatedVertexNumber();

      // Accessors for getting information relative to ijk errors get the error code from the three
      // ijk indexes into the argument variables
      void             getIJKErr(int& err0, int& err1, int& err2) const;
      // Return each error code from the three ijk indexes
      int              getIJKErr0() const;
      int              getIJKErr1() const;
      int              getIJKErr2() const;

      // Get track-id of the track causing a fatal error as indicated by the corresponding ijk error
      int              getErrTrackId() const;

      // Set new track reference point
      void             setTrackReferencePoint(const ThreeVector &ref);

      //--------------------------------------------------------------------------------------------
      // Overload operators
      //--------------------------------------------------------------------------------------------
      friend std::ostream& operator << (std::ostream& os, const MultiVertexFitter& vfit);

    protected:
      std::string      _expert;                 // string: name, email of MultiVertexFitter expert
      int              _stat;                   // status returned from fit.

      static const int _maxvtx = CTVMFT_MAXVTX; // maximum number of vertices
      static const int _maxmcn = CTVMFT_MAXMCN; // maximum number of mass constraints
      static const int _maxtrk = CTVMFT_MAXTRK; // maximum number of tracks
      static const int _maxitr = CTVMFT_MAXITR; // maximum number of iteration steps
      static const int _maxdim = 5*(CTVMFT_MAXVTX+1)+3*CTVMFT_MAXTRK+CTVMFT_MAXMCN;

      // FIDDLE must have access to MultiVertexFitter's (other) protected data like _maxvtx, etc. It
      // must therefore be a friend.
      struct        FIDDLE;
      friend struct FIDDLE;
      struct FIDDLE {
          int    excuse;
          float  chisqmax;
      };

      // CTVMQ must have access to MultiVertexFitter's (other) protected data like _maxvtx, etc. It
      // must therefore be a friend.
      struct        CTVMQ;
      friend struct CTVMQ;
      struct        CTVMQ {
          int    runnum;
          int    trgnum;
          int    iter;
          int    ntscut;
          int    nvertx;
          int    nmassc;
          int    ntrack;
          int    trkvtx[_maxvtx][_maxtrk];  // Logical in FORTRAN, but integer here to get the size
          int    trkmcn[_maxmcn][_maxtrk];  // Logical in FORTRAN, but integer here to get the size
          int    vtxpnt[2][_maxvtx];
          float  cmass [_maxmcn];
          int    cvtx  [_maxvtx];
          int    vtxvtx[_maxvtx][_maxvtx];  // Logical in FORTRAN, but integer here to get the size
          char   tkbank[_maxtrk][4];
          int    list  [_maxtrk];
          float  tmass [_maxtrk];
          int    matdim;
          int    tkerr [_maxtrk];
          int    ndof;
          float  chisqr[_maxitr+1];
          float  chit  [_maxtrk];
          float  chiv  [_maxvtx+1];
          float  chim  [_maxmcn];
          float  xyzpv0[3];
          float  exyzpv[3][3];
          float  xzslope;
          float  yzslope;
          float  xyzvrt[_maxvtx+1][3];
          float  dxyzpv[3];
          float  par   [_maxtrk][5];
          float  g     [_maxtrk][5][5];
          float  trkp4 [6][_maxtrk];
          float  vtxp4 [_maxvtx][4];
          float  mcnp4 [_maxmcn][4];
          float  dda   [8][_maxtrk];
          int    voff  [_maxvtx];
          int    toff  [_maxtrk];
          int    poff  [_maxvtx];
          int    coff  [_maxvtx];
          int    moff;
          float  par0  [_maxtrk][5];
          float  pardif[_maxtrk][5];
          float  fmcdif[_maxmcn];
          float  pcon  [2][_maxvtx];
          float  sang  [2][_maxvtx];
          float  drmax;
          float  rvmax;
          float  dzmax;
          float  trnmax;
          float  dsmin;
          int    ijkerr[3];
          float  pscale;
      };

      // CTVMFR must have access to MultiVertexFitter's (other) protected data like _maxdim, etc. 
      // It must therefore be a friend.
      struct        CTVMFR;
      friend struct CTVMFR;
      struct        CTVMFR {
          double vmat[_maxdim+1][_maxdim];
      };

      // TRKPRM must have access to MultiVertexFitter's (other) protected data like _maxtrk, etc.
      // It must therefore be a friend.
      struct        TRKPRM;
      friend struct TRKPRM;
      struct        TRKPRM {
          float  trhelix[_maxtrk][5];
          float  trem   [_maxtrk][5][5];
      };

      CTVMQ    _ctvmq;
      CTVMQ*   _ctvmq_com;
      CTVMFR   _ctvmfr;
      CTVMFR*  _ctvmfr_com;
      FIDDLE   _fiddle;
      FIDDLE*  _fiddle_com;
      TRKPRM   _trkprm;
      TRKPRM*  _trkprm_com;

      //--------------------------------------------------------------------------------------------
      // Private functions used by class
      //--------------------------------------------------------------------------------------------

    private:

      // Moves reference point of track parameters and errors to _referencePoint
      //void moveReferencePoint(HepVector &v, HepSymMatrix &m);
      // Moves reference point of track parameters and errors to _referencePoint.  Additionally 
      //checks if track has already been moved by examining it's "derived" link
      //void moveReferencePoint(const CdfTrack *trk, HepVector &v, HepSymMatrix &m);

      //--------------------------------------------------------------------------------------------
      // Data members of class
      //--------------------------------------------------------------------------------------------
      double      _bField;                        // B field in Tesla

      int         _currentAllocatedVertexNumber;  // index to enum vertexNumber
      ThreeVector _referencePoint;                // reference point of track
      Hep3Vector  _primaryVertex;                 // primary vertex relative to _referencePoint
      Hep3Vector  _cdfPrimaryVertex;              // primary vertex in CDF coordinate system
      bool        _extrapolateTrackErrors;        // extrapolate track errors to _referencePoint
  };
}

//--------------------------------------------------------------------------------------------------
// Inline functions
//--------------------------------------------------------------------------------------------------
inline
std::ostream& operator << (std::ostream& os, const mithep::MultiVertexFitter& vfit)
{
  vfit.print(os);
  return    (os);
}
inline
void mithep::MultiVertexFitter::setExcuse()
{
  _fiddle.excuse = 1;   // the default
}
inline
void mithep::MultiVertexFitter::setNoExcuse()
{
  _fiddle.excuse = 0;   // crash on input error
}
inline
void mithep::MultiVertexFitter::setChisqMax(const float chisqmx)
{
  _fiddle.chisqmax = chisqmx;
}
#endif
