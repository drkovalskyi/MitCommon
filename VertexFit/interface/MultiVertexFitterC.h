//--------------------------------------------------------------------------------------------------
// $Id: MultiVertexFitterC.h,v 1.2 2009/03/20 13:33:03 loizides Exp $
//
// MultiVertexFitterC class header file
//
// Compact version of MultiVertexFitter
//
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_VERTEXFIT_MULTIVERTEXFITTERC_H
#define MITCOMMON_VERTEXFIT_MULTIVERTEXFITTERC_H

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
#include "MitCommon/Ctvmft/interface/cdimensions.hh"

//-------------------------------------------------------------------------------------------------
// Fortran routines to get address of the start of the ctvmq and ctvmfr common blocks
//-------------------------------------------------------------------------------------------------
/* extern "C" { */
/*    int cctvmq_address_ (void); */
/*    int cctvmfr_address_(void); */
/*    int cfiddle_address_(void); */
/*    int ctrkprm_address_(void); */
/* } */

namespace mithep {
  class MultiVertexFitterC {

    public:
      //--------------------------------------------------------------------------------------------
      // Enumerations
      //--------------------------------------------------------------------------------------------
      enum vertexNumber { PRIMARY_VERTEX,VERTEX_1,VERTEX_2,VERTEX_3,VERTEX_4,VERTEX_5,VERTEX_6 };
      enum vertexIndex  { X_INDEX=0, Y_INDEX, Z_INDEX, P1_INDEX, P2_INDEX };
      enum trackIndex   { CURVATURE_INDEX=0, PHI_INDEX, COTTH_INDEX };

      //--------------------------------------------------------------------------------------------
      // *structors
      //--------------------------------------------------------------------------------------------
      MultiVertexFitterC();
      ~MultiVertexFitterC() {}

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
      bool addTrack             (const CLHEP::HepVector &pars, const CLHEP::HepSymMatrix &cov, 
                                 int trackid, float mass, vertexNumber jv);

      bool addTrack             (const TVectorD &pars, const TMatrixDSym &cov, 
                                 int trackid, double mass, vertexNumber jv);

      bool vertexPoint_2d       (vertexNumber jv1, vertexNumber jv2);
      bool vertexPoint_3d       (vertexNumber jv1, vertexNumber jv2);
      bool vertexPoint_1track   (vertexNumber jv1, vertexNumber jv2);
      bool vertexPoint_0track   (vertexNumber jv1, vertexNumber jv2);
      bool conversion_2d        (vertexNumber jv);
      bool conversion_3d        (vertexNumber jv);
      bool massConstrain        (int ntrk, const int trkIds[], float mass);

      void setPrimaryVertex     (float xv, float yv, float zv);
      void setPrimaryVertex     (CLHEP::Hep3Vector pv);
      bool setPrimaryVertexError(const CLHEP::HepSymMatrix &xverr);
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
      bool beamlineConstraint(float xb, float yb, CLHEP::HepSymMatrix berr,
                              float xzbslope,float yzbslope);
      bool beamlineConstraint(CLHEP::Hep3Vector pv, CLHEP::HepSymMatrix berr, 
                              float xzbslope, float yzbslope);

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
      FourVector       getTrackP4       (const int trkId) const;

      // return fit mass and get error
      float            getMass          (int ntrk, const int trkIds[], float& dmass) const;

      // return decay length
      float            getDecayLength   (vertexNumber nv, vertexNumber mv, 
                                         const CLHEP::Hep3Vector& dir,
					 float& dlerr) const;
      float            getDecayLength   (vertexNumber nv, vertexNumber mv, const ThreeVector& dir,
					 float& dlerr) const;                                      
      float            getZDecayLength  (vertexNumber nv, vertexNumber mv,
					 const CLHEP::Hep3Vector& dir, float& dlerr) const;  
      float            getZDecayLength  (vertexNumber nv, vertexNumber mv,
					 const ThreeVector& dir, float& dlerr) const;                
      float            getImpactPar     (vertexNumber prdV, vertexNumber dcyV,
					 const CLHEP::Hep3Vector &v, float &dxyerr) const;     
      float            getImpactPar     (vertexNumber prdV, vertexNumber dcyV,
					 const ThreeVector &v, float &dxyerr) const;                 
      float            get_dr           (vertexNumber nv, vertexNumber mv, float& drerr) const;
      float            get_dz           (vertexNumber nv, vertexNumber mv, float& dzerr) const;

      // return location of vertex
      CLHEP::Hep3Vector   getVertexHep     (vertexNumber nv) const;
      ThreeVector         getVertex        (vertexNumber nv) const;

      // return error matrix element.
      ThreeSymMatrix      getErrorMatrix   (vertexNumber nv) const;
      double              getErrorMatrixHep(int j, int k) const;
      CLHEP::HepSymMatrix getErrorMatrixHep(vertexNumber nv) const;
      CLHEP::HepSymMatrix getErrorMatrixHep(const int trkId) const;
      void                getPosMomErr     (CLHEP::HepMatrix& errors) const;
      int              vOff             (vertexNumber jv) const;
      int              tOff             (const int trkId) const;
      int              pOff             (int lp) const;
      int              cOff             (int lc) const;
      int              mOff             () const;
      double           VMat             (int i, int j) const;
      float            getPtError       (const int trkId) const;
      MultiVertexFitterC::vertexNumber
	               allocateVertexNumber();
      void             resetAllocatedVertexNumber();

      // Accessors for getting information relative to ijk errors.
      // Get the error code from the three ijk indexes into the argument variables
      void             getIJKErr(int& err0, int& err1, int& err2) const;
      // Return each error code from the three ijk indexes
      int              getIJKErr0() const;
      int              getIJKErr1() const;
      int              getIJKErr2() const;

      // Get the track-id of the track causing a fatal error as indicated 
      // by the corresponding ijk error
      int              getErrTrackId() const;

      // Set new track reference point
      void             setTrackReferencePoint(const ThreeVector &ref);

      //--------------------------------------------------------------------------------------------
      // Overload operators
      //--------------------------------------------------------------------------------------------
      friend std::ostream& operator << (std::ostream& os, const MultiVertexFitterC& vfit);

    protected:
      std::string      _expert;                  // string: name and email of expert
      int              _stat;                    // status returned from fit

      static const int _maxvtx = CCTVMFT_MAXVTX; // Maximum number of vertices
      static const int _maxmcn = CCTVMFT_MAXMCN; // Maximum number of mass constraints
      static const int _maxtrk = CCTVMFT_MAXTRK; // Maximum number of tracks
      static const int _maxitr = CCTVMFT_MAXITR; // Maximum number of iteration steps
      static const int _maxdim = 5*(CCTVMFT_MAXVTX+1)+3*CCTVMFT_MAXTRK+CCTVMFT_MAXMCN;

      // FIDDLE must have access to protected data like _maxvtx, etc. It must therefore be a
      // friend
      struct        FIDDLE;
      friend struct FIDDLE;
      struct FIDDLE {
          int    excuse;
          float  chisqmax;
      };

      // CTVMQ must have access to MultiVertexFitterC's (other) protected data like _maxvtx, etc.
      // It must therefore be a friend.
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

      // CTVMFR must have access to MultiVertexFitterC's (other) protected data like _maxdim, etc. 
      // It must therefore be a friend.
      struct        CTVMFR;
      friend struct CTVMFR;
      struct        CTVMFR {
          double vmat[_maxdim+1][_maxdim];
      };

      // TRKPRM must have access to MultiVertexFitterC's (other) protected data like _maxtrk, etc.
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

    private:
      double            _bField;                       // B field in Tesla
      int               _currentAllocatedVertexNumber; // index to enum vertexNumber
      ThreeVector       _referencePoint;               // reference point of track
      CLHEP::Hep3Vector _primaryVertex;                // primary vertex relative to _referencePoint
      CLHEP::Hep3Vector _cdfPrimaryVertex;             // primary vertex in CDF coordinate system
      bool              _extrapolateTrackErrors;       // extrapolate track errors _referencePoint
  };
}

//--------------------------------------------------------------------------------------------------
// Inline functions
//--------------------------------------------------------------------------------------------------
inline
std::ostream& operator << (std::ostream& os, const mithep::MultiVertexFitterC& vfit)
{
  vfit.print(os);
  return    (os);
}

//--------------------------------------------------------------------------------------------------
inline
void mithep::MultiVertexFitterC::setExcuse()
{
  _fiddle.excuse = 1;   // the default
}

//--------------------------------------------------------------------------------------------------
inline
void mithep::MultiVertexFitterC::setNoExcuse()
{
  _fiddle.excuse = 0;   // crash on input error
}

//--------------------------------------------------------------------------------------------------
inline
void mithep::MultiVertexFitterC::setChisqMax(const float chisqmx)
{
  _fiddle.chisqmax = chisqmx;
}
#endif
