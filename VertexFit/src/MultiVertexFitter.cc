// $Id: MultiVertexFitter.cc,v 1.4 2009/03/20 13:33:03 loizides Exp $

#include "MitCommon/VertexFit/interface/MultiVertexFitter.h"
#include "MitCommon/Ctvmft/interface/common_blocks.hh"
#include <algorithm>
#include <math.h>
#include <iostream>
#include <csignal>
#include <csetjmp>
#include <TMath.h>
#include <CLHEP/Matrix/Matrix.h>

extern "C" void  ctvmft_(int&, int&, int&);
extern "C" bool  mcalc_ (int&, int*, float&, float&, double*);
extern "C" void  dcalc_ (int&, int&, float*, float&, float&, float*);

using namespace std;
using namespace mithep;
using namespace CLHEP;

jmp_buf env;

//--------------------------------------------------------------------------------------------------
extern "C" void MultiVertexFitterSetStatus(int i) { 
  cout << "Warning, you are handling a severe error in MultiVertexFitter" << endl;
  longjmp(env,-66);
}

//--------------------------------------------------------------------------------------------------
MultiVertexFitter::MultiVertexFitter() :
  _bField                      (3.8),   // default B field for running
  _currentAllocatedVertexNumber(0),     // facilitates CandNode recursion
  _referencePoint              (0,0,0), // set reference point to (0,0,0) initially
  _primaryVertex               (0,0,0), // set primary vertex to (0,0,0) initially
  _cdfPrimaryVertex            (0,0,0)  // set pv in CMS coords to (0,0,0) initially
{
  // Set name and email of MultiVertexFitter expert
  _expert="Christoph Paus (paus@mit.edu)";

  // First get pointers to various FORTAN common blocks
  // _ctvmq_com  = (CTVMQ*)  ctvmq_address_();  //printf(" Common:  _ctvmq_com   %p\n", _ctvmq_com );
  // _ctvmfr_com = (CTVMFR*) ctvmfr_address_(); //printf(" Common:  _ctvmfr_com  %p\n", _ctvmfr_com);
  // _fiddle_com = (FIDDLE*) fiddle_address_(); //printf(" Common:  _fiddle_com  %p\n", _fiddle_com);
  // _trkprm_com = (TRKPRM*) trkprm_address_(); //printf(" Common:  _trkprm_com  %p\n", _trkprm_com);
  _ctvmq_com = (CTVMQ*)ctvmfr_;
  _ctvmfr_com = (CTVMFR*)ctvmfd_;
  _fiddle_com = (FIDDLE*)dctvmfi_;
  _trkprm_com = (TRKPRM*)dctvmtr_;

  // Initialize various arrays
  init();
  // Don't bomb program on error.
  _fiddle.excuse          = 1;
  // Don't extrapolate track errors by default
  _extrapolateTrackErrors = false;
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::init(double bField)
{
  // Set internal variable which keeps track of the b field
  _bField = bField;

  // Now initialize CTVMQ common.  Run and trigger numbers are dummies - they are not used.
  _ctvmq.runnum    = 1;
  _ctvmq.trgnum    = 100;
  // Eventually, we have to get the magnetic field from the right place.
  // Origignal with      bmag = 14.116 [kGauss]:
  // _ctvmq.pscale = 0.000149896 * bmag;
  // New bfield in Tesla bmag = 1.4116 [T]:
  _ctvmq.pscale    = 0.00149896 * bField;
  // Set the default maximum chi-square per dof.
  _fiddle.chisqmax = 225.0;
  // We also need to get the primary vertex from the right place, but for now we put in (0,0,0).
  setPrimaryVertex(0.0,0.0,0.0);
  float xverr[3][3];
  for (int j = 0; j<3; j++) {
    for (int k = 0; k<3; k++) {
      xverr[j][k] = 0.;
    }
  }
  xverr[0][0]   = 0.005;
  xverr[1][1]   = 0.005;
  xverr[2][2]   = 1.0;
  setPrimaryVertexError(xverr);
  // Zero number of tracks, vertices and mass constraints
  _ctvmq.ntrack = 0;
  _ctvmq.nvertx = 0;
  _ctvmq.nmassc = 0;
  // Zero track list and arrays containing the vertex and mass constraint configuration
  for (int j=0; j<_maxtrk; ++j) {
    _ctvmq.list[j]=0;
    for (int jv=0; jv<_maxvtx; ++jv)
      _ctvmq.trkvtx[jv][j] = false;
    for (int jmc=0; jmc<_maxmcn; ++jmc)
      _ctvmq.trkmcn[jmc][j] = false;
  }
  // Initialize the conversion and vertex pointing arrays
  for (int jv=0; jv<_maxvtx; ++jv) {
    _ctvmq.cvtx[jv]      =  0;
    _ctvmq.vtxpnt[0][jv] = -1;
    _ctvmq.vtxpnt[1][jv] =  0;
  }
  _ctvmq.drmax      =  2.0;
  _ctvmq.dzmax      = 20.0;
  _ctvmq.rvmax      = 70.0;
  _ctvmq.trnmax     =  0.5;
  _ctvmq.dsmin      = -100.0; // for CDF used -2, for CMS fuzzing around.
  // Set _stat to -999  and chisqq to -1 to symbolize that no fit has yet been done.
  _stat              = -999;
  _ctvmq.chisqr[0]  = -1.0;
  
  _primaryVertex    = Hep3Vector(0,0,0);
  _cdfPrimaryVertex = Hep3Vector(0,0,0);
  _referencePoint   = ThreeVector(0,0,0);
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::addTrack(const HepVector &v, const HepSymMatrix &cov,
				 int trackId, float mass, vertexNumber jv)
{
  // Check that this vertex number is within the allowed range.
  if (jv<VERTEX_1 || jv>_maxvtx)
    return false;
  _ctvmq.nvertx = jv>_ctvmq.nvertx ? jv : _ctvmq.nvertx;
  
  // Add track extrapolation, if the user set it
  HepVector    w = v;
  HepSymMatrix m = cov;
  //if (_extrapolateTrackErrors)
  //  // This will move the reference point
  //  moveReferencePoint(w, m);
  
  // Check that we have not exceeded the maximum number of tracks.
  if (_ctvmq.ntrack>=_maxtrk)
    return false;
  
  // Add this track
  _ctvmq.list  [_ctvmq.ntrack]       = trackId;
  _ctvmq.tkbank[_ctvmq.ntrack][0]    = 'Q';
  _ctvmq.tkbank[_ctvmq.ntrack][1]    = 'T';
  _ctvmq.tkbank[_ctvmq.ntrack][2]    = 'R';
  _ctvmq.tkbank[_ctvmq.ntrack][3]    = 'K';
  _ctvmq.tmass [_ctvmq.ntrack]       = mass;
  _ctvmq.trkvtx[jv-1][_ctvmq.ntrack] = true;
  
  // Put this track's helix parameters and error matrix into a fortran common block so that they
  // can be accessed by gettrk.  This is a dummy for now.
  _trkprm.trhelix[_ctvmq.ntrack][0]  = w[0];
  _trkprm.trhelix[_ctvmq.ntrack][1]  = w[1];
  _trkprm.trhelix[_ctvmq.ntrack][2]  = w[2];
  _trkprm.trhelix[_ctvmq.ntrack][3]  = w[3];
  _trkprm.trhelix[_ctvmq.ntrack][4]  = w[4];
  
  for (int j=0; j<5; ++j) {
    for (int k=0; k<5; ++k)
      _trkprm.trem[_ctvmq.ntrack][j][k]=m[j][k];
  }
  _ctvmq.ntrack++;
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::addTrack(const TVectorD &v, const TMatrixDSym &cov,
				 int trackId, double mass, vertexNumber jv)
{
  // Check that this vertex number is within the allowed range.
  if (jv<VERTEX_1 || jv>_maxvtx)
    return false;
  _ctvmq.nvertx = jv>_ctvmq.nvertx ? jv : _ctvmq.nvertx;
  
  //// Add track extrapolation, if the user set it
  //HepVector    w = v;
  //HepSymMatrix m = cov;
  TVectorD    w(v);
  TMatrixDSym m(cov);
  ////if (_extrapolateTrackErrors)
  ////  // This will move the reference point
  ////  moveReferencePoint(w, m);
  
  // Check that we have not exceeded the maximum number of tracks.
  if (_ctvmq.ntrack>=_maxtrk)
    return false;
  
  // Add this track
  _ctvmq.list  [_ctvmq.ntrack]       = trackId;
  _ctvmq.tkbank[_ctvmq.ntrack][0]    = 'Q';
  _ctvmq.tkbank[_ctvmq.ntrack][1]    = 'T';
  _ctvmq.tkbank[_ctvmq.ntrack][2]    = 'R';
  _ctvmq.tkbank[_ctvmq.ntrack][3]    = 'K';
  _ctvmq.tmass [_ctvmq.ntrack]       = mass;
  _ctvmq.trkvtx[jv-1][_ctvmq.ntrack] = true;
  
  // Put this track's helix parameters and error matrix into a fortran common block so that they
  // can be accessed by gettrk.  This is a dummy for now.
  _trkprm.trhelix[_ctvmq.ntrack][0]  = w[0];
  _trkprm.trhelix[_ctvmq.ntrack][1]  = w[1];
  _trkprm.trhelix[_ctvmq.ntrack][2]  = w[2];
  _trkprm.trhelix[_ctvmq.ntrack][3]  = w[3];
  _trkprm.trhelix[_ctvmq.ntrack][4]  = w[4];
  
  for (int j=0; j<5; ++j) {
    for (int k=0; k<5; ++k)
      _trkprm.trem[_ctvmq.ntrack][j][k]=m[j][k];
  }
  _ctvmq.ntrack++;
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::vertexPoint_2d(vertexNumber jv1, vertexNumber jv2)
{
  // Check that these vertex numbers are within allowed range and that the vertices are unique.
  if (jv1>_maxvtx || jv1<VERTEX_1)
    return false;
  if (jv2>_maxvtx || jv2<PRIMARY_VERTEX)
    return false;
  if (jv1 <= jv2)
    return false;
  
  // Setup vertex pointing.
  _ctvmq.vtxpnt[0][jv1-1] = jv2;
  _ctvmq.vtxpnt[1][jv1-1] = 1;           // 2d pointing.
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::vertexPoint_3d(vertexNumber jv1, vertexNumber jv2)
{
  // Check that these vertex numbers are within allowed range and that the vertices are distinct
  if (jv1>_maxvtx || jv1<VERTEX_1)
    return false;
  if (jv2>_maxvtx || jv2<PRIMARY_VERTEX)
    return false;
  if (jv1 <= jv2)
    return false;
  
  // Setup vertex pointing
  _ctvmq.vtxpnt[0][jv1-1] = jv2;
  _ctvmq.vtxpnt[1][jv1-1] = 2;           // 3d pointing
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::vertexPoint_1track(vertexNumber jv1, vertexNumber jv2)
{
  // Check that these vertex numbers are within allowed range and are distinct
  if (jv1>_maxvtx || jv1<VERTEX_1)
    return false;
  if (jv2>_maxvtx || jv2<PRIMARY_VERTEX)
    return false;
  if (jv1 <= jv2)
    return false;
  
  // Setup vertex pointing
  _ctvmq.vtxpnt[0][jv1-1] = jv2;
  _ctvmq.vtxpnt[1][jv1-1] = 3;           // Point to 1 track vertex
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::vertexPoint_0track(vertexNumber jv1, vertexNumber jv2)
{
  // jv2 is the zero track vertex. jv1 is the multi track vertex which points to jv2

  // Note: You must call this routine at least twice in order for ctvmft_zerotrackvtx to work ie
  // there must be at least 2 vertices pointed at a zero track vertex. ctvmft for this, but the
  // error message may not make it to your log file (look at the local variables in the stack frame
  // especially IJKERR(2). The significance of this error code is documented at the top of whatever
  // routine chucked you out (ctvmf00 in this case)

  // see ctvmft.f source file for discussion. See especially comments at the top of subroutines:
  // ctvmft and ctvmfa

  // Check that these vertex numbers are within allowed range and are distinct.
  if (jv1>_maxvtx || jv1<VERTEX_1)
    return false;
  if (jv2>_maxvtx || jv2<PRIMARY_VERTEX)
    return false;
  if (jv1 <= jv2)
    return false;
  
  // Setup vertex pointing.
  _ctvmq.vtxpnt[0][jv1-1] = jv2;
  _ctvmq.vtxpnt[1][jv1-1] = 4;           // Point to 0 track vertex.
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::conversion_2d(vertexNumber jv)
{
  if (jv<VERTEX_1 || jv>_ctvmq.nvertx)
    return false;

  _ctvmq.cvtx[jv-1] = 1;

  return true;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::conversion_3d(vertexNumber jv)
{
   if (jv<VERTEX_1 || jv>_ctvmq.nvertx)
      return false;

   _ctvmq.cvtx[jv-1] = 2;

   return true;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::massConstrain(int ntrk, const int trkIds[], float mass)
{
  // Check that we have not exceeded the allowed number of mass constraints.
  if (_ctvmq.nmassc>=_maxmcn)
    return false;

  // Set constraint mass
  _ctvmq.cmass[_ctvmq.nmassc]=mass;
     
  // For each track in contraint, set trkmcn true.  Since the number in tracks[] is the track
  // number, we have to find each track in the list of tracks.
  for (int jt=0; jt<ntrk; ++jt) {
    bool found=false;
    for (int kt=0; kt<_ctvmq.ntrack; ++kt) {
      if (trkIds[jt] == _ctvmq.list[kt]) {
	_ctvmq.trkmcn[_ctvmq.nmassc][kt]=true;
	found=true;
      }
    }
    if (!found)
      return false;
  }

  // Increment number of mass constraints.   
  _ctvmq.nmassc++;
  
  return true;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::beamlineConstraint(float xb, float yb, HepSymMatrix berr,
					   float xzbslope, float yzbslope)
{
  // Set beam position at z=0
  setPrimaryVertex(xb,yb,0);
  //if (_extrapolateTrackErrors) {
  //  float newXb = xb - _referencePoint.x() + _referencePoint.z() * xzbslope;
  //  float newYb = yb - _referencePoint.y() + _referencePoint.z() * yzbslope;
  //  setPrimaryVertex(newXb, newYb, 0);
  //}
  
  bool success = setPrimaryVertexError(berr);

  // Set the beamline slope values
  _ctvmq.xzslope = xzbslope;
  _ctvmq.yzslope = yzbslope;

  // Turn ON beamline constraint
  _ctvmq.vtxpnt[0][0] = -100; 
  
  return success;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::beamlineConstraint(Hep3Vector pv, HepSymMatrix berr, float xzbslope,
					   float yzbslope)
{
  // Check if input beam position coordinates are at z=0
  if (pv.z() != 0)
    return false;

  return beamlineConstraint(pv.x(),pv.y(),berr,xzbslope,yzbslope);
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::setPrimaryVertex(float xv, float yv, float zv)
{
  // Set x,y,z position of the primary vertex.
  _ctvmq.xyzpv0[0] = xv;
  _ctvmq.xyzpv0[1] = yv;
  _ctvmq.xyzpv0[2] = zv;
  
  _primaryVertex = Hep3Vector( xv, yv, zv );
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::setPrimaryVertex(Hep3Vector pv)
{
  // Set x,y,z position of the primary vertex.
  _ctvmq.xyzpv0[0] = pv.x();
  _ctvmq.xyzpv0[1] = pv.y();
  _ctvmq.xyzpv0[2] = pv.z();
  
  _primaryVertex   = pv;
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::setPrimaryVertexError(const float xverr[3][3])
{
  // Set the error matrix for the primary vertex.
  for (int j=0; j<3; ++j) {
    for (int k=0; k<3; ++k)
      _ctvmq.exyzpv[j][k]=xverr[j][k];
  }
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::setPrimaryVertexError(const HepSymMatrix &xverr)
{
  // Set the error matrix for the primary vertex using a HepSymMatrix.  First check that the matrix
  // is the correct size.
  if (xverr.num_row() != 3)
    return false;
  for (int j=0; j<3; j++) {
    for (int k=0; k<3; k++)
      _ctvmq.exyzpv[j][k]=xverr[j][k];
  }
  return true;
}

//--------------------------------------------------------------------------------------------------
bool MultiVertexFitter::fit()
{
  // Check that the diagonal elements of all the track error matrices are positive
  bool mstat = true;
  for (int trk=0; trk<_ctvmq.ntrack; ++trk) {
    for (int j=0; j<5; ++j) {
      // Check diagonal elements of error matrix.
      if (_trkprm.trem[trk][j][j] < 0.) {
	// The covariance matrix could not be inverted:  Set the error codes and fail this fit
        mstat = false;
        _ctvmq.ijkerr[0] = 3;
        _ctvmq.ijkerr[1] = 2;
        _ctvmq.ijkerr[2] = trk + 1;
      }
    }
    // Check that curvature of track is reasonable: Pt is above ~10MeV/c.  If not, set the error
    // codes and fail this fit
    if (fabs(_trkprm.trhelix[trk][1]) > 0.1) {
      //if (fabs(_trkprm.trhelix[trk][1]) > 0.01) {
      mstat = false;
      _ctvmq.ijkerr[0] = 3;
      _ctvmq.ijkerr[1] = 5;
      _ctvmq.ijkerr[2] = trk + 1;
    }
  }
  // If there was a problem with any track, fail the fit
  if (!mstat) {
    _stat = 1;
    return false;
  }

  // First copy information into CTVMFT common blocks
  *_ctvmq_com  = _ctvmq;
  *_ctvmfr_com = _ctvmfr;
  *_fiddle_com = _fiddle;
  *_trkprm_com = _trkprm;
  // Do the vertex fit.
  int print    = 0;
  int level    = 0;
  
// #if ( defined(LINUX) && defined(__USE_BSD) ) || defined(OSF1)
//   struct sigaction myaction = {MultiVertexFitterSetStatus, 0, 0, 0}, oldaction;
//   sigaction(SIGFPE, &myaction, &oldaction);
//   if (setjmp(env)!=0) {
//     sigaction(SIGFPE, &oldaction,0);
//     return -999;
//   }
// #endif
  
  ctvmft_(print,level,_stat);

// #if ( defined(LINUX) && defined(__USE_BSD) ) || defined(OSF1)
//   sigaction(SIGFPE, &oldaction,0);
// #endif
  
  // Now copy information from CTVMFT common blocks to local storage
  _ctvmq  = *_ctvmq_com;
  _ctvmfr = *_ctvmfr_com;
  _fiddle = *_fiddle_com;
  _trkprm = *_trkprm_com;
  
  return (_stat == 0);
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::print() const
{
  print(cout);
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::print(ostream& os) const
{
  os << "****************************** MultiVertexFitter "
     << "******************************" << endl;
  os << "Number of tracks: " << _ctvmq.ntrack << endl;
  os << "   Tracks: ";
  for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
    if (jt != 0) os << ",  ";
    os << _ctvmq.list[jt];
  }
  os << endl;
  os << "Number of vertices: " << _ctvmq.nvertx << endl;
  for (int jv=0; jv<_ctvmq.nvertx; ++jv) {
    os << "   Vertex " << jv+1 << " tracks: ";
    for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
      if (_ctvmq.trkvtx[jv][jt]) {
	os << " " << _ctvmq.list[jt];
      }
    }
    os << endl;
  }
  for (int jv=0; jv<_ctvmq.nvertx; ++jv) {
    if (_ctvmq.vtxpnt[0][jv]==0) {
      os << "   Vertex " << jv+1 << " points to the primary vertex ";
    }
    else if (_ctvmq.vtxpnt[0][jv]>0) {
      os << "   Vertex " << jv+1 << " points to vertex "
	 << _ctvmq.vtxpnt[0][jv];
    }
    if (_ctvmq.vtxpnt[1][jv]==1) {
      os << " in 2 dimensions" << endl;
    }
    else if (_ctvmq.vtxpnt[1][jv]==2) {
      os << " in 3 dimensions" << endl;
    }
    else if (_ctvmq.vtxpnt[1][jv]==3) {
      os << ", a single track vertex" << endl;
    }
    if (_ctvmq.cvtx[jv]>0) {
      os << "   Vertex " << jv+1 << " is a conversion" << endl;
    }
  }
  os << "Number of mass constraints: " << _ctvmq.nmassc << endl;
  for (int jmc=0; jmc<_ctvmq.nmassc; ++jmc) {
    os << "   Tracks ";
    for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
      if (_ctvmq.trkmcn[jmc][jt]) {
	os << " " << _ctvmq.list[jt];
      }
    }
    os << " constrained to mass " << _ctvmq.cmass[jmc]
       << " Gev/c^2" << endl;
  }
  if (_stat==-999) {
    os << "No fit has been done." << endl;
  }
  else {
    os << "***** Results of Fit *****" << endl;
    printErr(os);
    os << "   Status = " << _stat << endl;
    os.precision(7);
    os << "   Chi-square = " << scientific << _ctvmq.chisqr[0]
       << " for " << _ctvmq.ndof << " degrees of freedom." << endl;
    os << "   => probability = " << prob() << endl;
    for (int jv=0; jv<_ctvmq.nvertx; ++jv) {
      os << "Vertex " << jv+1
	 << " position: " << scientific
	 << _ctvmq.xyzvrt[jv+1][0] << " "
	 << _ctvmq.xyzvrt[jv+1][1] << " "
	 << _ctvmq.xyzvrt[jv+1][2] << endl;
    }
    for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
      os << "Track " << _ctvmq.list[jt]
	 << " - P4: " << scientific
	 << _ctvmq.trkp4[0][jt] << " "
	 << _ctvmq.trkp4[1][jt] << " "
	 << _ctvmq.trkp4[2][jt] << " "
	 << _ctvmq.trkp4[3][jt] << " "
	 << " - PT: " << scientific
	 << sqrt(_ctvmq.trkp4[0][jt]*_ctvmq.trkp4[0][jt]+
		 _ctvmq.trkp4[1][jt]*_ctvmq.trkp4[1][jt]) << endl;
    }
  }
  os << "****************************************"
     << "**************************" << endl;
  
  return;
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::printErr() const
{
  printErr(cout);
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::printErr(ostream& os) const
{
  os << "MultiVertexFitter: IJKERR = " << _ctvmq.ijkerr[0] << ", "
     << _ctvmq.ijkerr[1] << ", "
     << _ctvmq.ijkerr[2] << endl;
  if (status()==0 && _ctvmq.ijkerr[0]==0) return;
  if (_ctvmq.ijkerr[0] == -1) {
    os << "   Problem with GETTRK:  track requested is not in list."
       << endl
       << "   This should not happen - Contact MultiVertexFitter expert "
       << _expert << "." <<endl;
  }
  else if (_ctvmq.ijkerr[0]==1) {
    os << "   Problem in CTVM00:" << endl;
    if (_ctvmq.ijkerr[1]==1) {
      os << "      Number of tracks is " << _ctvmq.ntrack
	 << "." << endl;
      if (_ctvmq.ntrack < 2) {
	os << ", which is too few (must be at least 2)." << endl;
      }
      else if (_ctvmq.ntrack > _maxtrk) {
	os << ", which is too many (maximum is " << _maxtrk
	   << ")." << endl;
      }
      else {
	os << "      Problem with number of tracks"
	   << " for unknown reasons." << endl;
      }
    }
    else if (_ctvmq.ijkerr[1]==2) {
      os << "      Number of vertices is " << _ctvmq.nvertx
	 << "." << endl;
      if (_ctvmq.nvertx < 1) {
	os << ", which is too few (must be at least 1)." << endl;
      }
      else if (_ctvmq.nvertx > _maxvtx) {
	os << ", which is too many (maximum is " << _maxvtx
	   << ")." << endl;
      }
      else {
	os << endl << "      Problem with number of vertices"
	   << " for unknown reasons." << endl;
      }
    }
    else if (_ctvmq.ijkerr[1]==3) {
      os << "      Number of mass constraints is " << _ctvmq.nmassc
	 << "." << endl;
      if (_ctvmq.nmassc < 0) {
	os << ", which is negative." << endl;
      }
      else if (_ctvmq.nmassc > _maxmcn) {
	os << ", which is too many (maximum is " << _maxmcn
	   << ")." << endl;
      }
      else {
	os << endl << "      Problem with number of mass"
	   << " constraints for unknown reasons." << endl;
      }
    }
    else if (_ctvmq.ijkerr[1]==11) {
      os << "      Vertex " << _ctvmq.ijkerr[2]
	 << " has less than one track." << endl;
    }
    else if (_ctvmq.ijkerr[1]==12) {
      os << "      Vertex " << _ctvmq.ijkerr[2]
	 << " is a conversion vertex with a number of tracks"
	 << " different than two." << endl;
    }
    else if (_ctvmq.ijkerr[1]==13) {
      os << "      Vertex " << _ctvmq.ijkerr[2]
	 << " is a one track vertex that has no multi-track"
	 << " descendents." << endl;
    }
    else if (_ctvmq.ijkerr[1]==14) {
      os << "      Vertex " << _ctvmq.ijkerr[2]
	 << " does not point at a vertex with a lower number."
	 << endl;
    }
    else if (_ctvmq.ijkerr[1]==15) {
      os << "      Vertex " << _ctvmq.ijkerr[2]
	 << " has a parent vertex that is a conversion." << endl;
    }
    else if (_ctvmq.ijkerr[1]==16) {
      os << "      Vertex " << _ctvmq.ijkerr[2]
	 << " does 1 track pointing to a vertex with"
	 << " more than 1 track." << endl;
    }
    else if (_ctvmq.ijkerr[1]==17) {
      os << "      Vertex " << _ctvmq.ijkerr[2]
	 << " does 0 track pointing to a vertex with"
	 << " more than 0 track (?)." << endl;
    }
    else if (_ctvmq.ijkerr[1]==19) {
      os << "      Primary vertex error matrix is singular." << endl;
    }
    else if (_ctvmq.ijkerr[1]==21) {
      os << "      Track with Id " << _ctvmq.ijkerr[2]
	 << "is not in any vertex." << endl;
    }
    else if (_ctvmq.ijkerr[1]==22) {
      os << "      Track with Id " << _ctvmq.ijkerr[2]
	 << "is in multiple vertices." << endl;
    }
    else if (_ctvmq.ijkerr[1]==23) {
      os << "      Track with Id " << _ctvmq.ijkerr[2]
	 << "occurs more than once." << endl;
    }
    else if (_ctvmq.ijkerr[1]==31) {
      os << "      A mass constraint has less than 2 tracks." << endl;
    }
    else if (_ctvmq.ijkerr[1]==32) {
      os << "      The sum masses of the tracks in a mass constraint"
	 << " exceeds the constraint mass." << endl;
    }
    else if (_ctvmq.ijkerr[1]==33) {
      os << "      Beamline constraint. Beam covariance not set properly."
	 << " Negative diagonal elements." << endl;
    }
    else if (_ctvmq.ijkerr[1]==34) {
      os << "      Beamline constraint. Beam covariance not set properly."
	 << " Off-diagonal elements not zero." << endl;
    }
    else if (_ctvmq.ijkerr[1]==36) {
      os << "      Beamline constraint. Number of vertices = " 
	 << _ctvmq.nvertx << " Should be 1." << endl;
    }
  }
  else if (_ctvmq.ijkerr[0] == 2) {
    if (_ctvmq.ijkerr[1] == 20) {
      os << "   Problem in CTVM00: " << endl;
      os << "      Track has negative Id = "
	 << _ctvmq.list[_ctvmq.ijkerr[2]-1] << "." << endl;
    }
    else {
      os << "   Problem in CTVMFA with vertex "
	 << _ctvmq.ijkerr[2] << ": " << endl;
      os << "      Failure in vertex first approximation." << endl;
      if (_ctvmq.ijkerr[1] == 1) {
	os << "      Tracks are concentric circles." << endl;
      }
      if (_ctvmq.ijkerr[1] == 2) {
	os << "      Conversion vertex has widely separated"
	   << " exterior circles at midpoint." << endl;
      }
      if (_ctvmq.ijkerr[1] == 3) {
	os << "      Conversion vertex has widely separated"
	   << " interior circles at midpoint." << endl;
      }
      if (_ctvmq.ijkerr[1] == 4) {
	os << "      Vertex has widely separated"
	   << " exterior circles at approximate vertex." << endl;
      }
      if (_ctvmq.ijkerr[1] == 5) {
	os << "      Vertex has widely separated"
	   << " interior circles at approximate vertex." << endl;
      }
      if (_ctvmq.ijkerr[1] == 6) {
	os << "      Rv is too large at the chosen"
	   << " intersection point." << endl;
      }
      if (_ctvmq.ijkerr[1] == 7) {
	os << "      Delta z is too large at the chosen"
	   << " intersection point." << endl;
      }
      if (_ctvmq.ijkerr[1] == 8) {
	os << "      A track's turning to the chosen vertex"
	   << " is too large." << endl;
      }
      if (_ctvmq.ijkerr[1] == 9) {
	os << "      There is no solution with an adequately"
	   << " positive arc length." << endl;
      }
      if (_ctvmq.ijkerr[1] == 21) {
	os << "      zero-track vertexing: either/both vertex "
	   << " momenta are too small (<0.01 MeV)." << endl;
      }
      if (_ctvmq.ijkerr[1] == 22) {
	os << "      zero-track vertexing: Two lines (tracks) are "
	   << " parallel/antiparallel." << endl;
      }
      
    }
  }
  else if (_ctvmq.ijkerr[0] == 3) {
    os << "   Problem in CTVM01 with track with Id = "
       << _ctvmq.list[_ctvmq.ijkerr[2]-1] << ": " << endl;
    if (_ctvmq.ijkerr[1] == 1) {
      os << "      GETTRK cannot find Id in list." << endl;
    }
    if (_ctvmq.ijkerr[1] == 2) {
      os << "      Covariance matrix could not be inverted."  << endl 
	 << "      Offending track number (in order addded) is "
	 << _ctvmq.ijkerr[2] << "." << endl;
    }
    if (_ctvmq.ijkerr[1] == 3) {
      os << "      Track turns through too large an angle"
	 << " to the vertex." << endl;
    }
    if (_ctvmq.ijkerr[1] == 4) {
      os << "      Track moves too far backward to vertex." << endl;
    }
    if (_ctvmq.ijkerr[1] == 5) {
      os << "      Track with curvature > 0.01." << endl
	 << "      Offending track number is "
	 << _ctvmq.ijkerr[2] << "." << endl;
    }
  }
  else if (status() == 9) {
    os << "   General fit problem: " << endl;
    if (_ctvmq.ijkerr[1] == 1) {
      os << "      Singular solution matrix." << endl;
    }
    if (_ctvmq.ijkerr[1] == 2 || _ctvmq.ijkerr[1] == 3) {
      os << "      Too many iterations ( "
	 << _ctvmq.ijkerr[2] << "(." << endl;
    }
    if (_ctvmq.ijkerr[1] == 4) {
      os << "      Convergence failure." << endl;
    }
    if (_ctvmq.ijkerr[1] == 5) {
      os << "      Bad convergence." << endl;
    }
    if (_ctvmq.ijkerr[1] == 9) {
      os << "      Ill-formed  covariance matrix." << endl;
    }
  }
  else {
    os << "   The error codes above are not recognized." << endl
       << "   Contact MultiVertexFitter expert " << _expert << "." << endl;
  }
  return;
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::getIJKErr(int& err0, int& err1, int& err2) const
{
  err0 = _ctvmq.ijkerr[0];
  err1 = _ctvmq.ijkerr[1];
  err2 = _ctvmq.ijkerr[2];
  return;
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::getIJKErr0() const
{
  return _ctvmq.ijkerr[0];
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::getIJKErr1() const
{
  return _ctvmq.ijkerr[1];
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::getIJKErr2() const
{
  return _ctvmq.ijkerr[2];
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::getErrTrackId() const
{
  if (status() == 0) return 0;
  int trkId = 0;
  // Problems with track in CTVM01 or track has negative id in CTVM00 See PrintErr() for a more
  // detailed list of error codes.
  if ((_ctvmq.ijkerr[0] == 2 && _ctvmq.ijkerr[1] == 20) ||
      _ctvmq.ijkerr[0] == 3) {
    trkId = _ctvmq.list[_ctvmq.ijkerr[2]-1];
  }
  
  return trkId;
}

//--------------------------------------------------------------------------------------------------
string MultiVertexFitter::expert() const
{
  return _expert;
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::status() const
{
  return _stat;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::chisq() const
{
  // Chi-square of fit
  return _ctvmq.chisqr[0];
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::ndof() const
{
  // Number of degrees of freedom of fit.
  if (_ctvmq.chisqr[0] >= 0)
    return _ctvmq.ndof;
  else
    return 0;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::prob() const
{
  // Probability of chi-square of fit
  if (_ctvmq.chisqr[0]>=0.) {
    float chisq = _ctvmq.chisqr[0];
    int   nd    = _ctvmq.ndof;
    return TMath::Prob(chisq,nd);
  }
  else
    return -999.;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::chisq(const int trkId) const
{
  // This method returns the chisquare contribution for one track If fit not successful or not done
  // yet, return -1.
  if (_ctvmq.chisqr[0] < 0)
    return -1.;
  // Look for this track in the list of tracks.
  for (int jt = 0; jt < _ctvmq.ntrack; ++jt) {
    if (trkId == _ctvmq.list[jt]) {
      // Found the track, so return its chisquare contribution.
      return _ctvmq.chit[jt];
    }
  }
  // If track is not in list, return -1.
  return -1.;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::chisq_rphi() const
{
  // This method returns the chisquare contribution in the r-phi plane.
  int index[3] = {0,1,3};
  // If fit not successful or not done yet, return -1.
  if (_ctvmq.chisqr[0] < 0)
    return -1.;
  // Loop over the tracks in the event.
  float chisq = 0.;
  for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
    // Double loop over the relevant parameter indices.
    for (int k1=0; k1<3; ++k1) {
      for (int k2=0; k2<3; ++k2)
	// Add contribution to chisquare.
	chisq += _ctvmq.pardif[jt][index[k1]] *
                 _ctvmq.g[jt][index[k1]][index[k2]] *
                 _ctvmq.pardif[jt][index[k2]];
    }
  }
  // Return the chisquare.
  return chisq;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::chisq_z() const
{
  // This method returns the chisquare contribution in the z direction.
  int index[2] = {2,4};
  // If fit not successful or not done yet, return -1.
  if (_ctvmq.chisqr[0] < 0)
    return -1.;
  // Loop over the tracks in the event.
  float chisq = 0.;
  for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
    // Double loop over the relevant parameter indices.
    for (int k1=0; k1<2; ++k1) {
      for (int k2=0; k2<2; ++k2)
	// Add contribution to chisquare.
	chisq += _ctvmq.pardif[jt][index[k1]] *
                 _ctvmq.g[jt][index[k1]][index[k2]] *
                 _ctvmq.pardif[jt][index[k2]];
    }
  }
  // Return the chisquare.
  return chisq;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::chisq_rphiz() const
{
  // This method returns the chisquare contribution of the cross
  // terms in the r-phi and z directions.
  int index1[2] = {2,4};
  int index2[3] = {0,1,3};
  // If fit not successful or not done yet, return -1.
  if (_ctvmq.chisqr[0] < 0)
    return -1.;
  // Loop over the tracks in the event.
  float chisq = 0.;
  for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
    // Double loop over the relevant parameter indices.
    for (int k1=0; k1<2; ++k1) {
      for (int k2=0; k2<3; ++k2)
	// Add contribution to chisquare.
	chisq += _ctvmq.pardif[jt][index1[k1]] *
	         _ctvmq.g[jt][index1[k1]][index2[k2]] *
	         _ctvmq.pardif[jt][index2[k2]];
     }
   }

   // Return the chisquare.
   return 2.0 * chisq;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::chisq_rphi(const int trkId) const
{
  // This method returns the chisquare contribution in the r-phi plane.
  int index[3] = {0,1,3};
  // If fit not successful or not done yet, return -1.
  if (_ctvmq.chisqr[0] < 0)
    return -1.;
  // Loop over the tracks in the event, looking for the one we want
  for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
    if (trkId == _ctvmq.list[jt]) {
      // Found the track, so calculate its chisquare contribution.
      float chisq = 0.;
      // Double loop over the relevant parameter indices.
      for (int k1=0; k1<3; ++k1) {
	for (int k2=0; k2<3; ++k2) {
	  // Add contribution to chisquare.
	  chisq += _ctvmq.pardif[jt][index[k1]] *
	           _ctvmq.g[jt][index[k1]][index[k2]] *
                   _ctvmq.pardif[jt][index[k2]];
	}
      }
      return chisq;
    }
  }

  // Track not found, return -1.
  return -1.;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::chisq_z(const int trkId) const
{
  // This method returns the chisquare contribution in the z direction.
  int index[2] = {2,4};
  // If fit not successful or not done yet, return -1.
  if (_ctvmq.chisqr[0] < 0)
    return -1.;
  // Loop over the tracks in the event, looking for the one we want.
  for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
    if (trkId == _ctvmq.list[jt]) {
      // Found the track, so calculate its chisquare contribution.
      float chisq = 0.;
      // Double loop over the relevant parameter indices.
      for (int k1=0; k1<2; ++k1) {
	for (int k2=0; k2<2; ++k2)
	  // Add contribution to chisquare.
	  chisq += _ctvmq.pardif[jt][index[k1]] *
        	   _ctvmq.g[jt][index[k1]][index[k2]] *
                   _ctvmq.pardif[jt][index[k2]];
      }
      return chisq;
    }
  }

  // Track not found - return -1.
  return -1.;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::chisq_rphiz(const int trkId) const
{
  // This method returns the chisquare contribution of the cross terms in the r-phi and z
  // directions
  int index1[2] = { 2,4 };
  int index2[3] = { 0,1,3 };
  // If fit not successful or not done yet, return -1
  if (_ctvmq.chisqr[0] < 0)
    return -1.;
  // Loop over the tracks in the event
  for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
    if (trkId == _ctvmq.list[jt]) {
      // Found the track, so calculate its chisquare contribution
      float chisq = 0.;
      // Double loop over the relevant parameter indices
      for (int k1=0; k1<2; ++k1) {
	for (int k2=0; k2<3; ++k2)
	  // Add contribution to chisquare
	  chisq += _ctvmq.pardif[jt][index1[k1]] *
	           _ctvmq.g[jt][index1[k1]][index2[k2]] *
	           _ctvmq.pardif[jt][index2[k2]];
      }
      return 2.0 * chisq;
    }
  }
  // Track not found so return -1.
  return -1.;
}

//--------------------------------------------------------------------------------------------------
FourVector MultiVertexFitter::getTrackP4(const int trkId) const
{
  if (_stat != 0)
    return FourVector(0,0,0,0);
  // return four momentum of fit track 
  for (int jt=0; jt<_ctvmq.ntrack; ++jt) {
    // Find which track matches this Id.
    if (trkId == _ctvmq.list[jt]) {
      FourVector p4((double)_ctvmq.trkp4[0][jt], (double)_ctvmq.trkp4[1][jt],
		    (double)_ctvmq.trkp4[2][jt], (double)_ctvmq.trkp4[3][jt]);
      return p4;
    }
  }
  return FourVector(0,0,0,0);
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::getMass(int ntrk, const int trkIds[], float &dmass) const
{
// #if (defined(LINUX) && defined(__USE_BSD)) || defined(OSF1)
//   struct sigaction myaction = {MultiVertexFitterSetStatus, 0, 0, 0}, oldaction;
//   sigaction(SIGFPE, &myaction, &oldaction);
//   if (setjmp(env)!=0) {
//     sigaction(SIGFPE, &oldaction,0);
//     return -999;
//   }
// #endif
  
  dmass = -999.;
  if (_stat!=0)
    return -999.;
  // Get fit invariant mass of ntrk tracks listed in array tracks.  dmass is the error on the mass.
  dmass=-999.;
  if (ntrk <= 0)
    return 0;
  int jtrks[_maxtrk];
  for (int jt=0; jt<ntrk; ++jt) {
    bool found = false;
    for (int kt=0; kt<_ctvmq.ntrack; ++kt) {
      if (trkIds[jt] == _ctvmq.list[kt]) {
	found = true;
	jtrks[jt]=kt+1;
      }
    }
    if (!found)
      return 0;
  }
  // Copy information into CTVMFT common blocks
  *_ctvmq_com  = _ctvmq;
  *_ctvmfr_com = _ctvmfr;
  int    ntr   = ntrk;
  float  mass;
  double p4[4];
  mcalc_(ntr, jtrks, mass, dmass, p4);
  
// #if (defined(LINUX) && defined(__USE_BSD)) || defined(OSF1)
//   sigaction(SIGFPE, &oldaction,0);
// #endif

  return mass;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::getDecayLength(vertexNumber nv, vertexNumber mv, 
					const Hep3Vector& dir, float& dlerr) const
{
  dlerr = -999.;
  if (_stat!=0)
    return -999.;

  // Get the signed decay length from vertex nv to vertex mv along the x-y direction defined by the
  // 3 vector dir.  dlerr is the error on the decay length including the full error matrix.  Check
  // that vertices are within range.
  if (nv<0 || nv>=_ctvmq.nvertx)
    return -999.;
  if (mv<1 || mv>_ctvmq.nvertx)
    return -999.;
  if (nv>=mv)
    return -999.;

  float dir_t = sqrt(dir.x()*dir.x()+dir.y()*dir.y());
  if (dir_t <= 0.)
    return -999.;

  Hep3Vector dv = getVertexHep(mv)-getVertexHep(nv);
  float      dl = (dv.x()*dir.x()+dv.y()*dir.y())/dir_t;
  // Set up the column matrix of derivatives
  HepMatrix A(4,1);
  A(1,1) =  dir.x()/dir_t;
  A(2,1) =  dir.y()/dir_t;
  A(3,1) = -dir.x()/dir_t;
  A(4,1) = -dir.y()/dir_t;
  // Check if first vertex (nv) is primary vertex.  If it is, check if it was used in the primary
  // vertex.  If not, all of the corresponding error matrix elements are those supplied for the
  // primary vertex.
  int nvf = 0;
  if (nv==0) {
    nvf=-1;
    for (int jv=0; jv<_ctvmq.nvertx; ++jv) {
      if (_ctvmq.vtxpnt[0][jv]==0)
	nvf=0;
    }
  }
  // Get the relevant submatrix of the full error matrix.
  HepMatrix V(4,4,0);
  if (nvf < 0) {
    V(1,1) = getErrorMatrixHep(_ctvmq.voff[mv-1]+1,_ctvmq.voff[mv-1]+1);
    V(1,2) = getErrorMatrixHep(_ctvmq.voff[mv-1]+1,_ctvmq.voff[mv-1]+2);
    V(2,1) = getErrorMatrixHep(_ctvmq.voff[mv-1]+2,_ctvmq.voff[mv-1]+1);
    V(2,2) = getErrorMatrixHep(_ctvmq.voff[mv-1]+2,_ctvmq.voff[mv-1]+2);
    V(3,3) = _ctvmq.exyzpv[0][0];
    V(3,4) = _ctvmq.exyzpv[0][1];
    V(4,3) = _ctvmq.exyzpv[1][0];
    V(4,4) = _ctvmq.exyzpv[1][1];
  }
  else {
    // Get the indices into the error matrix vmat
    int index[4] = { _ctvmq.voff[mv-1]+1,_ctvmq.voff[mv-1]+2,0,0 };
    if (nv==0) {
      index[2] = 1;
      index[3] = 2;
    }
    else {
      index[2] = _ctvmq.voff[nv-1]+1;
      index[3] = _ctvmq.voff[nv-1]+2;
    }
    for (int j=0; j<4; ++j) {
      for (int k=0; k<4; ++k) {
	V[j][k] = getErrorMatrixHep(index[j],index[k]);
      }
    }
  }

  // Calculate square of dlerr
  dlerr = (A.T()*V*A)(1,1);
  if (dlerr >= 0.)
    dlerr = sqrt(dlerr);
  else
    dlerr = -sqrt(-dlerr);
  
  return dl;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::getDecayLength(vertexNumber nv, vertexNumber mv, 
                                        const ThreeVector& dir, float& dlerr) const
{
  Hep3Vector dirHep(dir.x(),dir.y(),dir.z());
  return getDecayLength(nv, mv, dirHep, dlerr);
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::getZDecayLength(vertexNumber nv, vertexNumber mv,
                              const Hep3Vector& mom, float& dlerr) const
{
  //----------------------------------------------------------------------------
  // Get the signed decay length from vertex nv to vertex mv along the
  // z direction of the momentum vector, mom.
  // dlerr is the error on the decay length including the full error
  // matrix.
  //----------------------------------------------------------------------------
  // Start with good initialization
  dlerr = -999.;

  // Check that the fit worked
  if (_stat != 0)
    return -999.;

  // Check that vertices are within range.
  if (nv<0 || nv>=_ctvmq.nvertx)
    return -999.;
  if (mv<1 || mv> _ctvmq.nvertx)
    return -999.;
  if (nv >= mv)
    return -999.;

  // Calculate the vector length
  float length = fabs(mom.z());
  if (length <= 0.)
    return -999.;

  // Get the vector pointing from first vertex (nv) to second vertex (mv)
  Hep3Vector dv = getVertexHep(mv) - getVertexHep(nv);

  //----------------------------------------------------------------------------
  // Calculate the "decay distance"
  //----------------------------------------------------------------------------
  // Project the vertex vector onto the momentum vector direction
  float      dl = (dv.z()*mom.z())/length;

  //----------------------------------------------------------------------------
  // Calculate the error on that distance
  //----------------------------------------------------------------------------
  // Set up the column matrix of derivatives
  HepMatrix A(2,1);
  A(1,1) =  mom.z()/length;
  A(2,1) = -mom.z()/length;

  // Need to catch the special case if the first vertex is the primary
  int nvf = 0;
  if (nv == 0) {
    nvf = -1;
    for (int jv=0; jv<_ctvmq.nvertx; ++jv) {
      if (_ctvmq.vtxpnt[0][jv] == 0)
        nvf = 0;
    }
  }
  // Get the relevant submatrix of the full error matrix.
  HepMatrix V(2,2,0);
  if (nvf < 0) {
    // Geometric uncertainties (positions second vertex)
    V(1,1) = getErrorMatrixHep(_ctvmq.voff[mv-1]+3,_ctvmq.voff[mv-1]+3);
    // Geometric uncertainties (positions first vertex)
    V(2,2) = _ctvmq.exyzpv[2][2];
  }
  else {
    // Get the indices into the error matrix vmat
    int index[2] = { _ctvmq.voff[mv-1]+3,_ctvmq.voff[nv-1]+3 };
    // Handeling the case of the primary vertex
    if (nv == 0)
      index[1] = 3;
    // All right... copy
    for(int j=0; j<2; ++j)
      for(int k=0; k<2; ++k)
        V[j][k] = getErrorMatrixHep(index[j],index[k]);
  }
  // Calculate square of dlerr
  dlerr = (A.T() * V * A )(1,1);
  if (dlerr >= 0.)
    dlerr =  sqrt(dlerr);
  else
    dlerr = -sqrt(-dlerr);

  return dl;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::getZDecayLength(vertexNumber nv, vertexNumber mv,
                              const ThreeVector& mom, float& dlerr) const
{
  Hep3Vector momHep(mom.x(),mom.y(),mom.z());
  return getZDecayLength(nv, mv, momHep, dlerr);
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::getImpactPar(vertexNumber prdV, vertexNumber dcyV,
                                         const Hep3Vector &v, float &dxyerr) const
{
  Hep3Vector   PVtx   = getVertexHep     (prdV);
  Hep3Vector   DVtx   = getVertexHep     (dcyV);
  HepSymMatrix PVtxCv = getErrorMatrixHep(prdV);
  HepSymMatrix DVtxCv = getErrorMatrixHep(dcyV);

  double norma = v.perp();
  if (norma <= 0) {
    dxyerr = -999.;
    return -999.0;
  }
  double dxy = ((v.cross(DVtx-PVtx)).z())/norma;

  // Calculate error on B impact parameter:
  double cosPhi = cos(v.phi());
  double sinPhi = sin(v.phi());
  dxyerr = cosPhi * cosPhi * (DVtxCv[1][1] + PVtxCv[1][1])
    +      sinPhi * sinPhi * (DVtxCv[0][0] + PVtxCv[0][0])
    -      2.0 * cosPhi * sinPhi * (DVtxCv[0][1] + PVtxCv[0][1]);
  dxyerr = (dxyerr>0.0) ? sqrt(dxyerr) : -999.;

  return dxy;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::getImpactPar(vertexNumber prdV, vertexNumber dcyV,
                                         const ThreeVector &v, float &dxyerr) const
{
  Hep3Vector vHep(v.x(),v.y(),v.z());
  return getImpactPar(prdV, dcyV, vHep, dxyerr);
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::get_dr(vertexNumber nv, vertexNumber mv, float& drerr) const
{
   drerr = -999.;
   if (_stat!=0)
     return -999.;
   // Get the transvese distance between vertices nv and mv and return it as the function value.
   // drerr is the uncertainty on the transverse distance, calculated from the full error matrix
   // including correlations.
   float dxyz[3];
   float dr;
   float dz;
   float dl[3];
   
   int mvert    = mv;
   int nvert    = nv;
   // Copy information into CTVMFT common blocks
   *_ctvmq_com  = _ctvmq;
   *_ctvmfr_com = _ctvmfr;
   // Do calculation
   dcalc_(mvert,nvert,dxyz,dr,dz,dl);
   drerr = dl[0];
   return dr;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::get_dz(vertexNumber nv, vertexNumber mv, float& dzerr) const
{
  dzerr = -999.;
  if (_stat!=0)
    return -999.;
  // Get the longitudinal distance between vertices nv and mv and return it as the function value.
  // dzerr is the uncertainty on the longitudinal distance, calculated from the full error matrix
  // including correlations.
  float dxyz[3];
  float dr;
  float dz;
  float dl[3];
  
  int mvert    = mv;
  int nvert    = nv;
  // Copy information into CTVMFT common blocks
  *_ctvmq_com  = _ctvmq;
  *_ctvmfr_com = _ctvmfr;

  // Do calculation
  dcalc_(mvert,nvert,dxyz,dr,dz,dl);
  dzerr = dl[1];
  return dz;
}

//--------------------------------------------------------------------------------------------------
Hep3Vector MultiVertexFitter::getVertexHep(vertexNumber nv) const
{
  if (_stat!=0)
    return Hep3Vector(-999,-999,-999);
  // Return x,y,z position of vertex nv.
  if (nv<0 || nv>_ctvmq.nvertx)
    return Hep3Vector(0,0,0);
  // Check if first vertex (nv) is primary vertex.  If it is, check if it was used in the primary
  // vertex.  If not, all of the corresponding error matrix elements are those supplied for the
  // primary vertex.
  int nvf=0;
  if (nv==0) {
    nvf=-1;
    for (int jv=0; jv<_ctvmq.nvertx; ++jv) {
      if (_ctvmq.vtxpnt[0][jv]==0)
	nvf=0;
    }
  }
  Hep3Vector vertex;
  // If primary vertex was not part of fit, take vertex location as supplied.
  if (nvf < 0) {
    vertex.setX(_ctvmq.xyzpv0[0]);
    vertex.setY(_ctvmq.xyzpv0[1]);
    vertex.setZ(_ctvmq.xyzpv0[2]);
  }
  else {
    vertex.setX(_ctvmq.xyzvrt[nv][0]);
    vertex.setY(_ctvmq.xyzvrt[nv][1]);
    vertex.setZ(_ctvmq.xyzvrt[nv][2]);
  }
  //// If we have a different reference point, need to add it back in
  //vertex += _referencePoint;
  return vertex;
}

//--------------------------------------------------------------------------------------------------
ThreeVector MultiVertexFitter::getVertex(vertexNumber nv) const
{
  if (_stat!=0)
    return ThreeVector(-999,-999,-999);
  // Return x,y,z position of vertex nv.
  if (nv<0 || nv>_ctvmq.nvertx)
    return ThreeVector(0,0,0);
  // Check if first vertex (nv) is primary vertex.  If it is, check if it was used in the primary
  // vertex.  If not, all of the corresponding error matrix elements are those supplied for the
  // primary vertex.
  int nvf=0;
  if (nv==0) {
    nvf=-1;
    for (int jv=0; jv<_ctvmq.nvertx; ++jv) {
      if (_ctvmq.vtxpnt[0][jv]==0)
	nvf=0;
    }
  }
  ThreeVector vertex;
  // If primary vertex was not part of fit, take vertex location as supplied.
  if (nvf < 0) {
    vertex.SetX(_ctvmq.xyzpv0[0]);
    vertex.SetY(_ctvmq.xyzpv0[1]);
    vertex.SetZ(_ctvmq.xyzpv0[2]);
  }
  else {
    vertex.SetX(_ctvmq.xyzvrt[nv][0]);
    vertex.SetY(_ctvmq.xyzvrt[nv][1]);
    vertex.SetZ(_ctvmq.xyzvrt[nv][2]);
  }
  // If we have a different reference point, need to add it back in
  vertex += _referencePoint;
  return vertex;
}

//--------------------------------------------------------------------------------------------------
ThreeSymMatrix MultiVertexFitter::getErrorMatrix(MultiVertexFitter::vertexNumber nv) const 
{  
  // return errors for vertex nv
  ThreeSymMatrix cov;
  
  // if this is the primary vertex, return the error matrix the user supplied
  if (nv==PRIMARY_VERTEX) {
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++) 
        cov(j,k) = _ctvmq.exyzpv[j][k];
    return cov;
  }
  
  if (_stat!=0)
    return cov;
  // return x,y,z position of vertex nv
  if (nv<VERTEX_1 || nv>_ctvmq.nvertx)
    return cov;
  // get offset for vertex nv
  int voff = _ctvmq.voff[nv-1];
  // fill matrix
  for (int i = 0 ; i < 3 ; ++i)
    for (int j = i ; j < 3 ; ++j)
      cov(i,j) = _ctvmfr.vmat[voff+i][voff+j];
  return cov;
}

//--------------------------------------------------------------------------------------------------
double MultiVertexFitter::getErrorMatrixHep(int j, int k) const
{
  if (_stat!=0)
    return -999.;
  // Return the j,k element of the full error matrix VMAT.  Since the CTVMFT documentation assumes
  // the indices start from 1 (ala Fortran), we will also assume this and convert the C++ indices.
  // Note also the that order of Fortran and C++ indices is different.  We assume that j and k are
  // in the Fortran order.
  if (j<1 || k<1 || j>_maxdim+1 || k>_maxdim)
    return -999.;

  return _ctvmfr.vmat[k-1][j-1];
}

//--------------------------------------------------------------------------------------------------
HepSymMatrix MultiVertexFitter::getErrorMatrixHep(MultiVertexFitter::vertexNumber nv) const 
{  
  // return errors for vertex nv
  HepSymMatrix cov(3,0);
  
  // if this is the primary vertex, return the error matrix the user supplied
  if (nv==PRIMARY_VERTEX) {
    for (int j=0; j<3; j++)
      for (int k=0; k<3; k++) 
        cov[j][k] = _ctvmq.exyzpv[j][k];
    return cov;
  }
  
  if (_stat!=0)
    return cov;
  // return x,y,z position of vertex nv
  if (nv<VERTEX_1 || nv>_ctvmq.nvertx)
    return cov;
  // get offset for vertex nv
  int voff = _ctvmq.voff[nv-1];
  // fill matrix
  for (int i = 0 ; i < 3 ; ++i)
    for (int j = i ; j < 3 ; ++j)
      cov[i][j] = _ctvmfr.vmat[voff+i][voff+j];
  return cov;
}

//--------------------------------------------------------------------------------------------------
HepSymMatrix MultiVertexFitter::getErrorMatrixHep(const int trkId) const 
{
  HepSymMatrix cov(3,0);
  if (_stat != 0)
    return cov;
  
  for (int nt=0; nt<_ctvmq.ntrack; ++nt) {
    
    // Find which track matches this Id
    if (trkId == _ctvmq.list[nt]) {
      
      // Position of track nt get offset for track nt
      int toff = _ctvmq.toff[nt];
      
      // Fill matrix -- Crv,Phi,Ctg
      for (int i = 0 ; i < 3 ; ++i)
	for (int j = i ; j < 3 ; ++j)
	  cov[i][j] = _ctvmfr.vmat[toff+i][toff+j];
    }
  }
  return cov;
}

//--------------------------------------------------------------------------------------------------
float MultiVertexFitter::getPtError(const int trkId) const
{  
  if (_stat != 0)
    return 0;
  
  int   toff;
  float pt,curv,curvErr,ptErr;

  for (int nt=0; nt<_ctvmq.ntrack; ++nt) {

    // Find which track matches this Id
    if (trkId == _ctvmq.list[nt]) {

      // Position of track nt get offset for track nt
      toff = _ctvmq.toff[nt];

      // Curvature error
      pt      = sqrt(_ctvmq.trkp4[0][nt]*_ctvmq.trkp4[0][nt] +
		     _ctvmq.trkp4[1][nt]*_ctvmq.trkp4[1][nt]);
      curv    = _ctvmq.pscale/pt;
      curvErr = sqrt(_ctvmfr.vmat[toff+0][toff+0]);
      ptErr   = _ctvmq.pscale/curv/curv*curvErr;
      return ptErr;
    }
  }

  return 0;
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::getPosMomErr(HepMatrix& errors) const
{
  // A c++ rewrite of the FORTRAN MASSAG function The result of this function is an error matrix in
  // position-momentum basis.  A 7x7 matrix of errors where the rows/columns are x, y, z, px, py,
  // pz, e.

  double           cosph [_maxtrk];
  double           sinph [_maxtrk];
  double           cosdph[_maxtrk];
  double           sindph[_maxtrk];
  Hep3Vector       mom3  [_maxtrk];
  HepLorentzVector pmom  [_maxtrk];
  HepLorentzVector total_mom;

  for (int lvtx = 0; lvtx < _ctvmq.nvertx; lvtx++) {
    for (int ltrk = 0; ltrk < _ctvmq.ntrack; ltrk++) {
      if (!_ctvmq.trkvtx[lvtx][ltrk])
	continue;
      cosph[ltrk] = cos(_ctvmq.par[ltrk][1]);
      sinph[ltrk] = sin(_ctvmq.par[ltrk][1]);
      double dphi = 0;
      sindph[ltrk] = 2 * _ctvmq.par[ltrk][0] * 
       (_ctvmq.xyzvrt[lvtx + 1][0] * cosph[ltrk] + 
        _ctvmq.xyzvrt[lvtx + 1][1] * sinph[ltrk]);
      cosdph[ltrk] = sqrt(1.0 - sindph[ltrk] * sindph[ltrk]);
      if (fabs(sindph[ltrk]) <= 1.0){
        dphi = asin(sindph[ltrk]);
      }
      double pt = _ctvmq.pscale * fabs(1./_ctvmq.par[ltrk][0]);
      mom3[ltrk].setX(pt * cos(_ctvmq.par[ltrk][1] + dphi));
      mom3[ltrk].setY(pt * sin(_ctvmq.par[ltrk][1] + dphi));
      mom3[ltrk].setZ(pt * _ctvmq.par[ltrk][2]);
      double e  = sqrt(_ctvmq.tmass[ltrk] * _ctvmq.tmass[ltrk] 
                       + mom3[ltrk].mag2());
      pmom[ltrk].setVect(mom3[ltrk]);
      pmom[ltrk].setE(e);

      total_mom += pmom[ltrk];
    }
  }

  // Easy so far, but now it gets ugly: fill dp_dpar with the derivatives of the position and
  // momentum with respect to the parameters

  int       ctvmft_dim = 3 * (_ctvmq.nvertx + _ctvmq.ntrack);
  HepMatrix deriv(ctvmft_dim, 7, 0);
  //HepMatrix dp_dpar[_maxvtx] = { deriv, deriv, deriv };
  HepMatrix dp_dpar[_maxvtx] = { deriv };

  // Fill the x, y, z rows: 
  for (int nvtx = 0; nvtx < _ctvmq.nvertx; nvtx++) {
    for (int lcomp = 0; lcomp < 3; lcomp++){
      dp_dpar[nvtx][(3 * nvtx) + lcomp][lcomp] = 1.0;
    }
    
    // Fill the px, py, pz, e rows:
    for (int lvtx = 0; lvtx < _ctvmq.nvertx; lvtx++){
      for (int ltrk = 0; ltrk < _ctvmq.ntrack; ltrk++){
        if (!_ctvmq.trkvtx[lvtx][ltrk])
	  continue;
	
	// Find the derivatives of dphi with respect to x, y, curvature, and phi0:
        double dphi_dx = 2.0 * _ctvmq.par[ltrk][0] * cosph[ltrk]/cosdph[ltrk];
        double dphi_dy = 2.0 * _ctvmq.par[ltrk][0] * sinph[ltrk]/cosdph[ltrk];
        double dphi_dc = 2.0 * 
          (_ctvmq.xyzvrt[lvtx + 1][0] * cosph[ltrk] + 
           _ctvmq.xyzvrt[lvtx + 1][1] * sinph[ltrk])/cosdph[ltrk];
        double dphi_dp = 2.0 * _ctvmq.par[ltrk][0] *
          (-_ctvmq.xyzvrt[lvtx + 1][0] * sinph[ltrk] + 
	    _ctvmq.xyzvrt[lvtx + 1][1] * cosph[ltrk])/cosdph[ltrk];

	// Now load the derivative matrix
        int lvele = 3 * lvtx;
	// dPx/dx:
        dp_dpar[nvtx][lvele][3] += -pmom[ltrk].y() * dphi_dx;
	// dPy/dx:
        dp_dpar[nvtx][lvele][4] +=  pmom[ltrk].x() * dphi_dx;
	// dPz/dx:
        dp_dpar[nvtx][lvele][5]  = 0.; 
	// dPy/dx:
        dp_dpar[nvtx][lvele][6]  = 0.; 
	
        lvele++;
	// dPx/dy:
        dp_dpar[nvtx][lvele][3] += -pmom[ltrk].y() * dphi_dy;
	// dPy/dy:
        dp_dpar[nvtx][lvele][4] +=  pmom[ltrk].x() * dphi_dy;
	// dPz/dy:
        dp_dpar[nvtx][lvele][5]  = 0.; 
	// dE/dy:
        dp_dpar[nvtx][lvele][6]  = 0.; 
	
        lvele++;
	// dPx/dz:
        dp_dpar[nvtx][lvele][3]  = 0.;
	// dPy/dz:
        dp_dpar[nvtx][lvele][4]  = 0.;
	// dPz/dz:
        dp_dpar[nvtx][lvele][5]  = 0.; 
	// dE/dz:
        dp_dpar[nvtx][lvele][6]  = 0.; 
	
        lvele = 3 * (ltrk + _ctvmq.nvertx);
	// dPx/dcurv[ltrk]:
        dp_dpar[nvtx][lvele][3]  = -(pmom[ltrk].x()/_ctvmq.par[ltrk][0])
	  - pmom[ltrk].y() * dphi_dc;
	// dPy/dcurv[ltrk]:
        dp_dpar[nvtx][lvele][4]  = -(pmom[ltrk].y()/_ctvmq.par[ltrk][0])
	  + pmom[ltrk].x() * dphi_dc;
	// dPz/dcurv[ltrk]:
        dp_dpar[nvtx][lvele][5]  = -(pmom[ltrk].z()/_ctvmq.par[ltrk][0]); 
	// dE/dcurv[ltrk]:
        dp_dpar[nvtx][lvele][6]  = 
          -mom3[ltrk].mag2()/(_ctvmq.par[ltrk][0] * pmom[ltrk].e()); 
	
        lvele++;
	// dPx/dphi[ltrk]:
        dp_dpar[nvtx][lvele][3]  = -pmom[ltrk].y() * (1.0 + dphi_dp);
	// dPy/dphi[ltrk]:
        dp_dpar[nvtx][lvele][4]  =  pmom[ltrk].x() * (1.0 + dphi_dp);
	// dPz/dphi[ltrk]:
        dp_dpar[nvtx][lvele][5]  = 0.;
	// dE/dphi[ltrk]:
        dp_dpar[nvtx][lvele][6]  = 0.;
	
        lvele++;
	// dPx/dcot[ltrk]:
        dp_dpar[nvtx][lvele][3]  = 0;
	// dPy/dcot[ltrk]:
        dp_dpar[nvtx][lvele][4]  = 0;
	// dPz/dcot[ltrk]:
        dp_dpar[nvtx][lvele][5]  = pmom[ltrk].perp();
	// dE/dcot[ltrk]:
        dp_dpar[nvtx][lvele][6]  = 
	  pmom[ltrk].perp2() * _ctvmq.par[ltrk][2] / pmom[ltrk].e();
      }
    }
  }

  // -----------------------------------------------------------------------------------------------
  // Now calculate the new error matrix
  // -----------------------------------------------------------------------------------------------
  // Extract the interesting bits from VMAT
  HepMatrix vmat(ctvmft_dim,ctvmft_dim,0);
  for (int lpar = 0; lpar < ctvmft_dim; lpar++){
    int l = lpar%3;
    int lvele = 0;
    if (lpar < 3  * _ctvmq.nvertx){
      int lvtx = lpar/3;
      lvele = _ctvmq.voff[lvtx] + l;
    }
    else {
      int ltrk = (lpar - 3 * _ctvmq.nvertx)/3;
      lvele = _ctvmq.toff[ltrk] + l;
    }
    for (int kpar = 0; kpar < ctvmft_dim; kpar++) {
      int k = kpar%3;
      int kvele = 0;
      if (kpar < 3  * _ctvmq.nvertx) {
        int kvtx = kpar/3;
        kvele = _ctvmq.voff[kvtx] + k;
      }
      else{
        int ktrk = (kpar - 3 * _ctvmq.nvertx)/3;
        kvele = _ctvmq.toff[ktrk] + k;
      }
      vmat[kpar][lpar] = _ctvmfr.vmat[kvele][lvele];
    }
  }

  //  Just about done
  HepMatrix ans(7,7,0);
  //HepMatrix answer[_maxvtx] = { ans, ans, ans };
  HepMatrix answer[_maxvtx] = { ans };
  for (int nvtx = 0; nvtx < _ctvmq.nvertx; nvtx++)
    answer[nvtx] = (dp_dpar[nvtx].T() * vmat) * dp_dpar[nvtx];
  errors = answer[0];
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::vOff(vertexNumber jv) const
{
  if (jv < VERTEX_1 || jv > _maxvtx)
    return -999;
  else
    return _ctvmq.voff[jv-1];
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::tOff(const int trkId) const
{
  for (int kt=0; kt<_ctvmq.ntrack; ++kt) {
    if (trkId == _ctvmq.list[kt])
      return _ctvmq.toff[kt];
  }
  return -999;
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::pOff(int lp) const
{ 
  if (lp < 1)
    return -999;
  else
    return _ctvmq.poff[lp-1];
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::cOff(int lc) const
{
  if (lc < 1)
    return -999;
  else
    return _ctvmq.coff[lc-1];
}

//--------------------------------------------------------------------------------------------------
int MultiVertexFitter::mOff() const
{
  return _ctvmq.moff;
}

//--------------------------------------------------------------------------------------------------
double MultiVertexFitter::VMat(int i, int j) const
{
  if (i <0 || j < 0)
    return -999.;
  else
    return _ctvmfr.vmat[i][j];
}

// Facilitates CandNode recursion.  CandNodes have no way of deciding which vertex they are, and
// these trivial functions help them do that.
//--------------------------------------------------------------------------------------------------
MultiVertexFitter::vertexNumber MultiVertexFitter::allocateVertexNumber()
{
  if ((_currentAllocatedVertexNumber < PRIMARY_VERTEX) ||
      (_currentAllocatedVertexNumber > _maxvtx)) {
    cout << "MultiVertexFitter::allocateVertexNumber: out of range!" << endl;
    return PRIMARY_VERTEX;
  }
  return vertexNumber(++_currentAllocatedVertexNumber);
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::resetAllocatedVertexNumber()
{
  _currentAllocatedVertexNumber = 0;
}

//--------------------------------------------------------------------------------------------------
void MultiVertexFitter::restoreFromCommons()
{
  _stat       = 0;
  // _ctvmq_com  = (CTVMQ*)  ctvmq_address_();
  // _ctvmfr_com = (CTVMFR*) ctvmfr_address_();
  // _fiddle_com = (FIDDLE*) fiddle_address_();
  // _trkprm_com = (TRKPRM*) trkprm_address_();
  _ctvmq_com = (CTVMQ*)ctvmfr_;
  _ctvmfr_com = (CTVMFR*)ctvmfd_;
  _fiddle_com = (FIDDLE*)dctvmfi_;
  _trkprm_com = (TRKPRM*)dctvmtr_;
  _ctvmq      = *_ctvmq_com;
  _ctvmfr     = *_ctvmfr_com;
  _fiddle     = *_fiddle_com;
  _trkprm     = *_trkprm_com;
}
