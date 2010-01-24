c $Id: cctvmft.F,v 1.1 2008/09/17 04:01:49 loizides Exp $
#include "MitCommon/Ctvmft/interface/cdimensions.hh"

C DEC/CMS REPLACEMENT HISTORY, Element CTVMFT.CDF
C
C       16-Feb-2004 12:00:00  PAUS  use .h files and synchronize dimensions
C       05-JUN-2003 12:00:00  PAUS  allow for multiple single track pointing
C       09-JAN-2003 14:30  F. BEDESCHI  Introduce beam line constrain option
C       15-SEP-2002 17:50  J. GUIMARAES Increased maximum number of tracks
C                                       to 50 (MaxTrk)
C
C       09-AUG-2002 14:13  msmartin  : Zero track vertexing supported from
C       the work of Jonas Radamacker, Mat Martin and Bill Ashmanskas.
C *36   27-JUN-1995 16:30:55 BERGE "ooops???  do what I said I wanted to do. How
C       did I fail?"
C *35   27-JUN-1995 15:35:32 BERGE "improve the aesthetic aspects of I*2 handlin
C      g"
C *34   24-JUN-1995 16:59:17 LAMMEL "correct bad FORMAT statments and I*2 intrin
C      sic functions"
C *33   27-MAY-1995 16:56:37 BERGE "protect against negative mass conversion fit
C       results"
C *32   11-MAY-1995 13:43:30 BERGE "fix IJK name"
C *31    9-MAY-1995 12:02:02 BERGE "fix first approximation, STV's"
C *30    6-APR-1995 12:29:06 BERGE "minus sign screwup"
C *29    5-APR-1995 14:00:06 BERGE "bug fixes, error handling rationalization, D
C      CALC cahnges, debug dump formatting, etcetra"
C *28   28-FEB-1995 10:44:28 BERGE "bug fix in vertex/track initialization"
C *27   10-FEB-1995 21:13:26 BERGE "modify DCALC to return both Lxy and Lxyz"
C *26   12-AUG-1994 13:31:33 BERGE "fix bug in CTVMFA call of GETTRK"
C *25   10-FEB-1994 19:08:21 BERGE "fix setting of EXCUSE ME"
C *24   28-JAN-1994 14:10:28 BERGE "modify for D_KLINE debugging"
C *23   20-DEC-1993 10:41:29 BERGE "modify error handling, improve print routine
C      "
C *22    5-NOV-1993 14:52:42 BERGE "minor bug fixes, idiot-proofing, format impr
C      ovements, = ?improvements?"
C *21   15-SEP-1993 14:02:04 BERGE "another.  it's getting to be a bad habit."
C *20   15-SEP-1993 13:51:47 BERGE "TRIV"
C *19   14-SEP-1993 14:25:56 BERGE "TRIV"
C *18   13-SEP-1993 14:05:24 BERGE "assuage Stephan's worries."
C *17   12-SEP-1993 12:18:10 BERGE "buffer ERLOGR calls left in DCALC"
C *16   12-SEP-1993 11:11:55 BERGE "Bells and whistles.  Idiot-proofing."
C *15    9-SEP-1993 16:53:46 BERGE "changes to convergence testing, format fixes
C      "
C *14    8-SEP-1993 17:28:12 BERGE "new(er) include file, convergence testing/li
C      mits"
C *13    5-SEP-1993 18:57:50 BERGE "Modify to multi-vertex version"
C *12   26-JUL-1993 14:02:45 BERGE "fix GETTRK calling sequence"
C *11    6-JUL-1993 14:01:49 BERGE "prepare for new GETTRK calling sequence"
C *10    9-JUN-1993 16:49:57 BERGE "correct comment on bank name specification"
C *9    12-OCT-1992 15:51:02 BERGE "protect the print routine from nonsense cova
C      riances"
C *8     9-JUL-1992 16:26:35 BERGE "use GETTRK, change mass constraint conventio
C      n"
C *7    19-APR-1992 10:26:30 LAMMEL "re-re-port to rs6000: FORMAT"
C *6     3-APR-1992 18:39:51 BERGE "Fix pointing constraint, use DINV"
C *5    22-MAR-1992 02:46:43 LAMMEL "re-port to rs6000: FORMAT"
C *4    12-MAR-1992 11:10:48 BERGE "fix initialization bug (Gerry Lynch)"
C *3    25-FEB-1992 14:54:38 BERGE "bug fixes, error print suppression"
C *2    12-JAN-1992 22:24:01 LAMMEL "port to rs6000: FORMAT statement"
C *1    17-SEP-1991 13:25:19 BERGE "combined vertex/mass/pointing constrained fi
C      t"
C DEC/CMS REPLACEMENT HISTORY, Element CTVMFT.CDF
C===============================================================================
      SUBROUTINE CCTVMFT (PRINT,LEVPRT,IERR)
C===============================================================================
C  Authors:  John Marriner & J.P. Berge
C  Original Creation:  22 Oct   4004 B.C.E.  J.P. Marriner
C           revision:   1 Sep   2746 A.U.C.  J.P. Berge      Multiple vertices
C           revision:  10 Rabi` 1415 A.H.    W.J. Ashmanskas Conversion fitting,
C                                                            1-track vertex
C           revision:  May 2002              M.S. Martin 0-track vertex
C  CDF multiple vertex and mass constrained fit package.
C  This package is documented in CDF/DOC/SEC_VTX/PUBLIC/1996
C          (which needs revision)
C  Allows the specification of tracks at one or more vertices, and performs
C  (geometrical) fits such that at each vertex all tracks at that vertex are
C  constrained to pass through the same space point ("Vertex fits").
C  There is a special case of the vertex fit that is designed to treat vertices
C  with two tracks where the two tracks originate from a photon conversion.
C  In this case, the two tracks have such a small space opening angle that
C  the usual vertex fit becomes unstable. For a vertex flagged as a conversion
C  candidate, extra constraints are imposed. The r-phi opening at the point
C  where they meet in space. The further requirement that the r-z opening angle
C  is zero may be optionally imposed.
C  These fits are not optional.  They are always performed.
C  Allows constraining sub-groups of tracks to specified masses ("mass fit").
C  Tracks in such a mass fit may be attached to different vertices (consider
C  the case of the decay Bd -> (psj)(k0). where the Bd mass is constra=ined).
C  These fits are optional.
C  Allows constraints such that the momentum sum of tracks at a vertex "points"
C  to another vertex ("pointing fits").
C  "Target" vertices may be treated as "primary";  In these cases, vertex
C  coordinates and covariances are required input and are treated as
C  measurements in the fit. Such vertices have no attached tracks.
C  Alternatively, target vertices may be one of the vertices formed from
C  track intersection fitting in the current fit.  Such vertices are considered
C  "secondary".  A secondary vertex that is pointed at a secondary vertex is
C  a tertiary vertex.  There is allowed a special case of a "one-track" vertex
C  (consider decays Xi -> (lambda)(pi) or B* -> (D0)(pi)).  Here the fit
C  constrains the momentum vector from the tertiary vertex to intersect the
C  single track at the secondary vertex.
C  Note that "conversion vertices may point to parent vertices. Note that the
C  two tracks at such a vertex may appear in a multi-vertex mass-constrained
C  fit (consider the case Chi -> (psj)(gamma->ee)), but may not appear in any
C  mass constraint fit restricted only to this vertex.
C  These fits are optional.
C  NOTA BENE:  Considerable care has been exercized in designing a formatted
C  dump of the CTVMFT fit results. This dump is activated through setting the
C  imput parameter PRINT to some non-zero value. If this effort has achieved
C  its purpose, many points which may remain unclear to the first time user
C  should be clarified through inspection of dumps of a few fits of interesting
C  fits.  This formattewd output should be printed with 124 columns.
C
C  Additional beam constraint option
C  =================================
C  F. Bedeschi, January 2003
C
C  An option to use a beamline constraint in the fit of the primary vertex
C  has been added. This option does not make sense in the case of secondary
C  vertex fits and has to be turned on only when fitting the primary with
C  no additional constraints. The routine is protected and will not function
C  if those options are selected as well as the bea line constraint.
C  In case this option is chosen the fit will calculate a result also with
C  just one track in input.
C  In order to enable this option the user must do the following:
C  - set NVERTX = 1
C  - set Mv=VtxPnt(1,1) = -100 (no pointing constraints when <0. -100 is a new code to
C                               enable this type of fit.)
C                         NB. In this case PRIMARY should be .FALSE. in the code
C  - set Cvtx(1) = 0 (no mass constraints)
C  - provide the beamline parameters:
C    --> assume xv = xb + axb*z, yv = yb + ayb*z, where (xv, yv) is the transverse
C        beam position at z, (xb, yb) is the transverse beam position at z = 0
C        as returned by the beam line fit and axb, ayb are the slopes of the beam
C        line as returned by the beam line fit
C    --> set XYZPV0 to (xb, yb, 0)
C    --> set EXYZPV(1,1) = sigma**2 of beam spread in x (~30 - 50 um)**2
C            EXYZPV(2,2) = sigma**2 of beam spread in y (~30 - 50 um)**2
C            EXYZPV(3,3) = sigma**2 of beam spread in z (~30 cm)**2
C            all other elements of EXYZPV should be set to 0
C    --> set XZSLOPE = axb (this is new parameter in the COMMON/CTVMQ/
C            YZSLOPE = ayb (this is new parameter in the COMMON/CTVMQ/
C
C=== End beamline constraint description
C
C--Input parameters
C  PRINT   If non-zero, dumps formatted fit results to unit PRINT
C  LEVPRT  Controls the level of (optional) printing. Inoperative if PRINT=0.
C          Avoid this (other than 0,1) unless you know what you are doing.
C--Input data communicated through the include file CTVMFT (COMMON /CTVMFC/)
C  Gross fit specifications (INTEGER);
C  NTRACK  Number of tracks in this fit.  (must be at least 2)
C  NVERTX  Number of vertices in this fit (at least one is required)
C  NMASSC  Number of mass constraints     (may be zero)
C  Required topological and constraint information;
C  TrkVtx(MaxTrk,MaxVtx)    LOGICAL
C          Specifies vertex to which each track is associated.
C          TrkVtx(Nt,Nv) is .TRUE. if track Nt is attached to vertex Nv.
C  VtxPnt(MaxVtx,2)         INTEGER
C         Vertex association fit specification.
C         If Mv=VtxPnt(Nv,1) < 0, no pointing fit is done.
C         If Mv=VtxPnt(1,1) = -100 a primary vertex fit with beamline constraint
C                             is assumed.
C         If Mv=VtxPnt(Nv,1) (Nv>Mv) is .GE.0, vertex Nv "points" to vertex Mv.
C         If Mv=0, the target vertex Mv is supplied as a measurement; the vertex
C         (x,y,z) together with its 3x3 covariance matrix must be furnished by
C         the user through the REAL arrays XYZPV0,EXYZPV.  "Primary" vertices
C         are the most usual type of such vertices.
C         If Mv > 0, the target vertex Mv is a secondary vertex and its coord-
C         inates are parameters resulting from a track intersection fit.
C         The kind of "pointing" of vertex Nv to vertex Mv is specified by
C         setting VtxPnt(Nv,2) to 1, 2, 3 or 4 as desired.
C         Note that more than one vertex can "point" to the same target vertex.
C         The constraint is that the sum of the momenta of all the tracks at
C         vertex Nv together with the summed momenta of all tracks at vertices
C         which point to Nv must be parallel OR anti-parallel to the vector
C         from the vertex Nv to vertex Mv. Allowed values of VtxPnt(Nv,2) are;
C           1 = constrain pointing of Nv towards Mv in the (R,Phi) plane.
C           2 = ALSO constrain pointing in (R,Z) plane (3 dimensional pointing).
C           3 = point Nv toward single-track vertex Mv
C           4 = point Nv toward zero-track vertex Mv
C         Note that vertices connected via a charged track are not correctly
C         handled.
C  Cvtx(MaxVtx)             INTEGER
C          Conversion constraint fit specification.  (The conversion fit
C          constrains the two tracks at a vertex to be parallel at the vertex.)
C           0 = no conversion fit is performed, the vertex is treated normally.
C           1 = constrain the track directions in the r-phi plane to be the same
C               at the radius of the vertex (at the point where the tracks are
C               copunctual, they are clinear).
C           2 = ALSO constrain opening angle to zero in r-z.
C  TrkMcn(MaxTrk,MaxMcn)    LOGICAL
C          Specifies the tracks that participate in each mass constraint.
C          Track Nt participates in constraint Nm if TrkMcn(Nt,Nm) is .TRUE.
C  CMASS(MaxMcn)            REAL
C          The Nmth cluster of tracks is constrained to have mass CMASS(Nm).
C  Track bank and data specifications;
C  TKBANK(MaxTrk)           CHARACTER*4
C         Bank name (exempli gratia, 'QTRK','TRKS','SVXS','CTCS') from which
C         to get each track's helix parameter and covariance data.
C         These bank names are used in any substantive way only in the data
C         access subroutine.  In CDF, this routine is C$TRK:GETTRK.
C  LIST  (MaxTrk)           INTEGER
C         Bank numbers of the tracks to be used.
C  TMASS (MaxTrk)           REAL
C         Masses (GeV/c**2) to be associated with the tracks, for purposes of
C         dE/dX and multiple Coulomb scattering calculations, and for any mass-
C         constrained fitting.
C  Output parameter
C  IERR   =  0 (No error), or error code
C  Returned through the common /CTVMFC/
C  PAR(5,NTRACK)      = The track parameter fit results
C  XYZVRT(3,0:MaxVtx)   The vertex x,y,z fit results.
C  VMAT(Maxdim,Maxdim)  The covariance matrix on the fit parameters
C  FMCDIF(NMASSC)       Residuals of the mass constraints (if used)
C  PCON(MaxVtx,1)       Residual of the (R,Phi) pointing constraint (if used)
C  PCON(MaxVtx,2)       Residual of the (R,Z)   pointing constraint (if used)
C  SANG(MaxVtx,1)       sin(ang) between vertex direction and momentum in x-y
C  SANG(MaxVtx,2)       sin(ang) between vertex direction and momentum in rho-z
C  CHISQR(0:MaxItr)     fit Chi Square result(s)
C  ... other goodies ...  Go look at the include file.
C Convenient access to VMAT entries for the fit results for vertex and track
C parameters is provided by the offset pointer lists VOFF and TOFF.
C VOFF(Nv+i) (i=1:3) gives row/column indices for vertex Nv (x,y,z).
C TOFF(Nt+i) (i=1:3) gives row/column indices for track Nt (Crv,Phi,Ctg)
C  Subroutines;
C  CTVM00  Checks the input for completeness and consistency
C  CTVMFA  Makes the vertex first approximations (uses GETTRK).
C  CTVM01  Collects the track data (uses GETTRK).
C  CTVMVF  Calculates the derivative contributions for vertex fitting.
C  CTVMCF                                              conversion fitting.
C  CTVMPF                                              pointing fittting.
C  CTVMMF                                              mass constraint fitting.
C  CTVMPR  Prints the fit results, if desired, to output unit PRINT.
C-------------------------------------------------------------------------------
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
#include "MitCommon/Ctvmft/interface/cctvmfi.h"
#include "MitCommon/Ctvmft/interface/cctvmuv.h"
#include "MitCommon/Ctvmft/interface/cctvmtu.h"
      INTEGER PRINT, IERR
      LOGICAL PRIMARY
      REAL*8  VXI(MAXDIM)
      REAL*8  WXYZPV(3,3), SUM
      REAL*8  LMAT(3,3), LTV(3,3), LTVL(3,3)
      REAL*8  LTVLX0(3), LTVXV(3), DXBEAM(3)
      REAL    C,PHI,COTH, D,Z, PT, RADV, TWOC, S,TS,TRAD, T1,T2,T3,T4
      REAL    SP,CP, XCYS,XSYC,XC2YC2, SINCS,COSCS,TWOCS, PMAG, XMAG
      REAL    XSVT,YSVT, WORK(MAXDIM), ADDCHIV
      REAL    PARFAC,CUTFAC, DELVS(3,MaxVtx)
      REAL    PCSAVE(MaxVtx,2), MFSAVE(MaxMcn), MFTEST(MaxMcn)
      REAL    STPTK,STEP(5), PARDJF(5,MaxTrk)
      INTEGER Nv,Mv,Nt,Np,Nm, NvF,NtF,LpF,LmF,LcF
      INTEGER NITER, NMAXCT,ISCUT,NPC
      INTEGER I,J,K,L, II(2), I1,I2, IFAIL, IMAX, VERTEX,VRTX
      INTEGER KTRACK,NPOINT,NCONV, KPRINT
C  Local variables pertaining to single-track vertices
      INTEGER STVCOUNT,STVFLAG(MaxVtx)
      INTEGER BEAM
      LOGICAL FIRST
      INTEGER MRUN,MNEV
      CHARACTER*6  NAME
      CHARACTER*80 STRING
C LEVPRT = print level (0=print disabled)
      INTEGER LEVPRT
      REAL    ZERO,MINUS,SHOUT, TEST
      SAVE    FIRST,BEAM
      SAVE    NMAXCT, MRUN,MNEV
      SAVE    ZERO,MINUS
c      DATA    ME   /0/
      DATA    ZERO /0.0/, MINUS /-1.0/
C Maximum number of cut steps allowed in any step
      DATA    NMAXCT / 7 /
      DATA    FIRST /.TRUE./
      DATA    MRUN  /-1/, MNEV /-1/

CP>>
c     For tuning only
      LOGICAL LOST
      INTEGER NTOT,NERR,NERRC2,NERRLX,NERRLN
      INTEGER NFIT(0:MAXITR),NOPT(0:MAXITR)
      REAL    MAGN,DV(3),PTVTX,PTTK1,PTTK2,PTTK3
      REAL    DXYZ(3),Dr,Dz,Dl(3)
      REAL    CHI2(0:MAXITR),CRP2(0:MAXITR),LXY(0:MAXITR)

      INTEGER IXCRP2(3)
      DATA    IXCRP2 / 1, 2, 4 /
CP<<
C*******************************************************************************
      DBGPRT = IABS(LEVPRT)
C avoid terminating on error in CTVM00
C Comment out the following line so that instead of bombing on an error,
C CTVMFT returns to VertexFit, so that the error can be properly handled.
C Craig Blocker  July 19, 2001
C      EXCUSE = ME
      IF (LEVPRT.LT.0) THEN
        EXCUSE =-1
      END IF
      IF (PRINT.EQ.0) DBGPRT=0
      KPRINT = 0
      IF (PRINT.GT.0 .AND. DBGPRT.GT.10) KPRINT=PRINT
C First entrance initializations
      IF (FIRST) THEN
        FIRST =.FALSE.
C do NOT use beam constrained parameters
        BEAM  = 0
        WRITE(print,'(A)') ' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        WRITE(print,'(A)') ' ---- CCTVMFT First call in program -------'
        WRITE(print,'(A)') ' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
        WRITE(print,'(3(A9,I4),/,3(A9,I4))')
     +    '  VtxsMax:',MAXVTX,'  TrksMax:',MAXTRK,'  ConsMax:',MAXMCN,
     +    '  DimsMax:',MAXDIM,'  ItrsMax:',MAXITR,'  UDim:   ',UDIM
        WRITE(print,'(A)') ' XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX'
      END IF
      IF (MRUN.NE.RUNNUM .OR. MNEV.NE.TRGNUM) THEN
C New event, reset track "memory"
        DO I=1,MaxTrk
          TKERR(I) =-1.0
        END DO
        MRUN = RUNNUM
        MNEV = TRGNUM
      END IF
      IERR   = 0
      DO I=1,3
        IJKERR(I) = 0
      END DO
      NTSCUT = 0
      CALL VZERO(CHISQR,MaxItr+1)
      CHISQR(0) =-1.
      CALL VZERO(CHIV  ,MaxVtx+1)
      CALL VZERO(CHIT  ,MaxTrk)
      CALL VZERO(PARDIF,MaxTrk*5)
      CALL VZERO(PARDJF,MaxTrk*5)
      CALL VZERO(PARERR,MaxTrk*5)
      CALL VZERO(TrkP4, MaxTrk*6)
      CALL VZERO(DDA,   MaxTrk*8)
      CALL VZERO(XYZVRT,MaxVtx*3+3)
      CALL VZERO(DELVS, MaxVtx*3)
      CALL VZERO(PCSAVE,MaxVtx*2)
      CALL VZERO(LMAT,  2*3*3)
      CALL VZERO(LTV,   2*3*3)
      CALL VZERO(LTVL,  2*3*3)
      CALL VZERO(LTVLX0,2*3)
      CALL VZERO(LTVXV, 2*3)
      CALL VZERO(MFSAVE,MaxMcn)
      CALL VZERO(MFTEST,MaxMcn)
      CALL VZERO(STVFLAG,MaxVtx)
      STVCOUNT = 0
C  Check Input Conditions for illegal conditions
C------------------------------------------------------------------------------
C Check input specifications
      CALL CCTVM00 (IERR,WXYZPV)
      IF (IERR.NE.0) THEN
        RETURN
      END IF
      TRAD = ABS(XYZPV0(1)) + ABS(XYZPV0(2))
C flag the appearance of a primary vertex
      PRIMARY =.FALSE.
C check for a primary vertex
      DO Nv=1,NVERTX
        IF (VtxPnt(Nv,1).EQ.0) THEN
          PRIMARY =.TRUE.
        END IF
      END DO
C primary vertex coordinates - Joao 02/11/03
      IF (PRIMARY .OR. (TRAD.GT.0.0.AND.VTXPNT(1,1).NE.-100)) THEN
         DO I=1,3
            XYZVRT(I,0) = XYZPV0(I)
C primary vertex fit displacement
            DXYZPV(I) = 0.0
         END DO
      ENDIF
C beamline constraint initializations
      IF (VTXPNT(1,1) .EQ.-100) THEN
         DO I = 1,3
            LMAT(I,I) = 1.0D0
         ENDDO
         LMAT(1,3) = -XZSLOPE
         LMAT(2,3) = -YZSLOPE
C Calculate Lt*V**-1
         CALL VZERO(LTV,2*3*3)
         DO I = 1,3
            DO J = 1,3
               DO K = 1,3
                  LTV(I,J) = LTV(I,J) + LMAT(K,I)*WXYZPV(K,J)
               ENDDO
            ENDDO
         ENDDO
C Calculate Lt*V**-1*L
         CALL VZERO(LTVL,2*3*3)
         DO I = 1,3
            DO J = 1,3
               DO K = 1,3
                  LTVL(I,J) = LTVL(I,J) + LTV(I,K)*LMAT(K,J)
               ENDDO
            ENDDO
         ENDDO
D        IF (DBGPRT.GT.0) THEN
D        write(print,5110)
D5110    format(1x,'CTVMFT: Beamline initialization done')
D         write(print,5120)XYZPV0,XZSLOPE, YZSLOPE
D5120    format(1x,'Starting vertex x,y,z: ',3(E10.4,2x),
D    >        'Slopes: x,y: ',2(e10.4,2x))
D         do i=1,3
D            write(print,5130)(WXYZPV(i,j), j=1,3)
D         enddo
D5130    format(1x,'WXYZPV: ',3(e10.4,2x))
D         do i = 1,3
D            write(print,5140)(ltv(i,j),j=1,3),(ltvl(i,j),j=1,3)
D5140       format(1x,'LTV: ',3(E10.4,2x)'LTVL: ',3(e10.4,2x))
D         enddo
D         ENDIF
      ENDIF

D     IF (DBGPRT.GT.0) THEN
D     WRITE(PRINT,1000) RUNNUM,TRGNUM, NTRACK,NVERTX,NMASSC
D     IF (PRIMARY) THEN
D       WRITE(PRINT,1005) XYZPV0
D       DO I=1,3
D         WORK(I) = SQRT(EXYZPV(I,I))
D       END DO
D       DO I=1,3
D       DO J=1,3
D         VMAT(I,J) = EXYZPV(I,J)/(WORK(I)*WORK(J))
D       END DO
D       END DO
D       WRITE(PRINT,1006) (WORK(I),I=1,3)
D       WRITE(PRINT,1006) ((VMAT(I,J),I=1,3),J=1,3)
D     END IF
D     WRITE(PRINT,1001) (I, LIST(I),TkBank(I),Tmass(I),I=1,NTRACK)
D1000 FORMAT('1Event ',I5,'.',I6,';   CTVMFT fit results,'
D    &,'        (Tracks,Vertices,Mass Constraints) ',2(I2,','),I2, /)
D1001 FORMAT('  Track',I3,';',I5,A5,F8.3 )
D1005 FORMAT(  ' Primary vertex ',3F10.4)
D1006 FORMAT(16X,3F10.4)
D     END IF
C         Vertex fit matrix dimensions;
C       3 vertex parameters (x,y,z) per (possible) primary vertex
C       3 vertex parameters per secondary vertex
C     + 0, 1, or 2 "pointing" parameters
C     + 3 parameters (curvature,phi0,coth(theta0)) per track
C     + N mass constraints
C  VOFF(Nv),TOFF(Nt) are offset pointers to ((Vertex,Pointing),Track) parameters
C           in (the 1st derivative vector, 2nd derivative matrix)  VXI, VMAT
C  VOFF(Nv) points to the start of the 3, 4, or 5 variables for vertex Nv;
C                                                    vertex Nv (x,y,z, p1,p2)
C  TOFF(Nt) points to the start of the 3 parameters; track  Nt (Crv,Phi,Ctg)
C  POFF(Lp) points to the start of the Lagrangian multiplier values for the
C           pointing constraint Lp.
C  COFF(Lc) points to the start of the Lagrangian multiplier values for the
C           conversion constraint Lc.
C  MOFF     points to the start of the Lagrangian multiplier values for the
C           mass constraints
C initialize the parameter offset pointer
      NvF = 0
      IF (PRIMARY) NvF = 3
C offsets to vertex parameters
      DO Nv=1,NVERTX
        VOFF(Nv) = NvF
        NvF = NvF + 3
      END DO
      NtF = NvF
C offsets to track parameters
      DO Nt=1,NTRACK
        TOFF(Nt) = NtF
        NtF = NtF + 3
      END DO
C count the number of "pointing" constraints
      NPOINT = 0
      LpF = NtF
C offsets to pointing Lagrangian multipliers
      DO Nv=1,NVERTX
        IF (VtxPnt(Nv,1).GE.0) THEN
          POFF(Nv) = LpF
          IF (VtxPnt(Nv,2).EQ.1.OR.VtxPnt(Nv,2).EQ.2) THEN
            NPOINT = NPOINT + VtxPnt(Nv,2)
            LpF    = LpF    + VtxPnt(Nv,2)
C "one-track vertex"
          ELSE IF (VtxPnt(Nv,2).EQ.3) THEN
            IF (STVFLAG(VtxPnt(Nv,1)).EQ.0) THEN
C  Count number of unique single-track vertices
              STVFLAG(VtxPnt(Nv,1)) = 1
              STVCOUNT = STVCOUNT+1
            ENDIF
            NPOINT = NPOINT + 2
            LpF    = LpF    + 2
C  "zero-track vertex"
          ELSE IF (VtxPnt(Nv,2).EQ.4) THEN
            NPOINT = NPOINT + 2
            LpF    = LpF    + 2
          END IF
        END IF
      END DO
C count number of "conversion" constraints
      NCONV = 0
      LcF = LpF
C offsets to zero-opening-angle constraint
      DO Nv=1,NVERTX
C Lagrangian multipliers
        IF (Cvtx(Nv).GT.0) THEN
          COFF(Nv) = LcF
          NCONV    = NCONV + Cvtx(Nv)
          LcF      = LcF   + Cvtx(Nv)
        END IF
      END DO
C offset to mass cnst Lagrangian multipliers
      MOFF = LcF
C Fit dimension (number of degrees of freedom)
      NDOF = 2*NTRACK - 3*NVERTX + NPOINT+NMASSC+NCONV
      IF (VTXPNT(1,1) .eq. -100)NDOF = NDOF + 3
C Dimension of matrix equations
      MATDIM = NtF + NMASSC+NPOINT+NCONV
C Data Collection and Checking, Fit Initialization, Vertex first approximation
C------------------------------------------------------------------------------
C Loop (backward) over secondary vertices. Ignore for primary vertex fits with
C beamline constraint.
      DO Nv=NVERTX,1,-1
         IERR = 0
         CALL CCTVMFA (Nv, PRINT,IERR)
         IF (IERR.NE.0) THEN
            WRITE (STRING,2000) IJKERR
 2000       FORMAT(I4,2I3,' vertex first approximation failure')
            GO TO 990
         END IF
      END DO
C------------------------------------------------------------------------------
C-------------- Fit iteration start -------------------------------------------
      ITER  = 0
      NITER = MAXITR
C fit iteration (re)entrance point
  100 CONTINUE
D     IF (DBGPRT.GT.0) THEN
D     WRITE(PRINT,1020) ITER
D1020 FORMAT(/,' Start Iteration step',I2)
D     END IF
C (re)initialize VXI, VMAT
      CALL VZERO(VXI, 2*MAXDIM)
      CALL VZERO(VMAT,2*MAXDIM*(MAXDIM+1))
C stuff the primary Vtx covariance into VMAT
      IF (PRIMARY) THEN
        DO I=1,3
        DO J=1,3
          VMAT(I,J) = WXYZPV(I,J)
C (measured - this fit)
          VXI(I) = VXI(I) - VMAT(I,J)*DXYZPV(J)
        END DO
        END DO
      END IF
C  Caclulate terms for the vertex fits via track inentersection.
      DO Nv=NVERTX,1,-1
        CALL CCTVMVF (KPRINT, Nv, VXI)
      END DO
C  accumulate vertex momentum sums, for (possible) pointing constraint fitting
C  4-momentum sums for vertex pointing
      CALL VZERO(VtxP4,MaxVtx*4)
      DO Nv=NVERTX,1,-1
        DO Nt=1,NTRACK
        IF (TrkVtx(Nt,Nv)) THEN
          DO I=1,4
            VtxP4(I,Nv) = VtxP4(I,Nv) + TrkP4(Nt,I)
          END DO
        END IF
C end track loop
        END DO
C does Nv point?
        Mv = VtxPnt(Nv,1)
C tertiary; add P(Nv) to P(Mv)
        IF (Mv.GT.0) THEN
          DO I=1,4
            VtxP4(I,Mv) = VtxP4(I,Mv) + VtxP4(I,Nv)
          END DO
        END IF
C End vertex loop
      END DO
C  calculate terms for "conversion" (zero opening angle) constraint
      DO Nv=1,NVERTX
C conversion constraint, this vertex
      IF (Cvtx(Nv).GT.0) THEN
        CALL CCTVMCF (KPRINT, Nv,  VXI)
      END IF
C end of the vertex loop
      END DO
C  calculate terms for pointing from this vertex towards its parent
      DO Nv=NVERTX,1,-1
C pointing constraint, this vertex
      IF (VtxPnt(Nv,1).GE.0) THEN
        CALL CCTVMPF (KPRINT, Nv,  VXI)
C save the last step quality
        LpF = POFF(Nv)
        PCSAVE(Nv,1) = ABS(VXI(LpF+1))
        PCSAVE(Nv,2) = ABS(VXI(LpF+2))
      END IF
C end of the vertex loop
      END DO
C  accumulate momentum sums, for (possible) mass constraint fitting
C 4-momentum sums for mass constraints
      CALL VZERO(McnP4,MaxMcn*4)
      DO Nm=1,NMASSC
      DO Nt=1,NTRACK
      IF (TrkMcn(Nt,Nm)) THEN
        DO I=1,4
          McnP4(I,Nm) = McnP4(I,Nm) + TrkP4(Nt,I)
        END DO
      END IF
      END DO
      END DO
C loop over mass constraints
      DO Nm=1,NMASSC
        CALL CCTVMMF (KPRINT, Nm, VXI)
        LmF = MOFF
        MFSAVE(Nm) = ABS(VXI(LmF+Nm))
      END DO
C
C  Additional terms for beamline constraint
C
        IF (VTXPNT(1,1) .EQ. -100) THEN
D        IF (DBGPRT.GT.0) THEN
D           write(print,5150)NTRACK, (XYZVRT(i,1),i=1,3)
D5150      format(1x,'Ntracks: ',i3,
D    >          'Vertex before iteration: ',3(e10.4,2x))
D       ENDIF
C Update VMAT
          DO I = 1,3
            DO J = 1,3
                VMAT(I,J) = VMAT(I,J) + LTVL(I,J)
            ENDDO
          ENDDO
C Update VXI
          CALL VZERO(LTVXV, 2*3)
          CALL VZERO(LTVLX0,2*3)
          DO I = 1,3
            DO J = 1,3
                LTVXV(I)  = LTVXV(I)  + LTV(I,J)*XYZPV0(J)
                LTVLX0(I) = LTVLX0(I) + LTVL(I,J)*XYZVRT(J,1)
            ENDDO
            VXI(I) = VXI(I) + LTVXV(I) - LTVLX0(I)
          ENDDO
        ENDIF
C End of the derivative accumulation phase;  Now solve matrix equations
C Find solution without finding inverse
      IF (ITER.EQ.0) THEN
        CALL MYDEQN  (MATDIM,VMAT,MAXDIM,WORK,IFAIL,1,VXI)
C Find solution + covariance matrix
      ELSE
        CALL MYDEQINV(MATDIM,VMAT,MAXDIM,WORK,IFAIL,1,VXI)
      ENDIF
      IF (IFAIL.NE.0) THEN
        IERR = 9
        IJKERR(2) = 1
        WRITE(STRING,2030)
 2030   FORMAT(' singular solution matrix')
        RETURN
      END IF
C Check the step just taken for acceptability and for possible convergence
C Initialize cut steps counter
      ISCUT = 0
C Step size scale factor - normally 1
      CUTFAC = 1.0
C Fraction of the step already taken
      PARFAC = 0.0
C   Check that track turning angles to the new vertex points are acceptable
  110 CONTINUE
C x,y for this vertex, this step cut factor
      DO 150 Nv=1,NVERTX
      NvF = VOFF(Nv)
        XSVT = XYZVRT(1,Nv) + CUTFAC*VXI(NvF+1)
        YSVT = XYZVRT(2,Nv) + CUTFAC*VXI(NvF+2)
C track loop, checking turning angle
        DO 140 Nt=1,NTRACK
        IF (.NOT.TrkVtx(Nt,Nv)) GO TO 140
        NtF = TOFF(Nt)
          C   = PAR0(1,Nt) + PARDIF(1,Nt) + CUTFAC*VXI(NtF+1)
          PHI = PAR0(2,Nt) + PARDIF(2,Nt) + CUTFAC*VXI(NtF+2)
          TS  = 2.0*C*(XSVT*COS(PHI)+YSVT*SIN(PHI))
C track turns too much, cut step size
          IF (ABS(TS).GT.TRNMAX) THEN
            ISCUT = ISCUT+1
            NTSCUT = NTSCUT + 1
C maximum number of cut steps exceeded
            IF (ISCUT.GT.NMAXCT) THEN
              IERR = 9
              IJKERR(2) = 2
              IJKERR(3) = ITER
              WRITE(STRING,2050) IJKERR
 2050         FORMAT(' too many cut steps', I5,2I3)
              GO TO 990
            END IF
            CUTFAC = 0.5*(CUTFAC-PARFAC)
D           IF (DBGPRT.GT.1) THEN
D           WRITE(PRINT,1050) ITER,ISCUT,Nv,Nt,XSVT,YSVT
D1050 FORMAT(' Iter',I2,' Cut',I2,',  Nv,Nt ',I1,I4,2F10.4)
D           END IF
            GO TO 110
          END IF
  140   CONTINUE
  150 CONTINUE
      PARFAC = PARFAC+CUTFAC
C  Provisionally acceptable step;
C  Update all the fit parameters and calculate Chi Square, this step
C  Re-initialize, recompute momentum sums with the new step results
C Update primary vertex coordinates
      IF (PRIMARY) THEN
        DO I=1,3
          DXYZPV(I)   = DXYZPV(I) + CUTFAC*VXI(I)
          XYZVRT(I,0) = XYZPV0(I) + DXYZPV(I)
        END DO
      END IF
C Contribution from the "primary" vertex
      CHIV(0) = 0.0
      IF (PRIMARY) THEN
        DO I=1,3
        DO J=1,3
          CHIV(0)= CHIV(0) + DXYZPV(I)*WXYZPV(I,J)*DXYZPV(J)
        END DO
        END DO
      END IF
      CHISQR(ITER) = CHIV(0)
C Update secondary vertex coordinates
      DO Nv=1,NVERTX
      Nvf = VOFF(Nv)
        DO I=1,3
          XYZVRT(I,Nv) = XYZVRT(I,Nv) + CUTFAC*VXI(NvF+I)
          Mv = VtxPnt(Nv,1)
C vertex separation vector
          IF (Mv.GE.0) THEN
            DELVS(I,Nv) = XYZVRT(I,Mv) - XYZVRT(I,Nv)
          END IF
        END DO
      END DO
D           IF (DBGPRT.GT.0) THEN
D           write(print,5160)NTRACK, (XYZVRT(i,1),i=1,3)
D5160      format(1x,'Ntracks: ',i3,
D    >          ', Vertex after iteration: ',3(e10.4,2x))
D           END IF
      DO 210 Nv=1,NVERTX
      CHIV(Nv) = 0.0
C Update the track parameters
      DO 210 Nt=1,NTRACK
      IF (.NOT.TrkVtx(Nt,Nv)) GO TO 210
        NtF = TOFF(Nt)
C half curvature step
        PARDIF(1,Nt) = PARDIF(1,Nt) + CUTFAC*VXI(NtF+1)
C phi step
        PARDIF(2,Nt) = PARDIF(2,Nt) + CUTFAC*VXI(NtF+2)
C cot(theta) step
        PARDIF(3,Nt) = PARDIF(3,Nt) + CUTFAC*VXI(NtF+3)
C new half curvature
        C    = PAR0(1,Nt) + PARDIF(1,Nt)
C new phi
        PHI  = PAR0(2,Nt) + PARDIF(2,Nt)
C new coth(theta)
        COTH = PAR0(3,Nt) + PARDIF(3,Nt)
        SP = SIN(PHI)
        CP = COS(PHI)
C         Xsv*cos(phi) + Ysv*sin(phi), -Xsv*sin(phi) + Ysv*cos(phi), 2C
        XCYS = XYZVRT(1,Nv)*CP + XYZVRT(2,Nv)*SP
        XSYC =-XYZVRT(1,Nv)*SP + XYZVRT(2,Nv)*CP
        TWOC = 2.0*C
C 10/23/02 Protect against TWOC begin zero.  Craig Blocker
c??        if(twoc .eq. 0.) then
        IF (ABS(TWOC) .LT. 1E-24) THEN
          S = XCYS
        else
          S = ASIN(TWOC*XCYS)/TWOC
        endif
        TS   = TWOC*XCYS
        TRAD = SQRT((1.-TS)*(1.+TS))
C Constrained d0
C 10/23/02 Protect against C being zero. Craig Blocker
c??        if(c .eq. 0.) then
        IF (ABS(C) .LT. 1E-24) THEN
          D = XSYC
        else
          D = XSYC - SIN(C*S)**2/C
        endif
C Constrained z0
        Z = XYZVRT(3,Nv) - COTH*S
C fit d0 - CTC fit d0
        PARDIF(4,Nt) = D - PAR0(4,Nt)
C fit z0 - CTC fit z0
        PARDIF(5,Nt) = Z - PAR0(5,Nt)
C 1/20/04 Protect against C being zero.  Here it is not
C clear what to set PT to.  I chose 999. GeV since this is
C a very large Pt in CDF, this is roughly a standard deviation
C from infinite, and should be clearly recognizable.  Craig Blocker
c??        if(c .eq. 0.) then
        IF (ABS(C) .LT. 1E-24) THEN
          PT = 999.
        else
          PT = PSCALE/ABS(C)
        endif
C x component of momentum
        TrkP4(Nt,1) = PT*(CP*TRAD-SP*TS)
C y component of momentum
        TrkP4(Nt,2) = PT*(SP*TRAD+CP*TS)
C z component of momentum
        TrkP4(Nt,3) = PT*COTH
C Transverse momentum
        TrkP4(Nt,5) = PT
C total momentum
        TrkP4(Nt,6) = SQRT(PT**2+TrkP4(Nt,3)**2)
C energy
        TrkP4(Nt,4) = SQRT(TrkP4(Nt,6)**2+TMASS(NT)**2)
C Calculate contribution to chi-squared for this track
        CHIT(Nt) = 0.0
C Loop over rows of the weight matrix
        DO I=1,5
C Loop over columns of the weight matrix
        DO J=1,5
        CHIT(Nt) = CHIT(Nt) + PARDIF(I,Nt)*G(J,I,Nt)*PARDIF(J,Nt)
        END DO
        END DO
        CHIV(Nv) = CHIV(Nv) + CHIT(Nt)
        CHISQR(ITER) = CHISQR(ITER) + CHIT(Nt)
C End the loop on tracks
  210 CONTINUE
C Beamline constraint contribution
      IF(VTXPNT(1,1).EQ.-100)THEN
         DO I = 1,3
            DXBEAM(I) = -XYZPV0(I)
            DO J = 1,3
               DXBEAM(I) = DXBEAM(I) + LMAT(I,J)*XYZVRT(J,1)
            ENDDO
         ENDDO
         ADDCHIV = 0.0
         DO I = 1,3
            DO J = 1,3
               ADDCHIV = ADDCHIV + DXBEAM(I)*DXBEAM(J)*WXYZPV(I,J)
            ENDDO
         ENDDO
         CHIV(NV) = CHIV(NV) + ADDCHIV
         CHISQR(ITER) = CHISQR(ITER) + ADDCHIV
      ENDIF
      CALL VZERO(VtxP4,MaxVtx*4)
      CALL VZERO(McnP4,MaxMcn*4)
      DO Nt=1,NTRACK
        DO Nv=NVERTX,1,-1
        IF (TrkVtx(Nt,Nv)) THEN
C accumulate vertex sums
        DO I=1,4
          VtxP4(I,Nv) = VtxP4(I,Nv) + TrkP4(Nt,I)
        END DO
        END IF
        END DO
C accumulate mass constraint sums
        DO Nm=1,NMASSC
        IF (TrkMcn(Nt,Nm)) THEN
        DO I=1,4
          McnP4(I,Nm) = McnP4(I,Nm) + TrkP4(Nt,I)
        END DO
        END IF
        END DO
      END DO
      DO Nv=1,NVERTX
C does Nv point?
        Mv = VtxPnt(Nv,1)
C tertiary; add P(Nv) to P(Mv)
        IF (Mv.GT.0) THEN
        DO I=1,4
          VtxP4(I,Mv) = VtxP4(I,Mv) + VtxP4(I,Nv)
        END DO
        END IF
      END DO
C  Check the quality of satisfaction of the pointing constraints, this step
      DO 220 Nv=1,NVERTX
      Mv = VtxPnt(Nv,1)
      IF (Mv.LT.0) GO TO 220
C IMAX points to larger PSUM component (1=x,2=y)
        IMAX = 1
        IF (ABS(VtxP4(1,Nv)).LT.ABS(VtxP4(2,Nv))) IMAX = 2
        NPC = VtxPnt(Nv,2)
        IF (NPC.EQ.3) NPC = 2
        IF (NPC.EQ.4) NPC = 2
        DO Np=1,NPC
          IF (Np.EQ.1) THEN
C  Residual of first pointing constraint, this step
            SUM  = VtxP4(1,Nv)*DELVS(2,Nv) - VtxP4(2,Nv)*DELVS(1,Nv)
            PMAG = SQRT(VtxP4(1,Nv)**2+VtxP4(2,Nv)**2)
            XMAG = SQRT(DELVS(1,Nv)**2+DELVS(2,Nv)**2)
        ELSE
C  Residual of second pointing constraint, this step
          SUM  = DELVS(3,Nv)*VtxP4(IMAX,Nv) - DELVS(IMAX,Nv)*VtxP4(3,Nv)
          PMAG = SQRT(VtxP4(IMAX,Nv)**2 + VtxP4(3,Nv)**2)
          XMAG = SQRT(DELVS(IMAX,Nv)**2+DELVS(3,Nv)**2)
        END IF
        PCON(Nv,Np) = SUM
C Check for very small vertex separation
        IF (XMAG .GT. .001) THEN
C OK.  Get sin(angle)
          SANG(Nv,Np) = PCON(Nv,Np)/(PMAG*XMAG)
        ELSE
          SANG(Nv,Np) = 2.0
        ENDIF
      END DO
C end vertex loop
  220 CONTINUE
C  Check the quality of satisfaction pf the mass constraints, this step
      DO Nm=1,NMASSC
C  Residual of mass constraint.
        SUM = SQRT(McnP4(1,Nm)**2+McnP4(2,Nm)**2+McnP4(3,Nm)**2)
        SUM = (McnP4(4,Nm)+SUM) * (McnP4(4,Nm)-SUM)
        IF (SUM.LE.0.0) THEN
          SUM = 0.0
        ELSE
          SUM = SQRT(SUM)
        END IF
        FMCDIF(Nm) = (SUM + CMASS(Nm)) * (SUM - CMASS(Nm))
        MFTEST(Nm) = 0.5 * FMCDIF(Nm)
C  Convert to Delta m in MeV from Delta m**2 in GeV
        FMCDIF(Nm) = 500.*FMCDIF(Nm)/CMASS(Nm)
      END DO
C  Check the improvements(?) in pointing and mass constraints, this step
C  If step size is too large the constraints will not be satisfied.
C  If there is a problem, the step size will be halved
C Check pointing constraints
      DO 230 Nv=1,NVERTX
      Mv = VtxPnt(Nv,1)
      IF (Mv.LT.0) GO TO 230
        NPC = VtxPnt(Nv,2)
        IF (NPC.EQ.3) NPC = 2
        IF (NPC.EQ.4) NPC = 2
        DO I=1,NPC
          C = PCON(Nv,I) / (PCSAVE(Nv,I) + 0.001)
          IF ((ABS(C).GT.1.0) .AND. (ABS(SANG(Nv,I)).GT..01)) THEN
D           IF (DBGPRT.GT.1) THEN
D           WRITE(PRINT,1051) ITER,ISCUT,Nv,Mv,I,C,SANG(Nv,I)
D1051 FORMAT(' Iter',I2,' Cut',I2,',  Nv,Mv,I',2I3,I4,2F10.5)
D           END IF
            GO TO 250
          END IF
        END DO
C End vertex loop
  230 CONTINUE
C Check mass constraints
      DO Nm=1,NMASSC
        C = ABS(MFTEST(Nm) / (MFSAVE(Nm) + 0.1))
        IF (C.GT.1.0) THEN
D           IF (DBGPRT.GT.1) THEN
D           WRITE(PRINT,1052) ITER,ISCUT,Nm,C,MFTEST(Nm)
D1052 FORMAT(' Iter',I2,' Cut',I2,',  Nm',I2,2F10.5)
D           END IF
          GO TO 250
        END IF
      END DO
      GO TO 300
  250 CONTINUE
C Count a cut step
      ISCUT = ISCUT + 1
      NTSCUT = NTSCUT + 1
C maximum number of cut steps exceeded
      IF (ISCUT .GT. NMAXCT) THEN
        IERR = 9
        IJKERR(2) = 3
        IJKERR(3) = ITER
        WRITE(STRING,2050) IJKERR
        GO TO 990
      END IF
C Prepare to subtract off 1/2 step size taken
      CUTFAC = -0.5*PARFAC
      GO TO 110
C Accept this step
  300 CONTINUE
      CHISQR(0) = CHISQR(ITER)
      STPTK = 1.0
      IF (ITER.GT.0) THEN
      STPTK = 0.0
      DO NT=1,NTRACK
C make the difference between this step and the previous
        DO I=1,5
            PARDJF(I,NT) = PARDIF(I,NT) - PARDJF(I,NT)
            STEP(I) = PARDJF(I,NT)/PARERR(I,NT)
            STPTK = AMAX1(STPTK,ABS(STEP(I)))
        END DO
      END DO
      END IF
D     IF (DBGPRT.GT.1) THEN
D     C = AMIN1(CHISQR(0),9999.0)
D     L = 1
D     IF (PRIMARY) L=0
D     DO I=L,NVERTX
D       WORK(I) = AMIN1(CHIV(I),999.0)
D     END DO
D     IF (NVERTX.EQ.1 .AND. L.EQ.1) THEN
D       WRITE(PRINT,1070) stptk, C
D1070 FORMAT(/,' StpTk ',F8.3,';  ChiSqr =',F7.2)
D     ELSE
D       WRITE(PRINT,1071) stptk, C, (I,WORK(I),I=L,NVERTX)
D1071 FORMAT(/,' StpTk ',F8.3,';  ChiSqr =',F7.2, 4(3X,'ChiV',I2,F7.2) )
D     END IF
D     IF (ITER.EQ.0) GO TO 500
D       DO NT=1,NTRACK
D       WORK(Nt) = AMIN1(CHIT(Nt),999.0)
D       DO I=1,5
D         STEP(I) = PARDJF(I,NT)/PARERR(I,NT)
D       END DO
D     WRITE(PRINT,1075) nt, WORK(Nt), STEP
D1075 FORMAT(i5,F10.4,2x,5F10.6)
D       END DO
D     END IF
C Check for step acceptability
  500 CONTINUE
CP>>
c
c     Do the calculation myself here
c
      DO I=1,3                          ! vertex vector
        DV(I) = XYZVRT(I,1) - XYZPV0(I)
      ENDDO 
      MAGN = 0.0
      LXY(ITER) = 0.0
      DO I=1,2                          ! length of momentum, Xvtx dot Ptvtx
        MAGN      = MAGN      + VtxP4(I,1)*VtxP4(I,1)
        LXY(ITER) = LXY(ITER) + VtxP4(I,1)*DV(I)
      ENDDO
      LXY (ITER) = LXY(ITER)/SQRT(MAGN) ! project to unit length
      CHI2(ITER) = CHISQR(0)            ! fit quality in 3d
      CRP2(ITER) = 0.0          ! fit quality quasi in 2d
      DO I=1,NTRACK
        DO J=1,3
          DO K=1,3
            CRP2(ITER) = CRP2(ITER) +
     &        PARDIF(IXCRP2(J),I)*G(IXCRP2(K),IXCRP2(J),I)*
     &        PARDIF(IXCRP2(K),I);
          ENDDO
        ENDDO
      ENDDO
c
c     Implement early termination 
c
      IF (ITER.EQ.0) THEN
        IF ((C20MAX.GT.0.0 .AND. CRP2(0).GT.C20MAX) .OR.
     &      (                     LXY(0).LT.LXY0MN)     ) THEN
          IERR = 9
          GOTO 990
        ENDIF
      ENDIF
      IF (ITER.EQ.1) THEN
        IF ((C21MAX.GT.0.0 .AND. CRP2(1).GT.C21MAX) .OR.
     &      (                     LXY(1).LT.LXY1MN)     ) THEN
          IERR = 9
          GOTO 990
        ENDIF
      ENDIF
CP<<
      IF (ITER.EQ.0) GO TO 550

      IF (STPTK.LT.0.1) THEN
C accept the convergence
        GO TO 900
      ELSE IF (ITER.GE.MAXITR) THEN
C "convergence", sort of. Accept the stinker?
        IF (STPTK.LT.0.2) THEN
          GO TO 900
        ELSE
          STPTK = AMIN1(STPTK,9999.0)
          IERR = 9
          IJKERR(2) = 4
          WRITE(STRING,2090) ITER,StpTk,CHISQR(0)
 2090     FORMAT('; Convergence failure, Itr =',I2,'; StpTk ',F7.2
     &,         '   ChiSqr(i) ', 1P,E9.2)
          GO TO 990
        END IF
      ELSE
        IF (ITER.GT.2) THEN
C after the 2nd step, if the "ChiSquare"
          TEST = FLOAT(NDOF) * chisqmx
C  is >15 "standard deviations", and if
          IF (CHISQR(0).GT.TEST) THEN
C   the last step was did not do much good
          IF (STPTK.GT.0.5) THEN
            STPTK = AMIN1(STPTK,9999.0)
C         give up in disgust
            IERR = 9
            IJKERR(2) = 5
            WRITE(STRING,2095) ITER,StpTk,CHISQR(0)
 2095     FORMAT('; Bad convergence, Itr =',I2,'; StpTk ',F7.2
     &,          '   ChiSqr(i) ', 1P, 6E9.2)
            GO TO 990
          END IF
          END IF
        END IF
      END IF
  550 CONTINUE
C return for the next iteration step
      DO NT=1,NTRACK
C save the last step on this track
        DO I=1,5
          PARDJF(I,NT) = PARDIF(I,NT)
        END DO
      END DO
      ITER = ITER + 1
      GO TO 100
C  Successful Return
C------------------------------------------------------------------------------
  900 CONTINUE
C final track parameters
      DO Nt=1,NTRACK
      DO I=1,5
        PAR(I,NT) = PAR0(I,Nt) + PARDIF(I,Nt)
      END DO
D        IF (DBGPRT.GT.0) THEN
D        write(print,5170)NT,(par(i,nt),i=1,5)
D5170   format(1x,'Final track ',i3,' params ',5(1x,e10.4))
D       ENDIF
      END DO
      IF (PRINT.GT.0) THEN
        CALL CCTVMPR(PRINT,DBGPRT,PARERR)
      END IF
cccCP>>
ccc
ccc     Book keeping for optimization
ccc
cccc      WRITE(*,*)
cc      NTOT = NTOT + 1
cc      NFIT(ITER) = NFIT(ITER) + 1
ccc
ccc     check general quality and reject accordingly
cc      LOST = .FALSE.
cc      IF ((C20MAX.GT.0.0 .AND. CRP2(0).GT.C20MAX) .OR.
cc     &    (C21MAX.GT.0.0 .AND. CRP2(1).GT.C21MAX) .OR.
cc     &     LXY(0).LT.LXY0MN .OR.  LXY(1).LT.LXY1MN     )  THEN
cc        NOPT(0) = NOPT(0) + 1
cc        IF (CRP2(ITER).LT.30.0 .AND. LXY(ITER).GT.0.015) THEN
cc          PTVTX = SQRT(VtxP4(1,1)*VtxP4(1,1)+VtxP4(2,1)*VtxP4(2,1))
cc          PTTK1 = SQRT(TRKP4(1,1)*TRKP4(1,1)+TRKP4(1,2)*TRKP4(1,2))
cc          PTTK2 = SQRT(TRKP4(2,1)*TRKP4(2,1)+TRKP4(2,2)*TRKP4(2,2))
cc          PTTK3 = SQRT(TRKP4(3,1)*TRKP4(3,1)+TRKP4(3,2)*TRKP4(3,2))
cc          WRITE(*,'(A,12F11.5)') ' ERROR -- !! Lost Vertex !! ',
cc     &      CHI2(0),CHI2(1),CRP2(0),CRP2(1),LXY(0),
cc     &      CHI2(ITER),CRP2(ITER),LXY(ITER),PTVTX,PTTK1,PTTK2,PTTK3
cc          NERR = NERR   + 1
cc          LOST = .TRUE.
cc        ENDIF
cc      ELSE
cc        NOPT(ITER) = NOPT(ITER) + 1
cc      ENDIF
ccc
ccc     start with too large chi2
cc      IF     (CRP2(0).GT.C20MAX .OR. CRP2(1).GT.C21MAX) THEN
cc        NOPT(0) = NOPT(0) + 1
cc        IF ((CHI2(ITER).LT.30)) THEN
cc          NERRC2 = NERRC2 + 1
cc        ENDIF
ccc
ccc     next too small Lxy
cc      ELSEIF (LXY(0).LT.LXY0MN .OR. LXY(1).LT.LXY1MN) THEN
cc        NOPT(0) = NOPT(0) + 1
cc        IF (LXY(ITER).GT.0.015) THEN
cc          NERRLN = NERRLN + 1
cc        ENDIF
cc      ENDIF
ccc
ccc     Printout the stuff
ccc
ccc      WRITE(*,'(A,3F10.5)') 'VtxPv ',XYZPV0(1),XYZPV0(2),XYZPV0(3)
ccc      WRITE(*,'(A,3F10.5)') 'Vtx1  ',XYZVRT(1,1),XYZVRT(2,1),XYZVRT(3,1)
ccc      WRITE(*,'(A,3F10.5)') 'P     ',VTXP4 (1,1),VTXP4 (2,1),VTXP4 (3,1)
ccc
ccc     Test 'no fit left behind' policy
ccc
cc      IF (LOST) THEN
cc        WRITE(*,'(3(A,F10.5))')
cc     &    ' iChi2: ',CHI2(0),' iCrp2: ',CRP2(0),'  iLxy: ',LXY(0)
cc        DO I=1,ITER
cc          WRITE(*,'(3(A,F10.5))')
cc     &      ' dChi2: ',CHI2(I)-CHI2(0),' dCrp2: ',CRP2(I)-CRP2(0),
cc     &      '  dLxy: ',LXY(I)-LXY(0)
cc        ENDDO
cc        WRITE(*,'(3(A,F10.5))')
cc     &    ' fChi2: ',CHI2(ITER),' fCrp2: ',CRP2(ITER),
cc     &    '  fLxy: ',LXY (ITER)
cc      ENDIF
cc        
cc      IF (MOD(NTOT,10000) .EQ. 0) THEN
cc        WRITE(*,'(12I9)') NTOT,(NFIT(I),I=0,MAXITR)
cc        WRITE(*,'(12I9)') NERR,(NOPT(I),I=0,MAXITR)
cc        WRITE(*,'( 3I9)') NERRC2,NERRLX,NERRLN
cc      ENDIF
ccCP<<
C  Check for peculiar case of non-positive diagonal matrix elements, which
C  should never happen, but if it happens, better to report it than to
C  let it slip by unnoticed.  All errors of this sort should be reported
C  immediately to the guilty parties.
      DO I = 1,TOFF(NTRACK)+3
        IF (VMAT(I,I).LE.0) THEN
          STRING = 'SERIOUS problem with ill-formed covariance matrix'
          IERR = 9
          IJKERR(1) = 9
          IJKERR(2) = 9
          CALL CERROR('CTVMFT',IERR,STRING)
          IF (PRINT.GT.0) WRITE (PRINT,*) MRUN,MNEV,I,VMAT(I,I)
          GOTO 999
        END IF
      END DO
      RETURN
C  Error Return
C-------------------------------------------------------------------------------
  990 CONTINUE
      IF (IERR.NE.0)  THEN
      IJKERR(1) = IERR
      NAME = 'CTVMFT'
ccc      write(*,9001) name,ierr,string
ccc 9001 format(1x,A,": Ierr = ",i3," => ",A)
D     IF (PRINT.GT.0 .OR. IERR.GT.50) THEN
D       CALL CERROR (NAME,IERR,STRING)
D     END IF
      END IF
  999 CONTINUE
      RETURN
      END
      SUBROUTINE CERROR (NAME,IERR,STRING)
      INTEGER IERR
      CHARACTER*6  NAME
      CHARACTER*80 STRING
      INTEGER COUNT
      SAVE    COUNT
      DATA    COUNT /0/
ccc      WRITE(*,1000) NAME,IERR,STRING
ccc 1000 FORMAT(1X,'subroutine ',A6,', Err',I3,' comment; ',A80)
      COUNT = COUNT+1
      RETURN
      END
C===============================================================================
      SUBROUTINE CCTVM00 (IERR, WXYZPV)
C===============================================================================
C  Checks the current input specifications at the beginning of this call to
C  CTVMFT: (Vertex, conversion, Pointing, Mass Constraints) compound fit.
C  Any error detected in CTVM00 is a structural failure on the input conditions
C  and is (by default) treated as a fatal error. We attempt to crash the job
C  (after printing a specific error condition message). The studious who are
C  running interactively can discover how to avoid this crash; cf. the variables
C  "EXCUSE" and "ME", set in a local common block in CTVMFT.
C  If an error is detected, IERR = IJKERR(1) = 1
C                                  IJKERR(2) = the precise detected trouble.
C  WXYZPV is returned, it is the primary vertex double precision weight matrix
C  VtxVtx, the vertex "geneology", association structure is made and returned.
C--Input data in the include file CTVMFT (COMMON /CTVMFC/)
C  NVERTX  Number of vertices in this fit (MUST be at least 1)
C  NTRACK  Number of tracks in this fit.  NTRACK MUST be at least 2.
C  TrkVtx(MAXTRK,MAXVTX)    (Logical)  See the discussion in CTVMFT.
C  VtxPnt(MAXVTX,2)         (Integer)  See the discussion in CTVMFT.
C  TrkMcn(MAXTRK,NMASSC)    (Logical)  See the discussion in CTVMFT.
C  CMASS(NMASSC)            (Real)     See the discussion in CTVMFT.
C  TKBANK(MAXTRK)           (Character*4)
C         Bank name (exempli gratia, 'QTRK','TRKS','SVXS','CTCS') from which
C         to get each track`s helix parameter and covariance data (GETTRK).
C  LIST  (MAXTRK)           (Integer)
C         Bank numbers of the tracks to be used.
C         Note that if LIST(i) is a negative integer it is assumed that the
C         information for this track is present unchanged from a previous call.
C  Edit log:
C  ---- ---
C  11/28/94  WJA  Require exactly 2 tracks for conversion vertex; always
C                 require 1 or more tracks per vertex; If there is a "single
C                 track vertex", this vertex must be contain at least one
C                 multi-track exclusive vertex in its descendant chain.
C  30 Jan 2003 F. Bedeschi: add checks specific to beamline constraint option
C  04 FEB 2003 J. Guimaraes: Fixed bug in fiddle common block (missing chisqmx)
C==============================================================================
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
#include "MitCommon/Ctvmft/interface/cctvmfi.h"
      INTEGER IERR
      REAL*8 WXYZPV(3,3)
      INTEGER Nv,Mv,Kv, Nt,Kt, Np,Nm
      INTEGER I,J, II,JJ
      INTEGER NtrVtx(MaxVtx)
      REAL    SUM, WORK(MAXDIM)
      REAL    SCR
      CHARACTER*4  NAME
      CHARACTER*80 STRING
      REAL ZERO,MINUS,SHOUT
      SAVE ZERO,MINUS
      DATA ZERO /0.0/, MINUS /-1.0/
C  Local variables pertaining to zero-track vertices
      INTEGER NMULTI
C--------------- Check Input Conditions for illegal conditions -----------------
C

      IF ((NTRACK.LT.2 .AND. VTXPNT(1,1) .NE. -100) .OR.
     >    (NTRACK.LT.1 .AND. VTXPNT(1,1) .EQ. -100) .OR.
     >    NTRACK.GT.MAXTRK) THEN
        IJKERR(2) = 1
        WRITE(STRING,1001) NTRACK
 1001 FORMAT('Inadequate or improper number of tracks;',I2)
        GO TO 900
      END IF
      IF (NVERTX.LT.1 .OR. NVERTX.GT.MAXVTX) THEN
        IJKERR(2) = 2
        WRITE(STRING,1002) NVERTX
 1002 FORMAT('Inadequate or improper number of vertices;',I2)
        GO TO 900
      END IF
      IF (NMASSC.LT.0 .OR. NMASSC.GT.MAXMCN) THEN
        IJKERR(2) = 3
        WRITE(STRING,1003) NMASSC
 1003 FORMAT('Improper number of mass constraints;',I2)
        GO TO 900
      END IF
C count the number of measured tracks
      DO Nv=1,MaxVtx
C present at each vertex
        NtrVtx(Nv) = 0
        Do Nt=1,NTRACK
          IF (TrkVtx(Nt,Nv)) NtrVtx(Nv) = NtrVtx(Nv)+1
        END DO
      END DO
C make the vertex geneology table
      DO Kv=NVERTX,1,-1
        DO Nv =NVERTX,1,-1
          VtxVtx(Kv,Nv) =.FALSE.
        END DO
        VtxVtx(Kv,Kv) = .TRUE.
C Kv is a descendant of its parent,
    1   IF (Nv.GT.0) THEN
C  its parent`s parent, ...
          VtxVtx(Kv,Nv) = .TRUE.
          Nv = VtxPnt(Nv,1)
          GO TO 1
        END IF
      END DO
C  Check the logic of vertex specifications
      DO 10 Kv=NVERTX,1,-1
        IJKERR(3) = Kv
C exactly 2 tracks for conversion
        IF (Cvtx(Kv).NE.0) THEN
          IF (NtrVtx(Kv).NE.2) THEN
            IJKERR(2) = 12
            WRITE (STRING,'(A,I2,A,I2,A)')
     &        'Conversion vertex ',Kv,' with ',NtrVtx(Kv),' tracks'
            GO TO 900
          END IF
        END IF
C   If not a beamline constrained fit  -- Joao 02/09/03
        IF (VtxPnt(1,1) .NE. -100) THEN
C   a one-track vertex must have a
           IF (NtrVtx(Kv).EQ.1) THEN
C   descendant with at least 2 tracks
              Mv = Kv
              DO Nv=Kv+1,NVERTX
                IF (VtxPnt(Nv,1).EQ.Mv) THEN
C                  override the vertexpointing for 1track to only 1track
C                  IF (NtrVtx(Nv).GT.1) GO TO 5
                  IF (NtrVtx(Nv).GT.0) GO TO 5
C Craig Blocker 1/9/03 - Changed the following line from
C Mv = VtxPnt(Nv,1)
C to get proper test done here.
                    Mv = Nv
                 END IF
              END DO
              IJKERR(2) = 13
              WRITE(STRING,'(A,I2,A,I2,A)')
     &             'Vertex ',Kv,' with ',NtrVtx(Kv),
     &             ' tracks is not legally pointed at'
              GO TO 900
C   a zero-track vertex must have 2 or more
           ELSE IF (NtrVtx(Kv).EQ.0) THEN
C   descendants with at least 2 tracks.
              NMULTI = 0
              Mv = Kv
              DO Nv=Kv+1,NVERTX
                 IF (VtxPnt(Nv,1).EQ.Mv) THEN
                    IF (NtrVtx(Nv).GT.1) NMULTI = NMULTI + 1
                 END IF
              END DO
C   If there are only 2 vertices pointed to this one
C   AND they both have >1 tracks then this is a valid
C   zero track vertex.
              IF (NMULTI.GE.2) GO TO 5
              IJKERR(2) = 13
              WRITE(STRING,'(A,I2,A,I2,A)')
     &             'Vertex ',Kv,' with ',NtrVtx(Kv),
     &             ' tracks is not legally pointed at'

              GO TO 900
           END IF
        END IF
C Check that Kv`s pointing is valid
    5   CONTINUE
        Nv = VtxPnt(Kv,1)
        IF (Nv.LE.0) THEN
          IF (Nv.LT.0) VtxPnt(Kv,2) = -1
        ELSE
C Daughter vertex number must exceed
          IF (Nv.GE.Kv) THEN
C   mother vertex number
            IJKERR(2) = 14
            WRITE (STRING,'(A,I2,A,I2)')
     &        'Vertex ',Kv,' has illegal target vertex number ',Nv
            GO TO 900
          END IF
C conversion vertex cannot be a target
          IF (Cvtx(Nv).GT.0) THEN
            IJKERR(2) = 15
            WRITE (STRING,'(A,I2,A)')
     &        'Vertex ',Kv,' has illegal target vertex type'
            GO TO 900
          END IF
C  Single-track pointing if and only if to single-track vertex
          IF (VtxPnt(Kv,2).EQ.3.NEQV.NtrVtx(Nv).EQ.1) THEN
            IJKERR(2) = 16
            WRITE (STRING,'(I2,A,I2,A)')
     &        Kv,' ==> ',Nv,' is an illegal 1-track vertex pointing'
            GO TO 900
          END IF
C  zero-track pointing if and only if to zero-track vertex
          IF (VtxPnt(Kv,2).EQ.4.NEQV.NtrVtx(Nv).EQ.0) THEN
            IJKERR(2) = 17
            WRITE (STRING,'(I2,A,I2,A)')
     &        Kv,' ==> ',Nv,' is an illegal 0-track vertex pointing'
            GO TO 900
          END IF
        END IF
   10 CONTINUE
      IJKERR(3) = 0
C  If VtxPnt requires a primary vertex, checks that the input primary vertex
C  covariance matrix is "reasonable" (id est, non-singular).
      NP = 0
      DO NV=1,NVERTX
        IF (VtxPnt(Nv,1).EQ.0) NP = NP+1
      END DO
      IF (NP.GT.0) THEN
        DO I=1,3
C primary vertex input position
          XYZPV(I)  = XYZPV0(I)
          DO J=1,3
C primary vertex error matrix
            WXYZPV(J,I)= EXYZPV(J,I)
          END DO
        END DO
C make the primary vertex weight matrix
        CALL DINV(3,WXYZPV,3,WORK,I)
        IF (I.NE.0) THEN
          IJKERR(2) = 19
          WRITE(STRING,'(A)')
     &       'The input primary vertex covariance matrix is singular.'
          GO TO 900
        END IF
      END IF
C  Check the specification of tracks
C every track must be at a vertex
      DO 20 Nt=1,NTRACK
        IJKERR(3) = Nt
        IF (LIST(Nt).LE.0) THEN
          IJKERR(2) = 20
C some track failed TkSlct
          IF (LIST(Nt).LT.0) THEN
            IJKERR(1) = 2
            IERR = 2
            RETURN
          END IF
          WRITE(STRING,'(A,I2,A)') 'Track ',Nt,' is not specified'
          GO TO 900
        END IF
        IF (TKBANK(NT).EQ.'QTRK') GO TO 11
        IF (TKBANK(NT).EQ.'TRKS') GO TO 11
        IF (TKBANK(NT).EQ.'TRKV') GO TO 11
        IF (TKBANK(NT).EQ.'SVXS') GO TO 11
        IF (TKBANK(NT).EQ.'CTCS') GO TO 11
          IJKERR(2) = 20
          WRITE(STRING,1005) NT,LIST(NT),TKBANK(NT)
          GO TO 900
 1005 FORMAT('Track ',I2,',',I4,2X,A4,' is improperly specified.')
   11   CONTINUE
        II = 0
        DO NV=1,NVERTX
          IF (TrkVtx(NT,NV)) II=II+1
        END DO
C every track belongs to a vertex
        IF (II.LT.1) THEN
          IJKERR(2) = 21
          WRITE (STRING,'(A,I2,A)')
     &    'track ',Nt,' is not at a vertex.'
          GO TO 900
C and to only one vertex
        ELSE IF (II.GT.1) THEN
          IJKERR(2) = 22
          WRITE(STRING,'(A,I2,A)')
     &    'Track ',Nt,' appears at two vertices'
          GO TO 900
        END IF
C each track may appear only once
        DO Kt=Nt+1,NTRACK
          IF (IABS(LIST(Nt)).EQ.IABS(LIST(Kt))) THEN
          IF (TkBank(Nt).EQ.TkBank(Kt)) THEN
            IJKERR(2) = 23
            WRITE(STRING,'(A,I2,A,I2,A)')
     &       'Track ',Nt,' and Track ',Kt,' are identical'
            GO TO 900
          END IF
          END IF
        END DO
   20 CONTINUE
      IJKERR(3) = 0
C Check the Mass Constraint/Track/Vertex specifications
      DO NM=1,NMASSC
        II = 0
        SUM = 0.0
        DO NT=1,NTRACK
          IF (TrkMcn(NT,NM)) THEN
            II=II+1
            SUM = SUM + TMASS(NT)
          END IF
        END DO
C check for enough tracks
        IF (II.LT.2) THEN
          IJKERR(2) = 31
          WRITE(STRING,1021) NM
 1021 FORMAT('Mass Constraint',I2,' has too few tracks.')
          GO TO 900
        END IF
C check for possible mass constraint
        IF (CMASS(NM).LE.SUM) THEN
          IJKERR(2) = 32
          WRITE(STRING,1022) NM,CMASS(NM),SUM
 1022 FORMAT('Mass Constraint',I2,', mass',F6.3,' has track mass0',F6.3)
          GO TO 900
        END IF
      END DO
C
C  Check for beamline constraint specifications
C
      IF (VTXPNT(1,1) .EQ. -100) THEN
         SCR = 0
         DO I = 1,3
            DO J = 1,3
               WXYZPV(I,J) = EXYZPV(I,J)
               IF (I .NE. J)SCR = SCR + ABS(EXYZPV(I,J))
            ENDDO
         ENDDO
         IF (EXYZPV(1,1) .LE. 0.0 .OR.
     >       EXYZPV(2,2) .LE. 0.0 .OR.
     >       EXYZPV(3,3) .LE. 0.0) THEN
            IJKERR(2) = 33
            WRITE(STRING, 1023)EXYZPV(1,1), EXYZPV(2,2), EXYZPV(3,3)
 1023 FORMAT('Beamline constraint.COVARIANCE NOT SET PROPERLY. C(I,I)='
     >        ,3(F6.3,2X))
            GO TO 900
         ENDIF
         IF (SCR .NE. 0.0) THEN
            IJKERR(2) = 34
            WRITE(STRING, 1024)
 1024  FORMAT('Beamline constraint. COVARIANCE has non diagonal terms.')
            GO TO 900
         ENDIF
C     Make primary vertex weight matrix
         CALL DINV(3,WXYZPV,3,WORK,I)
         IF (I .NE. 0) THEN
            IJKERR(2) = 35
            WRITE(STRING, 1025)
 1025   FORMAT('Beamline constraint. PV covariance cannot be inverted.')
            GO TO 900
         ENDIF
         IF (NVERTX .NE. 1)THEN
            IJKERR(2) = 36
            WRITE(STRING, 1026)NVERTX
 1026       FORMAT('Beamline constraint. NVERTX not 0: = ',F6.3)
            GO TO 900
         ENDIF
      ENDIF
 100  CONTINUE
      IERR = 0
      IJKERR(1) = 0
      RETURN
C     Executive action; Terminate with extreme prejudice.
 900  CONTINUE
      IERR      = 1
      IJKERR(1) = 1
ccc      PRINT 1999, IJKERR
ccc 1999 FORMAT('Improper input specification, IJKERR=',3I3)
      CALL CERROR ('CTVM00',IERR,STRING)
C Escape hatch, for debuggery
      IF (EXCUSE.EQ.0) THEN
        SHOUT = SQRT(MINUS) / ZERO
        STOP
      ELSE
        RETURN
      END IF
      END
C==============================================================================
      SUBROUTINE CCTVMFA (Nv, PRINT,IERR)
C==============================================================================
C
C      Calculates the initial approximation to the vertex position.
C      This subroutine tries to intersect circles that are (R,Phi) projections
C      of the helical trajectories of the two tracks specified in the calling
C      sequence.
C      If the circles intersect, the intersection with smaller z difference
C      between the 2 helices extrapolated to the intersection is chosen.
C      If the circles do not intersect, the vertex approximation is taken as
C      the point of closest appraoch of the two circles.  Both cases of non-
C      intersection ("interiour", "exteriour") are treated.
C      Note that the tracks are required to travel "forwards" to the vertex,
C      and this requirement may override the "smallest delta Z" criterion.
C      Tests are made on the vertex radius and Z difference between the two
C      tracks, and error codes are returned for unacceptable vertices.
C      Note that, where possible, the calculation is completed, even if an
C      error condition is flagged.
C      Test parameter settings used to accept/reject vertex approximations;
C                  DRMAX,DZMAX,RVMAX, TRNMAX,DSMIN
C      appear as data statements in the include file CTVMFT.INC.
C      IERR = 2  = IJKERR(1)  marks failure in vertex first approximation
C  IJKERR(2)= 1  Concentric circles, in the first approximation
C             2  (conversion) widely separated exteriour circles at the midpoint
C             3  (conversion) widely separated interiour circles at the midpoint
C             4  widely separated exteriour circles at the approximate vertex
C             5  widely separated interiour circles at the approximate vertex
C             6  Rv       is too large at the chosen intersection point.
C             7  delta(Z) is too large at the chosen intersection point.
C             8  a track turning angle to the chosen vertex is too large
C            19  no acceptable solution with an adequately positive arc length
C            21  zero-track vertexing: either/both vertex momenta are too small (<0.01 MeV)
C            22  zero-track vertexing: Two lines are parallel/antiparallel
C  Edit log:
C  ---- ---
C  10/xx/94  WJA  Allow NTRK = 1 or 2, for possible 1-track vertex; add
C                 line-intersects-circle approximation for 1-track vertex;
C                 add split-the-difference approximation for conversions;
C                 detect negative JJ() values to trigger these special-case
C                 first approximations.
C  05/xx/02  MSM  Allow NTRK = 0, for possible zero-track vertex; add
C                 line-intersects-line approximation for zero-track vertex;
C                 Detect JJ(1)=JJ(2)=0 to trigger this.
C===============================================================================
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmco.h"
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
#include "MitCommon/Ctvmft/interface/cctvmuv.h"
      INTEGER Nv, PRINT, IERR
      INTEGER I,J,K,L,M,N, JJ(2),II(2), NSOL, KPRINT
      INTEGER Nt, Mv, NTRK
      REAL    HELIX(5,2), P(5),Q(5,5)
      REAL    A,B, RA(2), RVI(2), TANG(2,2)
      REAL    Xsvt,Ysvt, AA,BB,CC,DD,RR, Pt
      REAL    TMP,PZV,FITA,FITB
      COMPLEX XYP0,PXYV,XYVR,PXYVR,XYTMP
C  Local variables for accumulating first-approximation momentum sums
      REAL    SPh0,CPh0,SPhs,CPhs,SPhi,CPhi
      LOGICAL CONV
      REAL*8  XC(2),YC(2),RC(2)
      REAL*8  AB(2), XX,YY,Y0,YY2,COST,SINT,DX,DY,D,U,V
      INTEGER BEAM
C Map the parameter order to CTVMFT form:
      INTEGER MAP(5)
      SAVE    BEAM, MAP
      DATA    BEAM /0/
C (Ctg,Crv,Z0,D0,Phi)->(Crv,Phi,Ctg,D0,Z0)
      DATA    MAP / 3,1,5,4,2 /
C Variables for the first appx for zero track vertexing case
C (line intersects line)
C ZV="Zero-Vertex" :
C First the two un-rotated vertex positions and momenta (2D)
      REAL ZVX1,ZVY1,ZVX2,ZVY2,ZVPX1,ZVPY1,ZVPX2,ZVPY2
C Next the two rotated vertex postions and momenta (2D)
      REAL ZVRX1,ZVRY1,ZVRX2,ZVRY2,ZVRPX1,ZVRPY1,ZVRPX2,ZVRPY2
C Next some double versions for a magnitude check (see code)
      REAL*8 ZVDPX1,ZVDPY1,ZVDPX2,ZVDPY2,ZVDPS1,ZVDPS2
C Next some doubles for calculating sin and cos from dot and cross prod
      REAL*8 ZVDCS1,ZVDCS2,ZVDSN1,ZVDSN2,ZVDRSN
C Next some real versions of the dot and cross products
      REAL ZVCS1,ZVCS2,ZVSN1,ZVSN2
C Next the anticlockwise angles from right x axis to the vertex momentum
      REAL ZVANG1,ZVANG2
C Next the rotation angle, and its sine and cosine
      REAL ZVROTANG,ZVSIN,ZVCOS
C Next the two roated gradients (M) and intercepts (C)
      REAL ZVRM1,ZVRM2,ZVRC1,ZVRC2
C Next some swim variables for the Z calculation
      REAL ZVZSWIM1,ZVZSWIM2,ZVZ1,ZVZ2,ZVPZ1,ZVPZ2
C Finally the rotated point of intersection..and a temp.
      REAL ZVRINTX,ZVRINTY,ZVMDIFF
C local variables to determine the average z of the input track set
      REAL SUMZOS, SUMOOS, ZAVG
C Start of code
C------------------------------------------------------------------------------
      KPRINT = 0
      IF (PRINT.GT.0 .AND. DBGPRT.GT.10) KPRINT=PRINT
      IERR      = 0
      IJKERR(1) = 0
C Special case beamline constraint -Joao 02/11/03
      IF (VTXPNT(1,1) .EQ.-100) THEN
        DO I=1,3
           XYZVRT(I,Nv) = XYZPV0(I)
        END DO
        GOTO 50
      ENDIF
C Find first two tracks (if they exist)
      JJ(1) = 0
      JJ(2) = 0
      DO Nt=1,NTRACK
        IF (TrkVtx(Nt,Nv)) THEN
          IF (JJ(1).EQ.0) THEN
            JJ(1) = Nt
          ELSE
            JJ(2) = Nt
            GO TO 1
          END IF
        END IF
      END DO
  1   CONTINUE
      CONV = .FALSE.
      NTRK = 2
C special case; conversion vertex
      IF (Cvtx(Nv).GT.0) CONV=.TRUE.
C special case; "single track vertex"
      IF (JJ(2).EQ.0 .AND. JJ(1).NE.0) THEN
        DO Mv=Nv+1,NVERTX
          IF (VtxPnt(Mv,1).EQ.Nv) THEN
C For 1-track vertex, JJ(2) is a vertex
            NTRK = 1
C   that points to this vertex
            JJ(2) = Mv
          END IF
        END DO
      END IF
C special case: "zero track vertex"
      IF (JJ(2).EQ.0 .AND. JJ(1).EQ.0) THEN
        NTRK = 0
        DO Mv=Nv+1,NVERTX
          IF (VtxPnt(Mv,1).EQ.Nv) THEN
C For 0-track vertex, JJ(2) and JJ(1) are the vertices
C that point to this vertex
            IF (JJ(1).EQ.0) THEN
              JJ(1) = Mv
            ELSE
              JJ(2) = Mv
            END IF
          END IF
        END DO
      END IF
C Continue if this is a zero track vertex.
C Otherwise leapfrog to the right place.
      IF (NTRK.EQ.0) THEN
C Initialise the position and momenta of 2
C Vertices. This gives the 2 lines to intersect
        ZVX1 = XYZVRT(1,JJ(1))
        ZVY1 = XYZVRT(2,JJ(1))
        ZVX2 = XYZVRT(1,JJ(2))
        ZVY2 = XYZVRT(2,JJ(2))
        ZVPX1 = VtxP4(1,JJ(1))
        ZVPY1 = VtxP4(2,JJ(1))
        ZVPX2 = VtxP4(1,JJ(2))
        ZVPY2 = VtxP4(2,JJ(2))
C And Z:
        ZVZ1 = XYZVRT(3,JJ(1))
        ZVZ2 = XYZVRT(3,JJ(2))
        ZVPZ1 = VtxP4(3,JJ(1))
        ZVPZ2 = VtxP4(3,JJ(2))
C Initialise rotation angle to zero
        ZVROTANG = 0
C Calculate mag**2 of each momenta
        ZVDPX1=ZVPX1
        ZVDPY2=ZVPY1
        ZVDPX2=ZVPX2
        ZVDPY2=ZVPY2
        ZVDPS1=ZVDPX1**2 + ZVDPY1**2
        ZVDPS2=ZVDPX2**2 + ZVDPY2**2
C Require each mag**2>1.0E-10 which means
C the mag>1.0E-5 which is 0.01 MeV which is safe.
        IF(ZVDPS1.LE.1.0E-10 .OR. ZVDPS2.LE.1.0E-10) THEN
C One or both vertex momenta are too small!
          IJKERR(2) = 21
          GO TO 1000
        ENDIF
C Require abs(sin(relative angle))>0.00001 (ie 1/100 millirad)
C To catch parallel momenta
        ZVDRSN = (ZVDPX1*ZVDPY2 - ZVDPX2*ZVDPY1)
     &           /(SQRT(ZVDPS1)*SQRT(ZVDPS2))
        IF(ABS(ZVDRSN).LE.1.0E-5) THEN
C Vertex Momenta are too parallel
          IJKERR(2) = 22
          GO TO 1000
        ENDIF
C Calculate Angle between each momenta and x axis:
        ZVDCS1 = ZVDPX1/(SQRT(ZVDPS1))
        ZVDCS2 = ZVDPX2/(SQRT(ZVDPS2))
        ZVDSN1 = ZVDPY1/(SQRT(ZVDPS1))
        ZVDSN2 = ZVDPY2/(SQRT(ZVDPS2))
C Explicit cast to real
        ZVCS1 = ZVDCS1
        ZVCS2 = ZVDCS2
        ZVSN1 = ZVDSN1
        ZVSN2 = ZVDSN2
C Pick right trig solution:
C First Vertex Momentum:
        IF(ZVCS1.LT.0.0) THEN
          IF(ZVSN1.LT.0.0) THEN
            ZVANG1 = PI - ACOS(ZVCS1)
          ELSE
            ZVANG1 = ACOS(ZVCS1)
          ENDIF
        ELSE
          IF(ZVSN1.LT.0.0) THEN
            ZVANG1 = TWOPI - ACOS(ZVCS1)
          ELSE
            ZVANG1 = PI + ACOS(ZVCS1)
          ENDIF
        ENDIF
C And second:
        IF(ZVCS2.LT.0.0) THEN
          IF(ZVSN2.LT.0.0) THEN
            ZVANG2 = PI - ACOS(ZVCS2)
          ELSE
            ZVANG2 = ACOS(ZVCS2)
          ENDIF
        ELSE
          IF(ZVSN2.LT.0.0) THEN
            ZVANG2 = TWOPI - ACOS(ZVCS2)
          ELSE
            ZVANG2 = PI + ACOS(ZVCS2)
          ENDIF
        ENDIF
C Set rotation angle:
        ZVROTANG = ZVANG1 + 0.5*(ZVANG2 - ZVANG1)
C Do the rotation:
        ZVSIN = SIN(ZVROTANG)
        ZVCOS = COS(ZVROTANG)
        ZVRX1 = ZVX1*ZVCOS - ZVY1*ZVSIN
        ZVRY1 = ZVX1*ZVSIN + ZVY1*ZVCOS
        ZVRX2 = ZVX2*ZVCOS - ZVY2*ZVSIN
        ZVRY2 = ZVX2*ZVSIN + ZVY2*ZVCOS
        ZVRPX1 = ZVPX1*ZVCOS - ZVPY1*ZVSIN
        ZVRPY1 = ZVPX1*ZVSIN + ZVPY1*ZVCOS
        ZVRPX2 = ZVPX2*ZVCOS - ZVPY2*ZVSIN
        ZVRPY2 = ZVPX2*ZVSIN + ZVPY2*ZVCOS
C Now solve intersection in "safe" coordinates
C (ie ZVRPX1 and ZVRPX2 should both be non-zero)
C Calculate gradients and intercepts:
        ZVRM1 = ZVRPY1/ZVRPX1
        ZVRM2 = ZVRPY2/ZVRPX2
        ZVRC1 = ZVRY1 - ZVRM1*ZVRX1
        ZVRC2 = ZVRY2 - ZVRM2*ZVRX2
C Calculate the rotated intersection:
        ZVRINTX = (ZVRC2-ZVRC1)/(ZVRM1-ZVRM2)
        ZVRINTY = (ZVRM1*ZVRINTX + ZVRC1)
C Now rotate back:
        ZVROTANG=-1.0*ZVROTANG
        ZVSIN = SIN(ZVROTANG)
        ZVCOS = COS(ZVROTANG)
        XSVI(1) = ZVRINTX*ZVCOS - ZVRINTY*ZVSIN
        YSVI(1) = ZVRINTX*ZVSIN + ZVRINTY*ZVCOS
C Now calculate how many momentum vectors to swim
C to get to the right Z:
C First Line:
        ZVZSWIM1 = (XSVI(1) - ZVX1)/ZVPX1
C Second Line:
        ZVZSWIM2 = (XSVI(1) - ZVX2)/ZVPX2
C And calculate the two z points:
        ZZ(1,1) = ZVZ1 + ZVZSWIM1*ZVPZ1
        ZZ(2,1) = ZVZ2 + ZVZSWIM2*ZVPZ2
        J = 1
        RVI(J) = SQRT(XSVI(J)**2+YSVI(J)**2)
        DZSVI(J) = ZZ(2,J) - ZZ(1,J)
C Jump to end:
        GO TO 45
      ENDIF
C loop over the two vertex-tracks
      DO I=1,NTRK
C Track bank number
        J = IABS(LIST(JJ(I)))
C dummy out MASS,RADIUS
        A = 0.0
        CALL CGETTRK (J,TKBANK(JJ(I)), A,A,0
     >,              K,P,Q, IERR)
C Quit, if error finding track data
        IF (IERR .NE. 0) THEN
          TKERR(JJ(I)) = IJKERR(2)
          GO TO 1000
        END IF
        DO K=1,5
          HELIX(MAP(K),I) = P(K)
        END DO
      END DO
C  For convenience, order the circles so that the first has smaller radius...
      II(1) = 1
      II(2) = 2
      IF (NTRK.EQ.2) THEN
        IF (ABS(HELIX(1,1)).LT.ABS(HELIX(1,2))) THEN
          II(1) = 2
          II(2) = 1
        END IF
      END IF
C make the circle radii, centers
      DO I=1,NTRK
C there are two tracks (circles)
        RC(I) = 0.5D0/HELIX(1,II(I))
        U     = RC(I)+HELIX(4,II(I))
        XC(I) =-U*SIN(HELIX(2,II(I)))
        YC(I) = U*COS(HELIX(2,II(I)))
        RA(I) = ABS(RC(I))
      END DO
C  there are usually two intersections
      DO J=1,2
        XSVI(J) = 0.0
        YSVI(J) = 0.0
        RVI (J) = 2.0*RVMAX
        DZSVI(J)= 2.0*DZMAX
        DO I=1,2
          SS(I,J) = 0.0
          ZZ(I,J) = 0.0
          TANG(I,J) = PI
        END DO
      END DO
      IF (NTRK.GT.1) GOTO 10
C Special handling of single-track vertex
      TMP = 1.0/(1.0+2.0*HELIX(1,1)*HELIX(4,1))
C Track parametrization
      FITA = HELIX(1,1)*TMP
      FITB = HELIX(4,1)*(1+HELIX(1,1)*HELIX(4,1))*TMP
      XYP0 = CMPLX(COS(HELIX(2,1)),SIN(HELIX(2,1)))
C Daughter vertex momentum
      PXYV = CMPLX(VTXP4(1,JJ(2)),VTXP4(2,JJ(2)))
      TMP = 1.0/ABS(PXYV)
      PXYV = PXYV*TMP
      PZV = VTXP4(3,JJ(2))*TMP
C Position & momentum, rotated to PHI0=0 coordinate system
      XYVR = CMPLX(XYZVRT(1,JJ(2)),XYZVRT(2,JJ(2)))*CONJG(XYP0)
      PXYVR = PXYV*CONJG(XYP0)
C Solve quadratic equation for displacement DD from vertex to track
      AA = FITA
      BB = 0.5*IMAG(PXYVR)-FITA*REAL(XYVR*CONJG(PXYVR))
      CC = FITA*ABS(XYVR)**2+FITB-IMAG(XYVR)
      TMP = BB*BB-AA*CC
      IF (TMP.GT.0.) THEN
        NSOL = 2
        IF (BB.LT.0.) THEN
          AA = -AA
          BB = -BB
          CC = -CC
        ENDIF
        BB = BB+SQRT(TMP)
        DD = CC/BB
        DO I=1,2
          XSVI(I) = XYZVRT(1,JJ(2))+DD*REAL(PXYV)
          YSVI(I) = XYZVRT(2,JJ(2))+DD*IMAG(PXYV)
          ZZ(2,I) = XYZVRT(3,JJ(2))+DD*PZV
          TANG(2,I) = 0.0
          DD = BB/AA
        ENDDO
        GOTO 20
      ENDIF
C No intersection; find closest approach and use half the distance
      NSOL = 1
C Check for line and circle antiparallel at closest approach
      IF (REAL(PXYVR).LT.0) THEN
        PXYV = -PXYV
        PXYVR = -PXYVR
        PZV = -PZV
      ENDIF
C Use SIN(PHI)**2/(1+COS(PHI)) for (1-COS(PHI)) to reduce roundoff
      CC = IMAG(PXYVR)**2/(1.0+REAL(PXYVR))
C Store 1/(2c) to avoid recalculation
      AA = 0.5/HELIX(1,1)
C TMP is signed distance from circle to line, along (-SIN(PHI),COS(PHI))
      TMP = (IMAG(XYVR)-HELIX(4,1))*REAL(PXYVR)
     &  -REAL(XYVR)*IMAG(PXYVR)+CC*AA
      XYTMP = CMPLX(IMAG(PXYV),CC)*AA
     &  +CMPLX(0.0,HELIX(4,1))+0.5*TMP*CMPLX(-IMAG(PXYV),REAL(PXYV))
      DD = REAL((XYTMP-XYVR)*CONJG(PXYVR))
      XYTMP = XYTMP*XYP0
      XSVI(1) = REAL(XYTMP)
      YSVI(1) = IMAG(XYTMP)
      ZZ(2,1) = XYZVRT(3,JJ(2))+DD*PZV
      TANG(2,1) = 0
      GOTO 20
C two track (or more) vertex
   10 CONTINUE
C get the circle center separation
      DX = XC(2) - XC(1)
      DY = YC(2) - YC(1)
      D  = DX*DX+DY*DY
C the circles are concentric
      IF (D.LE.0.D0) THEN
        IJKERR(2) = 1
        GO TO 1000
      END IF
      D = DSQRT(D)
C Separation (signed quantity);  if >0, the
      U = D-RA(1)-RA(2)
C  two circles do not intersect,    <0, they may
      FMCDIF(1) = U
C Special handling of conversion vertex
      IF (CONV) THEN
C the circles are too far appart to accept
        IF (ABS(U).GT.DRMAX) THEN
          IF (U.GT.0.0) THEN
C ...exteriour
            IJKERR(2) = 2
          ELSE
C ...interiour
            IJKERR(2) = 3
          END IF
          GO TO 1000
        END IF
        NSOL = 1
C  Vertex is track radius plus half of separation away from each center
        XSVI(1) = (XC(1)*(RA(2)+0.5*U)+XC(2)*(RA(1)+0.5*U))/D
        YSVI(1) = (YC(1)*(RA(2)+0.5*U)+YC(2)*(RA(1)+0.5*U))/D
C Branch down to vertex acceptability checking
        GO TO 20
      END IF
C Rotate, translate to a system where the
      COST = DX/D
C   two circle centers lie on the X axis.
      SINT = DY/D
C Y` = Y1` = Y2`
      Y0 = (-XC(1)*YC(2)+XC(2)*YC(1))/D
C X1', X2'
      DO I=1,2
        AB(I) = COST*XC(I) + SINT*YC(I)
      END DO
      U  = (XC(2)+XC(1))*(XC(2)-XC(1)) + (YC(2)+YC(1))*(YC(2)-YC(1))
      V  = (RA(2)+RA(1))*(RA(2)-RA(1))
C the common circle X` (if they intersect)
      XX = 0.5D0 * (U-V)/D
      U  = DSQRT((XX-AB(1))**2)
C Y''**2 is positive if they intersect
      YY2= (RA(1)+U) * (RA(1)-U)
C   the circles intersect;
      IF (YY2.GT.0.0) THEN
C    two intersection points (+/- Y)
        YY = DSQRT(YY2)
C     invert the translation, rotation
        DO J=1,2
          U = YY+Y0
          XSVI(J) = COST*XX - SINT*U
          YSVI(J) = SINT*XX + COST*U
          YY =-YY
        END DO
        NSOL = 2
        GO TO 20
      END IF
C       We get here if the two circles do not intersect;
C       Find how close in the XY plane they approach each other,
C       and take the point half way between them as our vertex approximation.
      U = D - (RA(1)+RA(2))
C A is outside of B
      IF (U.GT.0.0D0) THEN
        V = U
        J = 2
        IF (AB(1).LT.AB(2)) J=1
        XX = AB(J) + RA(J)
C A is inside of B
      ELSE
        IF (AB(1).LT.AB(2)) THEN
          XX = AB(2) - RA(2)
          V  = AB(1) - RA(1) - XX
        ELSE
          XX = AB(1) + RA(1)
          V  = AB(2) + RA(2) - XX
        END IF
      END IF
      XX = XX + 0.5D0*V
C rotate back to the original system
      XSVI(1) = COST*XX - SINT*Y0
      YSVI(1) = SINT*XX + COST*Y0
      NSOL = 1
C the circles are too far appart to accept
      IF (V.GT.DRMAX) THEN
        IF (U.GT.0.0D0) THEN
C ...exteriour
          IJKERR(2) = 4
        ELSE
C ...interiour
          IJKERR(2) = 5
        END IF
        GO TO 1000
      END IF
   20 CONTINUE
C loop over solutions
      DO 30 J=1,NSOL
C loop over tracks
        DO I=1,NTRK
C point from circle center
          U = (XSVI(J)-XC(I))/RC(I)
C   to the intersection vertex
          V =-(YSVI(J)-YC(I))/RC(I)
C turning angle from the track origin
          U = ATAN2(U,V) - HELIX(2,II(I))
          IF      (U.LT.-PI) THEN
            U = U + TWOPI
          ELSE IF (U.GT. PI) THEN
            U = U - TWOPI
          END IF
          TANG(I,J) = U
C arc length
          SS(I,J) = RC(I)*U
          ZZ(I,J) = HELIX(5,II(I)) + SS(I,J)*HELIX(3,II(I))
        END DO
        RVI(J) = SQRT(XSVI(J)**2+YSVI(J)**2)
        DZSVI(J) = ZZ(2,J) - ZZ(1,J)
   30 CONTINUE
C  Check that there is at least one potentially acceptable solution!
      A = AMIN1(RVI(1),RVI(2))
C check the vertex radius is acceptable
      IF (A.GT.RVMAX) THEN
        IJKERR(2) = 6
        GO TO 1000
      END IF
      A = AMIN1(ABS(DZSVI(1)),ABS(DZSVI(2)))
C check the Z difference is acceptable
      IF (A.GT.DZMAX) THEN
        IJKERR(2) = 7
        GO TO 1000
      END IF
      A = AMAX1(ABS(TANG(1,1)),ABS(TANG(2,1)))
      B = AMAX1(ABS(TANG(1,2)),ABS(TANG(2,2)))
C check there is an acceptable TANG
      IF (AMIN1(A,B).GT.TRNMAX) THEN
        IJKERR(2) = 8
        GO TO 1000
      END IF
C minimum track arc length, sol`n 1
      A = AMIN1(SS(1,1),SS(2,1))
C minimum track arc length, sol`n 2
      B = AMIN1(SS(1,2),SS(2,2))
C limit the minimum arc length
      IF (AMAX1(A,B).LT.DSMIN) THEN
        IJKERR(2) = 9
        GO TO 1000
      END IF
C  there may be a possible acceptable solution, in (R,deltaZ,S,Tang)
      J = 1
      IF (NSOL.EQ.1) GO TO 40
      IF (ABS(DZSVI(2)).LT.ABS(DZSVI(1))) J=2
      FLIP = 0
      A = AMAX1(ABS(TANG(1,J)),ABS(TANG(2,J)))
      B = AMIN1(SS(1,J),SS(2,J))
      IF (RVI(J).GT.RVMAX)        FLIP = 1
      IF (ABS(DZSVI(J)).GT.DZMAX) FLIP = FLIP+2
      IF (A.GT.TRNMAX)            FLIP = FLIP+4
      IF (B.LT.DSMIN)             FLIP = FLIP+8
      IF (FLIP.NE.0) THEN
        IF (J.EQ.1) THEN
          J = 2
        ELSE
          J = 1
        END IF
      END IF
C  final checks of the final solution...
   40 CONTINUE
      A = AMAX1(ABS(TANG(1,J)),ABS(TANG(2,J)))
      B = AMIN1(SS(1,J),SS(2,J))
C check that TANG is acceptable
      IF (A.GT.TRNMAX) THEN
        IJKERR(2) = 8
        GO TO 1000
      END IF
C limit the minimum arc length
      IF (B.LT.DSMIN) THEN
        IJKERR(2) = 9
        GO TO 1000
      END IF
C Continue point from zero track case:
   45 CONTINUE
C check the vertex radius
      IF (RVI(J).GT.RVMAX) THEN
        IJKERR(2) = 6
        GO TO 1000
      END IF
C check the Z difference
      IF (ABS(DZSVI(J)).GT.DZMAX) THEN
        IJKERR(2) = 7
        GO TO 1000
      END IF
C  "acceptable" solution.  Set initial vertex approximation
      XYZVRT(1,Nv) = XSVI(J)
      XYZVRT(2,Nv) = YSVI(J)
      XYZVRT(3,Nv) = 0.5 * (ZZ(1,J) + ZZ(2,J))
C  Go off and collect (GETTRK) the track data, and make the first approximation
C  vertex momentum sums (which will be needed for single track vertex first
C  approximation calculations).
C  GETTRK failures or failures in transporting the tracks not used in the vertex
C  finding to the vertex  will return IERR = IJKERR(1) = 3.
   50 CONTINUE
D     IF (DBGPRT.GT.0) THEN
D       RR = SQRT(XYZVRT(1,Nv)**2+XYZVRT(2,Nv)**2)
D       WRITE (PRINT,1010) Nv,(XYZVRT(I,Nv),I=1,3),RR
D1010   FORMAT (/,' Vertex',I2,
D    &      ' Approximation; ',3F10.4,5X,'RR',F7.2)
D     END IF
C Initialize vertex 3-momentum to the sum
      VtxP4(1,Nv) = 0
C of daughter vertex 3-momenta
      VtxP4(2,Nv) = 0
C (needed for 1-track vertex CTVMFA)
      VtxP4(3,Nv) = 0
      DO Mv=Nv+1,NVERTX
        IF (VtxPnt(Mv,1).EQ.Nv) THEN
          VtxP4(1,Nv) = VtxP4(1,Nv)+VtxP4(1,Mv)
          VtxP4(2,Nv) = VtxP4(2,Nv)+VtxP4(2,Mv)
          VtxP4(3,Nv) = VtxP4(3,Nv)+VtxP4(3,Mv)
        END IF
      END DO
C Collect track information for this vertex
      SUMZOS = 0.0
      SUMOOS = 0.0
      DO Nt=1,NTRACK
        IF (TrkVtx(Nt,Nv)) THEN
          CALL CCTVM01(KPRINT,Nt,BEAM,PARERR(1,Nt),IERR)
          IF (IERR.NE.0) THEN
C GETTRK failure, or, track cannot reach vertex
            GO TO 1000
          END IF
          TKERR(Nt) = 0
          SPh0 = SIN(Par0(2,Nt))
          CPh0 = COS(Par0(2,Nt))
          Xsvt = XYZVRT(1,Nv)*CPh0+XYZVRT(2,Nv)*SPh0
          Ysvt = XYZVRT(2,Nv)*CPh0-XYZVRT(1,Nv)*SPh0
C Proceed blissfully
          SPhs = 2*Par0(1,Nt)*Xsvt
C 10/23/02 Add protection here against unphysical Sphs.
C Craig Blocker
          Sphs = max(-1.,min(1.,SPhs))
          CPhs = SQRT((1-SPhs)*(1+SPhs))
          SPhi = CPh0*Sphs+SPh0*CPhs
          CPhi = CPh0*CPhs-SPh0*SPhs
C add this track`s contribution to vertex momentum
          PT = PSCALE/ABS(Par0(1,Nt))
          VtxP4(1,Nv) = VtxP4(1,Nv)+PT*CPhi
          VtxP4(2,Nv) = VtxP4(2,Nv)+PT*SPhi
          VtxP4(3,Nv) = VtxP4(3,Nv)+PT*Par0(3,Nt)
C increment sums for average z0. Use only track with z measurement
          IF(PARERR(5,NT).GT.0.0.AND.VTXPNT(1,1).EQ.-100)THEN
             SUMZOS = SUMZOS + PAR0(5,NT)/(PARERR(5,NT)**2)
             SUMOOS = SUMOOS + 1.0/(PARERR(5,NT)**2)
D            IF (DBGPRT.GT.0) THEN
D             write(print,5161)NT, (par0(i,NT),i=1,5)
D5161        format(1x,'track ',i3,' Par: C,phi,ctg,D,Z',5(1x,e10.4))
D             write(print,5162)NT, (parerr(i,NT),i=1,5)
D5162        format(1x,'track ',i3,' Err: C,phi,ctg,D,Z',5(1x,e10.4))
D            END IF
          END IF
        END IF
      END DO
C Get z average to get good starting point in beamline fit
      IF (VTXPNT(1,1).EQ.-100) THEN
         ZAVG = 0.0
         IF (SUMOOS.GT.0.0)THEN
            ZAVG = SUMZOS/SUMOOS
            XYZVRT(3,NV) = ZAVG
         END IF
      END IF
 1000 CONTINUE
      IF (IERR.NE.0) THEN
        IJKERR(1) = IERR
      ELSE IF (IJKERR(2).NE.0) THEN
        IERR      = 2
        IJKERR(1) = 2
        IJKERR(3) = Nv
      END IF
      RETURN
      END
C==============================================================================
      SUBROUTINE CCTVM01 (PRINT, NT,BEAM, PARERR, IERR)
C==============================================================================
C===Description:
C   Collects/checks track data for the use of CTVMFT for track NT (uses GETTRK),
C   with bank number LIST(NT), bank type TKBANK(NT), and mass TMASS(NT).
C   The track/vertex configuration TrkVtx and vertex approximation
C   XYZVRT(xyz,NV) are communicated through the include file CTVMFT.INC.
C   The vertex approximation is used to find the RADIUS at which dE/dX and
C   Coulomb multiple scattering contributions are evaluated.
C===Input Arguments:
C   PRINT    If compiled DEBUG, outputs formatted output to unit PRINT
C   NT       Desired track, in the array LIST
C   BEAM     Beam constraint flag, for GETTRK
C===Output Arguments:
C   PARERR   Track helix parameter errors
C   IERR     = IJKERR(1) = 3   flags an error getting the track parameters
C              IJKERR(2) = 1   GETTRK returns an error for this track
C                          2   the track covariasnce matrix is uninvertable
C                          3   turns through too large an angle to the vertex
C                          4   moves too far backwards to the vertex
C===Author:
C   see CTVMFT
C-------------------------------------------------------------------------------
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmco.h"
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
C===Local Declarations:
      INTEGER NT, PRINT, IERR
      REAL    PARERR(5)
      LOGICAL BEAM
      INTEGER I,J,K,L, NV
      REAL    P(5),Q(5,5), HELIX(5,MaxTrk),HWGT(5,5,MaxTrk)
      REAL    RADV, TS,S
      REAL    WORK(MAXDIM)
      REAL*8  RC,XYC(2),U,V, ELM(5,5)
      CHARACTER*80 STRING
C Map the parameter order to CTVMFT form:
      INTEGER MAP(5)
C (Ctg,Crv,Z0,D0,Phi)->(Crv,Phi,Ctg,D0,Z0)
      DATA    MAP / 3,1,5,4,2 /
C===Start of Code:
C-------------------------------------------------------------------------------
C get the vertex for this track
        DO NV=1,NVERTX
          IF (TrkVtx(NT,NV)) GO TO 11
        END DO
   11   RADV = SQRT(XYZVRT(1,NV)**2 + XYZVRT(2,NV)**2)
C call GETTRK for track data
        IF (LIST(NT).GT.0) THEN
          CALL CGETTRK (LIST(NT),TKBANK(NT),TMASS(NT),RADV, BEAM
     >                 ,K,P,Q, IERR)
C Quit on error finding track data bank
          IF (IERR.NE.0) THEN
            TKERR(NT) = IERR
            IJKERR(2) = 1
            GO TO 100
          END IF
C save the parameter vector
          DO I=1,5
            K = MAP(I)
            HELIX(K,NT) = P(I)
            PARERR(K) = SQRT(Q(I,I))
            DO J=1,5
              L = MAP(J)
              ELM(L,K) = Q(J,I)
            END DO
          END DO
C make -and save- the weight matrix
          CALL DINV(5,ELM,5,WORK,IERR)
          IF (IERR.NE.0) THEN
            TKERR(NT) = IERR
            IJKERR(2) = 2
            GO TO 100
          END IF
          DO I=1,5
          DO J=1,5
            HWGT(J,I,NT) = ELM(J,I)
          END DO
          END DO
        END IF
C load the fit parameter vector,
        DO I=1,5
          PAR0(I,NT) = HELIX(I,NT)
C and the corresponding weight matrix
          DO J=1,5
            G(J,I,NT) = HWGT(J,I,NT)
          END DO
        END DO
        RC = 0.5D0/HELIX(1,NT)
        U     = RC + HELIX(4,NT)
        XYC(1)=-U*SIN(HELIX(2,NT))
        XYC(2)= U*COS(HELIX(2,NT))
        U = (XYZVRT(1,NV)-XYC(1))/RC
        V =-(XYZVRT(2,NV)-XYC(2))/RC
        TS = ATAN2(U,V) - HELIX(2,NT)
        IF      (TS.LT.-PI) THEN
          TS = TS + TWOPI
        ELSE IF (TS.GT. PI) THEN
          TS = TS - TWOPI
        END IF
        S = RC*TS
C it turns too much
        IF (ABS(TS).GT.TRNMAX) THEN
          IJKERR(2) = 3
C ?what to do?
        ELSE IF (S.LT.DSMIN) THEN
          IJKERR(2) = 4
        END IF
D       IF (PRINT.GT.0) THEN
D       IF (IJKERR(2).NE.0) THEN
D         WRITE(STRING,1023) NT,LIST(NT), S
D1023 FORMAT(' negative arc length, track',I2,I6 ,F7.2)
D       END IF
D       IF (PRINT.GT.0) THEN
D         IF (NT.EQ.1) WRITE(PRINT,1025)
D         WRITE(PRINT,1026) NT,(HELIX(I,NT),I=1,5), S
D1025 FORMAT(/,' Track  Initial helix parameters,  arc length to Vtx')
D1026 FORMAT(I6,2X,1P,E11.3,0P,F10.4,3F10.4, F12.2)
D       END IF
D       END IF
C===Return to Caller:
  100 CONTINUE
      IF (IJKERR(2).NE.0) THEN
        IERR = 3
        IJKERR(1) = 3
        IJKERR(3) = NT
      END IF
      RETURN
      END
C==============================================================================
      SUBROUTINE CCTVMVF (PRINT, Nv,  VXI)
C==============================================================================
C===Description:
C   Calculates contributions from track Nt (at vertex Nv) to the derivative
C   vector VXI and derivitive matrix VMAT for the vertex fit section of CTVMFT.
C  This portion of the code is somewhat awkward to follow.  See the derivation
C  of these formulae in CDF note 19zz (J.P. Marriner)
C  NB:  this code does not properly handle tracks that turn through more than
C       90 degrees between the distance of closest approach to the z axis and
C       the vertex.  These tracks are considered pathological.  The code will
C       be correct everywhere, however, if the proper branch of the arcsin
C       function is taken.  (The initial approximation may be so bad that the
C       desired solution will not be found, however).
C===Input Arguments:
C   Nv      vertex
C===Output Arguments:
C   VXI     first derivative contributions for this track to the vertex fit
C   VMAT    second derivative matrix (in common, include file CTVMFT.INC)
C------------------------------------------------------------------------ ------
C===Implicit None Declaration:
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmco.h"
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
C===Local Declarations:
      INTEGER PRINT
      INTEGER NV
      REAL*8  VXI(MAXDIM)
      INTEGER I,J,K,L
      INTEGER  NT, NtF,NvF, TI,TJ, VI,Vj
      REAL    C,PHI,COTH, D,Z, PT
      REAL    TWOC, S,TS,TRAD
      REAL    SP,CP, XCYS,XSYC,XC2YS2, SINCS,COSCS,TWOCS
      REAL    DPXDT, DPYDT, DSDC,DSDP,DSDX,DSDY
      REAL    EMAT(3,2), FMAT(3,2)
C------------------------------------------------------------------------ ------
C===Start of Code:
      DO 50 Nt=1,NTRACK
      IF (.NOT.TrkVtx(Nt,Nv)) GO TO 50
C current half curvature
      C    = PAR0(1,NT) + PARDIF(1,NT)
C current phi
      PHI  = PAR0(2,NT) + PARDIF(2,NT)
C current coth(theta)
      COTH = PAR0(3,NT) + PARDIF(3,NT)
      TWOC = 2.0*C
C        Xsv*cos(phi) + Ysv*sin(phi),  -Xsv*sin(phi) + Ysv*cos(phi)
C      (Xsv*cos(phi))**2 + (Ysv*sin(phi))**2
      SP     = SIN(PHI)
      CP     = COS(PHI)
      XCYS   =  XYZVRT(1,Nv)*CP + XYZVRT(2,Nv)*SP
      XSYC   = -XYZVRT(1,Nv)*SP + XYZVRT(2,Nv)*CP
      XC2YS2 = (XYZVRT(1,Nv)*CP)**2 + (XYZVRT(2,Nv)*SP)**2
C       t = 2.0 * c * (Xsv*cos(phi) + Ysv*sin(phi)) = argument of arcsin
C      s = projected distance in xy plane along helix to vertex
C       sin ( c*s), cos ( c*s ), 2.0 * sin(c*s) * cos(c*s)
C      d(arcsin(t))/dt = 1 / sqrt (1 - t**2 )
      TS    = TWOC*XCYS
C 1/20/04 Protect against TWOC begin zero.  Craig Blocker
      if(twoc .eq. 0.) then
        s = xcys
      else
        S = ASIN(TS)/TWOC
      endif
      SINCS = SIN(C*S)
      COSCS = COS(C*S)
      TWOCS = 2.0*SINCS*COSCS
      TRAD  = SQRT((1.-TS)*(1.+TS))
C       Helix parameters d0 and z0 are dependent variables.
C       They are functions of the vertex and track parameters.
C 1/20/04 Protect against C being zero. Craig Blocker
      if(c .eq. 0.) then
        D = xsyc
      else
        D = XSYC - SINCS**2/C
      endif
      Z = XYZVRT(3,Nv) - COTH*S
C       Difference between fitted and measured d0, z0
      PARDIF(4,NT) = D - PAR0(4,NT)
      PARDIF(5,NT) = Z - PAR0(5,NT)
C      Remember all flavors of momentum for this track
C 1/20/04 Protect against C being zero.  Here it is not
C clear what to set PT to.  I chose 999. GeV since this is
C a very large Pt in CDF, this is roughly a standard deviation
C from infinite, and should be clearly recognizable.  Craig Blocker
        if(c .eq. 0.) then
          PT = 999.
        else
          PT = PSCALE/ABS(C)
        endif
C x component of momentum
      TrkP4(NT,1) = PT*(CP*TRAD-SP*TS)
C y component of momentum
      TrkP4(NT,2) = PT*(SP*TRAD+CP*TS)
C z component of momentum
      TrkP4(NT,3) = PT*COTH
C Transverse momentum
      TrkP4(NT,5) = PT
C total momentum
      TrkP4(NT,6) = PT*SQRT(1.0+COTH**2)
C energy
      TrkP4(NT,4) = SQRT(TrkP4(NT,6)**2+TMASS(NT)**2)
C dPx/dt
      DPXDT =-PT*(SP+TS*CP/TRAD)
C dPy/dt
      DPYDT = PT*(CP-TS*SP/TRAD)
C dPx/dc
      if(c.ne.0.) then
        DDA(NT,1) =-TrkP4(NT,1)/C + 2.0*XCYS*DPXDT
      else
        dda(nt,1) = 999.
      endif
C dPx/dPhi
      DDA(NT,2) =-PT*(TRAD*SP+TS*CP) + TWOC*XSYC*DPXDT
C dPx/dXsv
      DDA(NT,3) = TWOC*CP*DPXDT
C dPx/dYsv
      DDA(NT,4) = TWOC*SP*DPXDT
C dPy/dc
      if(c .ne. 0.) then
        DDA(NT,5) =-TrkP4(NT,2)/C + 2.0*XCYS*DPYDT
      else
        dda(nt,5) = 999.
      endif
C dPy/dPhi
      DDA(NT,6) = PT*(TRAD*CP-TS*SP) + TWOC*XSYC*DPYDT
C dPy/dXsv
      DDA(NT,7) = TWOC*CP*DPYDT
C dPy/dYsv
      DDA(NT,8) = TWOC*SP*DPYDT
C      Get the derivatives of s with respect to the fit parameters
C         See if argument of arcsin is small
C         (track with a small turning angle to the vertex point)
C Yes. use power series approximation
      IF (ABS(TS).LT.5.E-3) THEN
        DSDC =-(0.666667 - 0.3*TS**2)*TS**3
C No.  Use full functional form
      ELSE
        DSDC = TS/TRAD - TWOC*S
      ENDIF
C Get full derivative ds/dc
      if(twoc .ne. 0.) then
        DSDC = DSDC/(TWOC*C)
      else
        dsdx = 999.
      endif
C ds/dphi
      DSDP = XSYC/TRAD
C ds/dXsv
      DSDX = CP/TRAD
C ds/dYsv
      DSDY = SP/TRAD
C      d(d0)/dc, d(d0)/dphi, d(d0)/dcot(theta)
      if(c .ne. 0.) then
        EMAT(1,1) = -SINCS*(TWOC*S*COSCS-SINCS)/C**2 - TWOCS*DSDC
      else
        emat(1,1) = 999.
      endif
      EMAT(2,1) = -XCYS - TWOCS*DSDP
      EMAT(3,1) = 0.0
C      d(z0)/dc, d(z0)/dphi, d(z0)/dcot(theta)
      EMAT(1,2) = -COTH*DSDC
      EMAT(2,2) = -COTH*DSDP
      EMAT(3,2) = -S
C      d(d0)/dXsv, d(d0)/dYsv, d(d0)/dZsv
      FMAT(1,1) = -SP - TWOCS*CP/TRAD
      FMAT(2,1) =  CP - TWOCS*SP/TRAD
      FMAT(3,1) = 0.0
C      d(z0)/dXsv, d(z0)/dYsv, d(z0)/dZsv
      FMAT(1,2) = -COTH*DSDX
      FMAT(2,2) = -COTH*DSDY
      FMAT(3,2) = 1.0
C Index into matrix equations for track NT
      NtF = TOFF(NT)
C Index into matrix equations for vertex NV
      NvF = VOFF(Nv)
      DO 20 I=1,3
        TI = NtF + I
        VI = NvF + I
C loop over (d0,z0)
        DO K=1,2
C            Ft * G2 * Xi2
          VXI(Vi) = VXI(Vi) - FMAT(I,K)* (G(K+3,4,NT)*PARDIF(4,NT)
     &                                  + G(K+3,5,NT)*PARDIF(5,NT))
C           (Gd + Et * G2 ) * Xi2
          VXI(Ti) = VXI(Ti) - PARDIF(K+3,NT) *
     &    (G(I,K+3,NT) + EMAT(I,1)*G(4,K+3,NT) + EMAT(I,2)*G(5,K+3,NT))
        END DO
        DO 10 J=1,3
          TJ = NtF + J
          VJ = NvF + J
C            Ft * Gd * Xi3
          VXI(Vi) = VXI(Vi)
     &    - (FMAT(I,1)*G(4,J,NT)+FMAT(I,2)*G(5,J,NT))*PARDIF(J,NT)
C            (G3 + Et * Gd) * Xi3
          VXI(Ti) = VXI(Ti) - (G(I,J,NT)
     &    + EMAT(I,1)*G(4,J,NT)+EMAT(I,2)*G(5,J,NT))*PARDIF(J,NT)
C            G3
          VMAT(Ti,Tj) = G(I,J,NT)
          DO K=1,2
C            Gd * F
            VMAT(Vi,Tj) = VMAT(Vi,Tj)+FMAT(I,K)*G(K+3,J,NT)
C            Et * Gdt + Gd * E
            VMAT(Ti,Tj) = VMAT(Ti,Tj)
     &      + EMAT(I,K)*G(K+3,J,NT)+ G(I,K+3,NT)*EMAT(J,K)
            DO L=1,2
C            Ft * G2 * F
              VMAT(Vi,Vj) = VMAT(Vi,Vj)
     &        + FMAT(I,K)*G(K+3,L+3,NT)*FMAT(J,L)
C            Et * G2 * F
              VMAT(Vi,Tj) = VMAT(Vi,Tj)
     &        + FMAT(I,K)*G(K+3,L+3,NT)*EMAT(J,L)
C            Et * G2 * E
              VMAT(Ti,Tj) = VMAT(Ti,Tj)
     &        + EMAT(I,K)*G(K+3,L+3,NT)*EMAT(J,L)
C end of loop on L
            END DO
C end of loop on K
          END DO
C end of loop on J
   10   CONTINUE
C end of loop on I
   20 CONTINUE
D     IF (PRINT.GT.0)  THEN
D     TS = TS/C
D     WRITE(PRINT,1045) NT,(PARDIF(I,NT),I=1,5), D,Z,TS
D1045 FORMAT(/,' Track',I3,1X,1P,5E11.3,3X,0P,3F11.4)
D     WRITE(PRINT,1047) (DDA(Nt,J),J=1,8)
D1047 FORMAT(10X,1P,8E11.3)
D     WRITE(PRINT,1046) (VXI(I),I=1,MATDIM)
D     DO I=1,MATDIM
D     WRITE(PRINT,1046) (VMAT(J,I),J=1,I)
D1046 FORMAT(5X,0P,10E11.3)
D     END DO
D     END IF
C end of the track loop, this vertex
   50 CONTINUE
C------------------------------------------------------------------------ ------
C===Return to Caller:
C Symmetrize the derivative matrix
      DO I=1,MATDIM-1
      DO J=I+1,MATDIM
         VMAT(J,I) = VMAT(I,J)
      END DO
      END DO
      RETURN
      END
C===============================================================================
      SUBROUTINE CCTVMCF (PRINT, NV,  VXI)
C===============================================================================
C  CTVMCF calculates "conversion constraint" contributions to the derivative
C         vector and matrix for vertex Nv (08/29/94 JPB and WJA)
C--Input parameters
C        Input data in the include file CTVMFT (COMMON /CTVMFC/)
C  A physically reasonable constraint to impose in fitting a vertex at which
C  a photon conversion is believed to have occurred is that the invariant
C  mass of the two tracks be as small as possible, i.e. that the opening angle
C  of the two tracks at the vertex be zero.  Where the tracks are copunctual,
C  they should be collinear.
C  In r-phi, a convenient way to express this constraint is that at the radius
C  of conversion, the sum of the azimuth and the turning angle is equal for
C  both tracks, i.e. Phi1+2*Crv1*S1(r) = Phi2+2*Crv2*S2(r), where S_i(r) is
C  the arc length (projected into r-phi) of track i at radius r, and Phi_i
C  and Crv_i are the azimuth and half-curvature of track i.  In other words,
C  we equate the r-phi momentum directions of the two tracks at the vertex.
C  In r-z, the constraint is clearly that Ctg1 = Ctg2.
C-------------------------------------------------------------------------------
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmco.h"
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
      INTEGER PRINT
      INTEGER Nv
      REAL*8  VXI(MaxDim)
      INTEGER Trk1,Trk2,NvF,NcF,I,J,IT1F,IT2F
      REAL    XV,YV,Crv1,Crv2,Phi1,Phi2,Ctg1,Ctg2,Phs1,Phs2,S1,S2
      REAL    DS1DX,DS2DX,DS1DY,DS2DY,DSDCrv1,DSDCrv2,DSDPhi1,DSDPhi2
      REAL    CPhi1,CPhi2,SPhi1,SPhi2,CPhs1,CPhs2,SPhs1,SPhs2
      REAL    SHOUT,MINUS/-1/,ZERO/0/
      CHARACTER*80 STRING
C-------------------------------------------------------------------------------
C  Find the first two tracks for this vertex
      Trk1 = 0
      Trk2 = 0
      DO I=1,NTRACK
        IF (TrkVtx(I,Nv)) THEN
          IF (Trk1.EQ.0) THEN
            Trk1 = I
          ELSE
            Trk2 = I
            GO TO 100
          END IF
        END IF
      END DO
 100  CONTINUE
C  Save some possibly useful offsets into VMAT
      IT1F = TOFF(Trk1)
      IT2F = TOFF(Trk2)
      NvF = VOFF(Nv)
      NcF = COFF(Nv)
C  Vertex position
      XV = XYZVRT(1,Nv)
      YV = XYZVRT(2,Nv)
C  Track parameters
      Crv1 = PAR0(1,Trk1)+PARDIF(1,Trk1)
      Crv2 = PAR0(1,Trk2)+PARDIF(1,Trk2)
      Phi1 = PAR0(2,Trk1)+PARDIF(2,Trk1)
      Phi2 = PAR0(2,Trk2)+PARDIF(2,Trk2)
      Ctg1 = PAR0(3,Trk1)+PARDIF(3,Trk1)
      Ctg2 = PAR0(3,Trk2)+PARDIF(3,Trk2)
C  Simple functions of track parameters
      SPhi1 = SIN(Phi1)
      SPhi2 = SIN(Phi2)
      CPhi1 = COS(Phi1)
      CPhi2 = COS(Phi2)
C  Track turning angles (to get momentum direction)
      SPhs1 = 2*Crv1*(XV*CPhi1+YV*SPhi1)
      SPhs2 = 2*Crv2*(XV*CPhi2+YV*SPhi2)
      CPhs1 = SQRT((1+SPhs1)*(1-SPhs1))
      CPhs2 = SQRT((1+SPhs2)*(1-SPhs2))
      Phs1 = ASIN(SPhs1)
      Phs2 = ASIN(SPhs2)
C  Arc length in XY plane
      S1 = 0.5*Phs1/Crv1
      S2 = 0.5*Phs2/Crv2
C  d(arc length)/d(vertex position)
      DS1DX = CPhi1/CPhs1
      DS2DX = CPhi2/CPhs2
      DS1DY = SPhi1/CPhs1
      DS2DY = SPhi2/CPhs2
C  d(arc length)/d(track parameters)
      IF (ABS(CPhs1).GT.0.005) THEN
        DSDCrv1 = 0.5*(SPhs1/CPhs1-Phs1)/(Crv1*Crv1)
      ELSE
        DSDCrv1 = 0.33333333*SPhs1**3*(0.9*SPhs1*SPhs1-1)/(Crv1*Crv1)
      END IF
      IF (ABS(CPhs2).GT.0.005) THEN
        DSDCrv2 = 0.5*(SPhs2/CPhs2-Phs2)/(Crv2*Crv2)
      ELSE
        DSDCrv2 = 0.33333333*SPhs2**3*(0.9*SPhs2*SPhs2-1)/(Crv2*Crv2)
      END IF
      DSDPhi1 = (-XV*SPhi1+YV*CPhi1)/CPhs1
      DSDPhi2 = (-XV*SPhi2+YV*CPhi2)/CPhs2
C  Constraint in r-phi, C_rphi = (Phi1+2*Crv1*S1)-(Phi2+2*Crv2*S2)
      VXI(NcF+1) = (Phi2+Phs2)-(Phi1+Phs1)
      IF (VXI(NcF+1).LT.-PI) VXI(NcF+1) = VXI(NcF+1)+2D0*PI
      IF (VXI(NcF+1).GT.+PI) VXI(NcF+1) = VXI(NcF+1)-2D0*PI
C  Optional constraint in r-z, C_rz = Ctg1-Ctg2
      IF (Cvtx(Nv).EQ.2) THEN
        VXI(NcF+2) = Ctg2-Ctg1
      END IF
C  Derivatives of the r-phi constraint C_rphi with respect to the fit
C  parameters (X,Y,Z,Crv1,Phi1,Ctg1,Crv2,Phi2,Ctg2)
C d/dXv
      VMAT(NVF+1,NCF+1) = 2D0*(Crv1*DS1DX-Crv2*DS2DX)
C d/dYv
      VMAT(NVF+2,NCF+1) = 2D0*(Crv1*DS1DY-Crv2*DS2DY)
C d/dCrv1
      VMAT(IT1F+1,NCF+1) = +(1D0+2D0*(S1+Crv1*DSDCrv1))
C d/dCrv2
      VMAT(IT2F+1,NCF+1) = -(1D0+2D0*(S2+Crv2*DSDCrv2))
C d/dPhi1
      VMAT(IT1F+2,NCF+1) = +(1D0+2D0*Crv1*DSDPhi1)
C d/dPhi2
      VMAT(IT2F+2,NCF+1) = -(1D0+2D0*Crv2*DSDPhi2)
C  Optional derivatives of the r-z constraint C_rz with respect to
C  the fit parameters
      IF (Cvtx(Nv).EQ.2) THEN
C d/dCtg1
        VMAT(IT1F+3,NcF+2) = +1D0
C d/dCtg2
        VMAT(IT2F+3,NcF+2) = -1D0
      END IF
C  Symmetrize
      DO I=1,MATDIM-1
        DO J=I+1,MATDIM
          VMAT(J,I) = VMAT(I,J)
        END DO
      END DO
C  All done
      RETURN
      END
C===============================================================================
      SUBROUTINE CCTVMPF (PRINT, NV,  VXI)
C===============================================================================
C  CTVMPD calculates "pointing constraint" contributions to the derivative
C         vector and matrix for vertex Nv, "pointing" at vertex Mv =VtxPnt(Nv,2)
C--Input parameters
C        Input data in the include file CTVMFT (COMMON /CTVMFC/)
C  Let dX, dY, dZ = difference between the first and the second vertex
C  Px = VtxP4(1,Nv), etc
C  dX = DELV(1), etc
C  The first pointing constraint is  Px * dY - Py * dX = 0,
C  The second is                     Pi * dZ - Pz * di = 0
C                         (where i is x if |Px| > |Py|, otherwise i is y)
C  The first pointing constraint constrains the projections of the two vectors
C  P and d onto the x,y plane to be parallel (or antiparallel).
C  The second constraint constrains the vectors to be parallel or antiparallel
C  in 3 dimensions.
C  NB:  The program makes this explicit choice of axes so that one can apply an
C  x,y pointing constraint or (optionally) a 3 dimensional constraint.
C  This choice of axes will not produce the desired  results for the case
C  Px = Py = 0.  This case is thought to be unlikely; It might occur for massive
C  objects (Z`s, eg) produced with no transverse momentum.  For these types of
C  events one is probably better advised to use the "beam constraint" anyway.
C  Edit log:
C  ---- ---
C  10/xx/94  WJA  To handle 1-track vertex, add NPC variable and set it to
C                 the number of Lagrange multipliers needed to effect the
C                 desired pointing constraint: normally equal to VtxPnt(Nv,2)
C                 but equal to 2 if VtxPnt(Nv,2) is zero, since 1-track point
C                 points in both projections.
C  05/xx/02  MSM  Similar for zero-track vertexing.
C-------------------------------------------------------------------------------
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
      INTEGER PRINT
      INTEGER Nv
      REAL*8  VXI(MaxDim)
      REAL    TEMP(5), C, DELV(3)
      INTEGER I1,I2, IMAT,JMAT
      INTEGER NvF,MvF,Kv,KvF,NtF, LpF,Lp, I,J, Nt,Mv, Np
      INTEGER NPC
C-------------------------------------------------------------------------------
C   IMAT points to larger PSUM component (1=x,2=y)
      IF (ABS(VtxP4(1,Nv)).GT.ABS(VtxP4(2,Nv))) THEN
        IMAT = 1
        JMAT = 0
      ELSE
        IMAT = 2
        JMAT = 4
      END IF
C vertex towards which Nv "points"
      Mv = VtxPnt(Nv,1)
C vertex displacement vector
      DO I=1,3
        DELV(I) = XYZVRT(I,Mv) - XYZVRT(I,Nv)
      END DO
C offset to 2nd vertex
      NvF = VOFF(Nv)
C offset to this pointing constraint
      LpF = POFF(Nv)
C offset to target (primary) vertex
      MvF = 0
C the target vertex is NOT primary
      IF (Mv.GT.0) MvF=VOFF(Mv)
      NPC = VtxPnt(Nv,2)
C Use both constraints for 1-track vertex
      IF (NPC.EQ.3) NPC = 2
      IF (NPC.EQ.4) NPC = 2
C Loop over pointing constraints
      DO 50 Np=1,NPC
C offset to the current pointing constraint
        Lp = LpF + Np
C   First constraint.  1=x, 2=y
        IF (Np.EQ.1) THEN
          I1 = 1
          I2 = 2
C otherwise 2nd constraint
        ELSE
C 1=bigger of px,py (IMAT), 2=z
          I1 = IMAT
          I2 = 3
        ENDIF
C "residual" (Pi)
        VXI(Lp) =-VtxP4(I1,Nv)*DELV(I2) +VtxP4(I2,Nv)*DELV(I1)
C dPi/dI1p (1st vtx)
        VMAT(MvF+I1,Lp) =-VtxP4(I2,Nv)
C dPi/dI2p
        VMAT(MvF+I2,Lp) = VtxP4(I1,Nv)
C dPi/dI1s (2nd vtx)
        VMAT(NvF+I1,Lp) =-VMAT(MvF+I1,Lp)
C dPi/dI2s
        VMAT(NvF+I2,Lp) =-VMAT(MvF+I2,Lp)
D     IF (PRINT.GT.0)  THEN
D     WRITE(PRINT,2040) Lp, Nv,Mv, I1,I2, VXI(Lp)
D    @,  VMAT(MVF+I1,Lp),VMAT(MvF+I2,Lp),VMAT(NvF+I1,Lp),VMAT(NvF+I2,Lp)
D2040 FORMAT(/,' Pnt Lp',I3, 2X,2I2,I4,I3, 1P,E11.3,3X,4E11.3)
D     END IF
C Scan over ancestor vertices
        DO 40 Kv=NVERTX,Nv,-1
        IF (.NOT.VtxVtx(Kv,Nv)) THEN
          GO TO 40
        END IF
        KvF = VOFF(Kv)
C Loop over tracks
        DO 20 Nt=1,NTRACK
C   check vertex association
        IF (.NOT.TrkVtx(Nt,Kv)) GO TO 20
C offset for track Nt
          NtF = TOFF(Nt)
          IF (Np .EQ. 1) THEN
C  (R,Phi) pointing constraint,  dP1/dc, dP1/dphi, dP1/dcotg, dP1/dXsv, dP1/dYsv
            TEMP(1) = DDA(Nt,1)*DELV(2) - DDA(Nt,5)*DELV(1)
            TEMP(2) = DDA(Nt,2)*DELV(2) - DDA(Nt,6)*DELV(1)
            TEMP(3) = 0.0
            TEMP(4) = DDA(Nt,3)*DELV(2) - DDA(Nt,7)*DELV(1)
            TEMP(5) = DDA(Nt,4)*DELV(2) - DDA(Nt,8)*DELV(1)
          ELSE
C  (R,Z) pointing constraint,    dP2/dc, dP2/dphi, dP2/dcotg, dP2/dXsv, dP2/dYsv
            C = PAR0(1,Nt) + PARDIF(1,Nt)
            TEMP(1) = DDA(Nt,JMAT+1)*DELV(3)
            TEMP(1) = TEMP(1) + TrkP4(Nt,3)*DELV(IMAT)/C
            TEMP(2) = DDA(Nt,JMAT+2)*DELV(3)
            TEMP(3) =-DELV(IMAT)*TrkP4(Nt,5)
            TEMP(4) = DDA(Nt,JMAT+3)*DELV(3)
            TEMP(5) = DDA(Nt,JMAT+4)*DELV(3)
          END IF
          VMAT(NtF+1,Lp) = TEMP(1)
          VMAT(NtF+2,Lp) = TEMP(2)
          VMAT(NtF+3,Lp) = TEMP(3)
          VMAT(KvF+1,Lp) = TEMP(4) + VMAT(KvF+1,Lp)
          VMAT(KvF+2,Lp) = TEMP(5) + VMAT(KvF+2,Lp)
D     IF (PRINT.GT.0)  THEN
D     WRITE(PRINT,2041) NT,LP, TEMP
D2041 FORMAT(7X,I3,I5,22X,7E11.3)
D     END IF
C End track loop
   20   CONTINUE
C End ancestor vertex loop
   40   CONTINUE
C End pointing constraints loop
   50 CONTINUE
C===Return to Caller:
  100 CONTINUE
C Symmetrize the derivative matrix
      DO I=1,MATDIM-1
      DO J=I+1,MATDIM
         VMAT(J,I) = VMAT(I,J)
      END DO
      END DO
      RETURN
      END
C===============================================================================
      SUBROUTINE CCTVMMF (PRINT, NM, VXI)
C===============================================================================
C  Derivative contributions for the CTVM mass constraint fitting
C--Input parameter
C   Nm           The mass constraint index
C-------------------------------------------------------------------------------
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
      INTEGER PRINT
      INTEGER NM
      REAL*8  VXI(MAXDIM)
      REAL*8  SUM
      INTEGER I,J, Nt,Nv, LmF,NvF,NtF
      REAL    C
C-------------------------------------------------------------------------------
      LmF = MOFF + Nm
C  Difference in m**2 of tracks and constraint / 2
      SUM = SQRT(McnP4(1,Nm)**2 + McnP4(2,Nm)**2 + McnP4(3,Nm)**2)
      SUM = (McnP4(4,Nm) + SUM) * (McnP4(4,Nm) - SUM)
      SUM = DSQRT(SUM)
      VXI(LmF) = 0.5D0 * (CMASS(Nm)+SUM) * (CMASS(Nm)-SUM)
C  Loop over tracks contributing to this mass constraint
      DO 50 Nv=1,NVERTX
      NvF = VOFF(Nv)
      DO 40 NT=1,NTRACK
      IF (.NOT.TrkVtx(Nt,Nv)) GO TO 40
      IF (.NOT.TrkMcn(Nt,Nm)) GO TO 40
C Index into matrix equations for track NT
        NtF = TOFF(Nt)
C  dM**2/dc, dM**2/dphi, dM**2/dcot(theta)
        C   = PAR0(1,Nt) + PARDIF(1,Nt)
        SUM = -McnP4(4,Nm)*TrkP4(Nt,6)**2/TrkP4(Nt,4)
        SUM = (SUM + McnP4(3,Nm)*TrkP4(Nt,3)) / C
        SUM = SUM - McnP4(1,Nm)*DDA(Nt,1) - McnP4(2,Nm)*DDA(Nt,5)
        VMAT(NtF+1,LmF) = SUM
        SUM =-McnP4(1,Nm)*DDA(Nt,2) - McnP4(2,Nm)*DDA(Nt,6)
        VMAT(NtF+2,LmF) = SUM
        SUM = (McnP4(4,Nm)*TrkP4(Nt,3) / TrkP4(Nt,4)-McnP4(3,Nm))
        SUM = TrkP4(Nt,5) * SUM
        VMAT(NtF+3,LmF) = SUM
C  dM**2/dXsv, dM**2/dYsv
        VMAT(NvF+1,LmF) = VMAT(NvF+1,LmF)
     &              - McnP4(1,Nm)*DDA(Nt,3)-McnP4(2,Nm)*DDA(Nt,7)
        VMAT(NvF+2,LmF) = VMAT(NvF+2,LmF)
     &              - McnP4(1,Nm)*DDA(Nt,4)-McnP4(2,Nm)*DDA(Nt,8)
D     IF (PRINT.GT.0)  THEN
D     WRITE(PRINT,1045) NT, VXI(LmF)
D    @,     (VMAT(NvF+J,LmF),J=1,2),(VMAT(NtF+J,LmF),J=1,3)
D1045 FORMAT(/,' Mcn, Tk',I3, 1P,E11.3,2X,5E11.3)
D     END IF
   40 CONTINUE
C Continue loop on tracks
   50 CONTINUE
C-------------------------------------------------------------------------------
C===Return to Caller:
C Symmetrize the derivative matrix
      DO I=1,MATDIM-1
      DO J=I+1,MATDIM
         VMAT(J,I) = VMAT(I,J)
      END DO
      END DO
      RETURN
      END
C===============================================================================
      SUBROUTINE CCTVMPR(LSUNIT,DBGPRT,PARERR)
C===============================================================================
C     Author:  John Marriner, CDF, FNAL  (with a little help from his friends)
C     print the fit results of CTVMFT
C     Input parameters
C     LSUNIT = fortran logical unit for write
C     DBGPRT = print level
C     PARERR   Original (input) track parameters
C     Preconditions
C     CTVMFT must have been called and returned no error
C     CTVMPR assumes fit information is stored in common defined by CTVMFT.INC
C-------------------------------------------------------------------------------
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
      INTEGER DBGPRT, LSUNIT
      REAL    PARERR(5,MaxTrk)
      INTEGER I,J,K,L,M,N, IFAIL, IM, IN, INDEX(5)
      INTEGER Nt,Nv,Mv, Nc,Np,Nm, NtF,NvF, NPAR, NVPAR, NPC
      INTEGER NtrVtx(0:Maxvtx), TMcn(MaxTrk)
      LOGICAL PRIMARY
      REAL    CHIP, VALUE,SIGMA, U,V,W, ERRRAT(5), PULL(3)
      REAL    DXYZ(3), Dr,Dz, Dl(3)
      REAL    ERR(MAXDIM), WORK(MAXDIM), CORR(MAXDIM,MAXDIM)
      REAL*8  P4(4), EPAR(5,5)
      INTEGER MRUN,MNEV
      INTEGER MTRACK(MaxTrk)
      LOGICAL FIRST
      CHARACTER*9 VN, VP, VS, VB
      SAVE    MRUN,MNEV, VP,VS,VB
      DATA    VP / 'primary  ' /, VS / 'secondary' /, VB / '         ' /
      DATA    MRUN /-1/, MNEV /-1/
C-------------------------------------------------------------------------------
C            Write header, chi-squared, # deg of freedom
      IF (RunNum.NE.MRUN .OR. TrgNum.NE.MNEV) THEN
        WRITE(LSUNIT,1000)
 1000   FORMAT('1',123('=') )
        MRUN = RunNum
        MNEV = TrgNum
      END IF
      DO I=1,MaxTrk
        MTRACK(I) = I
      END DO
      Nc = 0
      Np = 0
      PRIMARY =.FALSE.
      NtrVtx(0) = 0
      DO Nv=1,NVERTX
        NtrVtx(Nv) = 0
        Do Nt=1,NTRACK
          IF (TrkVtx(Nt,Nv)) NtrVtx(Nv) = NtrVtx(Nv)+1
        END DO
        IF (VtxPnt(Nv,1).EQ.0) PRIMARY =.TRUE.
        Nc  = Nc + Cvtx(Nv)
        IF (VtxPnt(Nv,1).GE.0) THEN
          NPC = VtxPnt(Nv,2)
          IF (NPC.EQ.3) NPC = 2
          IF (NPC.EQ.4) NPC = 2
          Np = Np+NPC
        END IF
      END DO
C           make the fit errors and the correlation matrix
      NPAR = TOFF(NTRACK)+3
C Loop over fitted parameters
      DO I=1,NPAR
C Get error = sqrt(diagonal element)
        ERR(I) = DSQRT(VMAT(I,I))
      END DO
C Correlation matrix
      DO I=1,NPAR
      DO J=1,NPAR
        CORR(I,J) = VMAT(I,J)/(ERR(I)*ERR(J))
      END DO
      END DO
      I = RunNum
      J = TrgNum
      CHIP = AMIN1(CHISQR(0),999.0)
      CALL CMCALC(NTRACK,MTRACK, Value,SIGMA, P4)
      VALUE = AMIN1(VALUE,999.999)
      SIGMA = AMIN1(SIGMA, 99.999)
      WRITE(LSUNIT,1001) I,J,CHIP,NDOF,ITER,NTSCUT
     &,                  NTRACK,NVERTX,Np,Nc,NMASSC,Value,SIGMA
 1001 FORMAT
     &(1X,123('-'),  /,' Event ',I5,'.',I6,';  Ct_VM fit results'
     &, 5X,'Chi square =',F7.2,'  Ndof',I3,'  (Iter',I3,', NtCut',I3,')'
     &, 20X,'Mass  Sigma'
     &, /, 18X,I2,' Tracks,',I2,' Vertices;    Constraints (Pointing',I3
     &, ',  Conversion',I2,',  Mass',I2,')',16X,F8.3,F7.3
     & )
C            write vertex fit Chi-Square results
      WRITE(LSUNIT,1010)
 1010 FORMAT(' Vertex coordinate fit results' ,/
     &,   11X,'Vtx  Ntr',5X,'Xv',7X,'Yv',8X,'Zv',7X,'Chisq'
     &,    3X,'Np  Nc    DLr       Dlz'
     &,   12X,'SumPx   SumPy   SumPz   Sum_E')
      L = 1
      IF (PRIMARY) L=0
      DO Nv=L,NVERTX
        IF (Nv.GT.0) THEN
          NvF = VOFF(Nv)
        END IF
        CHIP = AMIN1(CHIV(Nv),999.0)
        Dr = 0
        Dz = 0
        IF (Nv.EQ.0) THEN
          NvF = 0
          Np =-1
          Nc =-1
          WRITE(LSUNIT,1015) VP,Nv,NtrVtx(Nv),(XYZVRT(I,Nv),I=1,3),CHIP
          WRITE(LSUNIT,1016) (ERR(NvF+I),I=1,3)
        ELSE
          VN = VS
          IF (Nv.GT.1) VN=VB
          Mv = VtxPnt(Nv,1)
          Nc = Cvtx(Nv)
          Np = VtxPnt(Nv,2)
          IF (Np.LT.0) Np = 0
          IF (NP.EQ.3) NP = 2
          IF (NP.EQ.4) NP = 2
          IF (Np.LT.1) PCON(Nv,1)=0.0
          IF (Np.LT.2) PCON(Nv,2)=0.0
          NvF = VOFF(Nv)
          IF (Np.LE.0) THEN
            WRITE(LSUNIT,1015) VN,Nv,NtrVtx(Nv),(XYZVRT(I,Nv),I=1,3)
     &,            CHIP,Np,Nc,      (VtxP4(I,Nv),I=1,4)
            WRITE(LSUNIT,1016) (ERR(NvF+I),I=1,3)
          ELSE
            CALL CDCALC(Nv,Mv,DXYZ,Dr,Dz,Dl)
            WRITE(LSUNIT,1017) VN,Nv,Mv,NtrVtx(Nv),(XYZVRT(I,Nv),I=1,3)
     &,            CHIP,Np,Nc,Dr,Dz,(VtxP4(I,Nv),I=1,4)
            WRITE(LSUNIT,1018) (ERR(NvF+I),I=1,3),Dl
          END IF
        END IF
 1015 FORMAT(1X,A9,I3,I5,       F11.4,F9.4,F10.4,F9.2
     &,             I5,I4,20X   ,F14.3,3F8.3)
 1016 FORMAT(20X,2F9.4,F10.4      , 58X,3F6.3)
 1017 FORMAT(1X,A9,I3,',',I1,I3,2X, F9.4,F9.4,F10.4,F9.2
     &,     I5,I4,F9.4,F10.4,5X,F10.3,3F8.3)
 1018 FORMAT(20X,2F9.4,F10.4,15X,2F7.4,F8.4, 20X,3F6.3)
      END DO
      WRITE(LSUNIT,*)
C            Write information on mass constraints
      IF (NMASSC.GT.0) THEN
        DO Nm=1,NMASSC
        L = 0
        DO Nt=1,NTRACK
          IF (TrkMcn(Nt,Nm)) THEN
            L = L+1
            TMcn(L) = Nt
          END IF
        END DO
        IF (Nm.EQ.1) THEN
          WRITE(LSUNIT,1020) Nm, CMASS(Nm),FMCDIF(Nm),(TMcn(I),I=1,L)
 1020 FORMAT(' Mass Constraints; Mcn',I2,' Mass ',F7.4,' Gev/c**2'
     &,      '  Residual ' ,1P,E9.1,' Mev/c**2   Tracks',10I3)
        ELSE
          WRITE(LSUNIT,1021) Nm, CMASS(Nm),FMCDIF(Nm),(TMcn(I),I=1,L)
 1021 FORMAT(I24, F13.4,17X,1P,E12.1 ,18X               ,10I3)
        END IF
        END DO
        WRITE(LSUNIT,*)
      END IF
C            Now, print the track fit results
      WRITE(LSUNIT,1050)
 1050 FORMAT(' Track parameter fit results -'
     &,/,5X,'Vtx  Bank     Mass     Crv',10X,'Phi',7X,'Ctg',8X,'D0'
     &,  7X,'Z0',8X,'Chisq',22X,'ErrFit / ErrMst')
      DO 50 Nv=1,Nvertx
      WRITE(LSUNIT,*)
      DO 49 Nt=1,Ntrack
      IF (.NOT.TrkVtx(Nt,Nv)) GO TO 49
        NtF = TOFF(Nt)
        VALUE = AMIN1(CHIT(Nt),999.0)
        DO I=1,3
          ERRRAT(I) = ERR(NtF+I) / PARERR(I,Nt)
        END DO
        WRITE(LSUNIT,1055) Nt,Nv, TKBANK(Nt),LIST(Nt),TMASS(Nt)
     &,                        (PAR(J,Nt),J=1,5), VALUE
        WRITE(LSUNIT,1056) (ERR(NtF+J),J=1,3), (ERRRAT(J),J=1,3)
 1055  FORMAT(I3,')',I3,3X,A4,I3,F8.4,1P,E12.3,0P,2F10.5,2F10.4,F9.2)
 1056  FORMAT(25X,1P,E12.3,0P,2F10.5, 49X, 3F6.3)
      IF (DBGPRT.LT.987654) GO TO 49
C            Loop over weight matrix
      DO I = 1, 5
      DO J = 1, 5
        EPAR(I,J) = G(I,J,NT)
      END DO
      END DO
C            Invert weight matrix to get error matrix
      CALL DINV(5,EPAR,5,WORK,IFAIL)
C            Get index into fit error matrix for track NT
C            Loop over fitted parameters (first 3 in a track)
      DO I = 1, 3
C            Get sigma = denominator for pull
      SIGMA = EPAR(I,I) - VMAT(NtF+I,NtF+I)
C            Allow modest underflow (rounding errors?)
      IF (ABS(SIGMA) .LT. 0.001D0*EPAR(I,I)) SIGMA = 0.001D0*EPAR(I,I)
C            Test for negative argument
      IF (SIGMA .LT. 0.0) GO TO 45
C            Pull = (Fitted value - CTC fit)/sigma
      PULL(I) = PARDIF(I,NT)/SQRT(SIGMA)
C            Continue loop on fitted parameters
      END DO
C            Write pulls, contribution to chi-squared for this track
      WRITE (LSUNIT,1062) (PULL(I),I=1,3)
 1062 FORMAT(6X,'  Pulls',0P,3F9.4)
      GO TO 49
   45 CONTINUE
      WRITE (LSUNIT,21)
   21 FORMAT (3X,'  Pulls. Error in calculating pulls.',
     &     '  Fit is suspect.')
C            Continue loop on tracks
   49 CONTINUE
   50 CONTINUE
C Loop over fitted parameters
   60 CONTINUE
      IF (DBGPRT.LE.99) GO TO 900
C Write header for error matrix print
      WRITE (LSUNIT,1060)
 1060 FORMAT(/,' Correlation matrix')
      DO I=1,NPAR
        WRITE(LSUNIT,1065) (CORR(I,J),J=1,I)
 1065 FORMAT( 8(1X,3F5.2) )
        IF (MOD(I,3).EQ.0) WRITE(LSUNIT,*)
      END DO
C            Done. Exit
  900 CONTINUE
      RETURN
      END
C===============================================================================
      SUBROUTINE CCTVMsv (WHICH, WHERE)
C===============================================================================
C  Silly, trivial little routine to save, interchange, or restore the results
C  of a CTVMFT fit.
C  WHICH < 0   Copy the CTVMFT fit results to local storage for preservation
C        = 0   Interchange the stored fit result with that in the CTVMFT common
C        > 0   Overwrite the CTVMFT fit result with the stored result
C  WHERE       An index specifying where to store or where to extract fit
C              information.
C-------------------------------------------------------------------------------
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
      INTEGER   IJKLMN
      PARAMETER(IJKLMN = MAXDIM*(MAXDIM+1))
      INTEGER   WHICH, WHERE
      INTEGER   I
      INTEGER   SAVRSL(UDIM,4), SLUSH(UDIM)
      REAL*8    UMAT (MAXDIM*(MAXDIM+1),4)
      REAL*8    BLUSH(MAXDIM*(MAXDIM+1))
      REAL      SCREAM, MINUS, ZERO
      SAVE      SAVRSL, UMAT, MINUS,ZERO
      DATA      MINUS /-1.0/, ZERO /0.0/
C============ save or restore existing constrained fit results ===============
      IF (WHERE.LE.0 .OR. WHERE.GT.4) THEN
        PRINT 1000, WHERE
 1000   FORMAT (' CTVMsave; nonsense WHERE!' )
        SCREAM = SQRT(MINUS) / ZERO
        STOP
      END IF
      IF (WHICH) 10,20,30
C save the existing fit from CTVMFq,r
   10 CONTINUE
      DO I=1,UDIM
        SAVRSL(I,WHERE) = UVWXYZ(I)
      END DO
      DO I=1,IJKLMN
        UMAT(I,WHERE) = VMAT(I,1)
      END DO
      GO TO 100
C flip the stored fit with that in CTVMFq,r
   20 CONTINUE
      DO I=1,UDIM
        SLUSH (I) = SAVRSL(I,WHERE)
        SAVRSL(I,WHERE) = UVWXYZ(I)
        UVWXYZ(I) = SLUSH (I)
      END DO
      DO I=1,IJKLMN
        BLUSH(I)  = UMAT(I,WHERE)
        UMAT (I,WHERE)  = VMAT(I,1)
        VMAT(I,1) = BLUSH(I)
      END DO
      GO TO 100
C restore the CTVMFq,r contents
   30 CONTINUE
      DO I=1,UDIM
        UVWXYZ(I) = SAVRSL(I,WHERE)
      END DO
      DO I=1,IJKLMN
        VMAT(I,1) = UMAT(I,WHERE)
      END DO
  100 CONTINUE
      RETURN
      END
C======================================================================= =======
      SUBROUTINE CMCALC (NTRKQ,MTRACK, FMASS,DMASS, P4)
C======================================================================= =======
C
C=======Description
C
C  Input parameters:
C  Requires the common in the CTVMFT include file to have relevant information
C  stored from a CTVMFT entrance.
C  NTRKQ = Number of tracks comprising the invariant mass.
C  MTRACK(i) selects the track specified in LIST(i) in the CTVMFT call
C  Output parameters:
C  FMASS = invariant mass of track combination
C  DMASS = computed error on the mass (will be returned as negative if
C           DMASS**2 was computed to be negative. This should not happen
C           if the matrix VMAT is valid)
C=======Author
C  John Marriner
C  6 April 1992
C----------------------------------------------------------------------- -------
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
      INTEGER      NTRKQ
      INTEGER      MTRACK(MAXTRK)
      REAL         FMASS, DMASS
      INTEGER      I,J, IERR, Nv,NvF, Nt,NtF, NDIM, Nvx,NVRTX(MaxVtx)
      REAL         DMASSQ
      REAL*8       SUM,QMASS, P4(4), DM(MAXDIM)
      CHARACTER*80 STRING
      CHARACTER*6  NAME
      DATA         NAME /'MCALC '/
C===============================================================================
      CALL VZERO(P4,2*4)
      CALL VZERO(DM,2*MAXDIM)
      FMASS = 0.0
      DMASS = 0.0
C complain
      IF (NTRKQ.GT.MAXTRK) THEN
        IERR = 91
        WRITE(STRING,901) NAME,IERR,NTRKQ
  901   FORMAT(A6,I3,'; ',I2,' is too many tracks')
        CALL CERROR('MCALC ',IERR,STRING)
        FMASS =-1.0
        GO TO 100
      END IF
      DO J=1,NTRKQ
        Nt = MTRACK(J)
C energy
        TrkP4(NT,4) = SQRT(TrkP4(NT,6)**2+TMASS(NT)**2)
        DO I=1,4
          P4(I) = P4(I) + TrkP4(Nt,I)
        END DO
      END DO
      NDIM = 0
      Nvx  = 0
      CALL VZERO(NVRTX,MAXVTX)
      DMASSQ = 0.0
C Loop over tracks in the fit, find dM
      DO 10 J=1,NTRKQ
        Nt = MTRACK(J)
        NtF = TOFF(Nt)
C find the vertex for this track
        Nv = 1
        DMASSQ = DMASSQ + TMASS(Nt)
        DO WHILE (.NOT.TrkVtx(Nt,Nv))
          Nv = Nv+1
        END DO
C count number of vertices involved
        IF (NVRTX(Nv).EQ.0) THEN
          Nvx = Nvx+1
          NVRTX(Nv) = Nv
        END IF
        NvF = VOFF(Nv)
C maximum matrix dimension in VMAT
        NDIM = MAX0(NDIM,NtF+3)
C first,dM**2/dXsv and  dM**2/dYsv, treat the displacement to this vertex
        DM(NvF+1) = DM(NvF+1)-2.0*(P4(1)*DDA(Nt,3)+P4(2)*DDA(Nt,7))
        DM(NvF+2) = DM(NvF+2)-2.0*(P4(1)*DDA(Nt,4)+P4(2)*DDA(Nt,8))
C then, dM**2/dc, dM**2/dphi, dM**2/dcot(theta), helix parameters for this track
        SUM = -P4(4)*TrkP4(Nt,6)**2/TrkP4(Nt,4) + P4(3)*TrkP4(Nt,3)
        SUM = SUM/PAR(1,Nt) - (P4(1)*DDA(Nt,1) + P4(2)*DDA(Nt,5))
        DM(NtF+1) = 2.0*SUM
        DM(NtF+2) =-2.0*(P4(1)*DDA(Nt,2) + P4(2)*DDA(Nt,6))
        SUM = (P4(4)*TrkP4(Nt,3)/TrkP4(Nt,4) - P4(3)) * TrkP4(Nt,5)
        DM(NtF+3) = 2.0*SUM
   10 CONTINUE
C               mass, mass**2
      QMASS = DSQRT(P4(1)**2+P4(2)**2+P4(3)**2)
      IF (P4(4).LE.QMASS) THEN
        IF (Nvx.EQ.1 .AND. Cvtx(Nv).GT.0) THEN
          QMASS = DMASSQ
        ELSE
C ? crunch
C 10/23/02  Change next line by adding minus signs so that code
C does not bomb is energy is less than momentum due to roundoff error.
C Craig Blocker
          QMASS = -DSQRT(-(P4(4)+QMASS)*(P4(4)-QMASS))
        END IF
      ELSE
        QMASS = DSQRT((P4(4)+QMASS)*(P4(4)-QMASS))
      END IF
      FMASS = QMASS
C Set error in mass squared = 0 for summing
      DMASSQ = 0.0
      DO I=1,NDIM
      DO J=1,NDIM
        DMASSQ = DMASSQ + DM(I)*VMAT(J,I)*DM(J)
      END DO
      END DO
C check for, and flag, negative error**2
      IF (DMASSQ.GE.0.) THEN
         DMASS = SQRT(DMASSQ)
      ELSE
        DMASS = -SQRT(-DMASSQ)
      ENDIF
      IF (QMASS.EQ.0.0) THEN
        FMASS = -0.001
      ELSE
        DMASS = DMASS/(2.*QMASS)
      END IF
  100 CONTINUE
      RETURN
      END
C======================================================================= =======
      SUBROUTINE CEMASSQ(NTRKQ,QMASS, PXYZE,IERR)
C======================================================================= =======
C  EMASSQ   calculates M**2, sigma(M**2), via a call th MCALC
C           This is furnished as a convenience for people who have calls to this
C           obsolescent routine. Note the more restrictive set of track pointers
C           as compared with MCALC.
C----------------------------------------------------------------------- -------
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
      INTEGER   NTRKQ, IERR
      REAL      QMASS(MAXTRK), PXYZE(6)
      INTEGER   I,N, MTRACK(MAXTRK)
      REAL      FMASS, DMASS
      REAL*8    P4(4)
C----------------------------------------------------------------------- -------
      IERR = 0
      CALL VZERO(P4,   2*4)
      CALL VZERO(PXYZE,6)
      N = MIN0(MaxTrk,NTRKQ)
      DO I=1,N
        MTRACK(I) = I
        TMASS(I) = QMASS(I)
      END DO
      CALL CMCALC (NTRKQ,MTRACK,FMASS,DMASS, P4)
      IF (FMASS.LT.0.0) THEN
        IERR = FMASS
        GO TO 900
      END IF
      DO I=1,4
        PXYZE(I) = P4(I)
      END DO
      PXYZE(5) = FMASS**2
      PXYZE(6) = 2.0 * FMASS * DMASS
  900 CONTINUE
      RETURN
      END
C===============================================================================
      SUBROUTINE CDCALC (Nv,Mv,DXYZ,Dr,Dz,Dl)
C===============================================================================
C      Auxialliary subroutine for use with CTVMFT.
C      Calculates the distance between the two vertices Nv, Mv and its error.
C      Mostly useful only if Nv "points at" Mv.
C      Input parameters:  Nv, Mv  The vertex indices.
C          The rest of the Input data comes from the common block /CTVMFT/
C          which must contain a valid fit
C      Output parameters:
C      DXYZ = (unrotated) x, y, and z components of vertex displacement
C      Dr   = (signed) radial vertex separation (dot product of Dxy * VtxPxy)
C      Dz   = (signed) vertex separation (dot product of Dz * VtxPz)
C      Dl   = uncertainty of the three-vector (Dr,Dphi,Dz) (note Dphi~0)
C-------------------------------------------------------------------------------
        IMPLICIT NONE
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"
      INTEGER      Nv,Mv
      REAL         Dxyz(3), Dr,Dz, Dl(3), FL(3), U
      INTEGER      I,J,K, IERR, NvF,MvF
      REAL*8       S, ROT(3,3), E(3,3),F(3,3)
      CHARACTER*80 STRING
C-------------------------------------------------------------------------------
C  Setup vertex offset pointers, check for primary
      IF (Nv.LT. 1 .OR. Nv.GT.NVERTX) THEN
        IERR = 99
        WRITE(STRING,900) Nv
  900   FORMAT(' Vertex Nv',I3,' is invalid')
        CALL CERROR ('DCALC ',IERR,STRING)
ccc        STOP
      ELSE IF (Mv.LT.0 .OR. Mv.GE.Nv) THEN
        IERR = 98
        WRITE(STRING,901) Mv
  901   FORMAT(' Vertex Mv',I3,' is invalid')
        CALL CERROR ('DCALC ',IERR,STRING)
ccc        STOP
      END IF
      DO I=1,9
        ROT(I,1) = 0.0D0
      END DO
C            Pointers to information in VMAT for the two vertices
      NvF = VOFF(Nv)
      IF (Mv.GT.0) THEN
        MvF = VOFF(Mv)
      ELSE
        MvF =-1
        DO I=1,NVERTX
C the primary participated in this fit
          IF (VtxPnt(I,1).EQ.0) THEN
            MvF = 0
          END IF
        END DO
      END IF
C covariance matrix on vertex separation
      DO I=1,3
      DO J=I,3
        E(I,J) = VMAT(NvF+I,NvF+J)
        IF (MvF.GE.0) THEN
          E(I,J) = E(I,J) + VMAT(MvF+I,MvF+J)
          E(I,J) = E(I,J) - VMAT(NvF+J,MvF+I) * 2.0D0
        ELSE
          E(I,J) = E(I,J) + EXYZPV(I,J)
        END IF
        E(J,I) = E(I,J)
      END DO
      END DO
C vertex separation vector
      IF (MvF.GE.0) THEN
        DO I=1,3
        DXYZ(I) = XYZVRT(I,Nv) - XYZVRT(I,Mv)
        END DO
      ELSE
        DO I=1,3
        DXYZ(I) = XYZVRT(I,Nv) - XYZPV0(I)
        END DO
      END IF
      S = SQRT(DXYZ(1)**2 + DXYZ(2)**2)
c If the two vertices have exactly the same x and y, then S is
c zero and we need to protect against dividing by S.  C. Blocker 11/20/03
c      if(s.le.0.) then
      IF (ABS(S).LE.1E-42) THEN
        dr = 0.
        dz = 0.
        dl(1) = -999.
        dl(2) = -999.
        dl(3) = -999.
        return
      endif
      ROT(1,1) = DXYZ(1)/S
      ROT(1,2) = DXYZ(2)/S
      ROT(2,1) =-ROT(1,2)
      ROT(2,2) = ROT(1,1)
      ROT(3,3) = 1.0D0
      do i=1,2
        s = 0.0D0
        do j=1,2
          s = s + rot(i,j)*DXYZ(J)
        end do
        fl(i) = s
      end do
C pre and post multiply xy part
      DO I=1,2
      DO J=1,2
        S = 0.0D0
        DO K=1,2
          S = S + ROT(I,K) * E(K,J)
        END DO
        F(I,J) = S
      END DO
      END DO
      DO I=1,2
      DO J=1,2
        S = 0.0D0
        DO K=1,2
          S = S + F(I,K) * ROT(J,K)
        END DO
        E(I,J) = S
      END DO
      END DO
      U  = DXYZ(1)*VtxP4(1,Nv) + DXYZ(2)*VtxP4(2,Nv)
C (signed) xy vertex separation
      Dr = SIGN(SQRT(DXYZ(1)**2+DXYZ(2)**2),U)
      U  = DXYZ(3)*VtxP4(3,Nv)
C (signed)  z vertex separation
      Dz = SIGN(ABS(DXYZ(3)),U)
      DO I=1,3
        IF (E(I,I).GT.0.0) THEN
          Dl(I) = SQRT(E(I,I))
        ELSE
          Dl(I) = 0.0
        END IF
      END DO
      RETURN
      END
