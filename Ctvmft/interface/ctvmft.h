C========================= Include ctvmft.inc ==========================
C============ include file for the constrained fit module CTVMFT =======

C ----------------------------------------------------------------------
C Single precision stuff

C Counters of the used numbers of Vtx ..

      INTEGER     NVERTX                ! Number of vertices
      INTEGER     NMASSC                ! Number of mass constraints
      INTEGER     NTRACK                ! Number of tracks, this fit

C Vertex track associations

      LOGICAL     TRKVTX(MAXTRK,MAXVTX) ! (Track,Vertex) assoc. table
      INTEGER     VTXPNT(MAXVTX,2)      ! Vertex association information

C Conversion constraints info

      INTEGER     CVTX  (MAXVTX)        ! Conversion constraint specs
      LOGICAL     TRKMCN(MAXTRK,MAXMCN) ! (Track,Mass_Constraint) table
      REAL        CMASS (MAXMCN)        ! Constraint mass

C  Required Input Track Specifications

      CHARACTER*4 TKBANK(MAXTRK)        ! Track bank type (e.g. 'TRKS')
      INTEGER     LIST  (MAXTRK)        ! Track bank number
      REAL        TMASS (MAXTRK)        ! Track mass assignment

C  Fit input, collected using the input specifications

      REAL        PAR0  (5,MAXTRK)      ! Original trk parameter values
      REAL        G     (5,5,MAXTRK)    ! Fitted trk parameter covar.

C  Optional Input, Required Data for Pointing Constraints

      REAL        XYZPV0(3)             ! Primary vertex: XYZ
      REAL        EXYZPV(3,3)           ! Primary vertex: covariance
      REAL        XZSLOPE,YZSLOPE       ! beamline slopes vs z direction
      INTEGER     RUNNUM,TRGNUM         ! Run and Event Number
      INTEGER     ITER                  ! Number of iterative steps
      INTEGER     NTSCUT                ! Total number of cut steps

C  Output; Results of the Fit

      INTEGER     MATDIM                ! Dimension of matrix, this fit
      INTEGER     NDOF                  ! Number of degrees of freedom
      LOGICAL     VTXVTX(MAXVTX,MAXVTX) ! Vertex geneology
      REAL        CHISQR(0:MAXITR)      ! Overall fit Chi Square result
      REAL        CHIT  (MAXTRK)        ! Track contr. to Chi Square
      REAL        CHIV  (0:MAXVTX)      ! Vertex fit Chi Squares
      REAL        CHIM  (MAXMCN)        ! Mass constraint Chi Squares
      REAL        PAR   (5,MAXTRK)      ! Fitted trk helix parameters
      REAL        PARDIF(5,MAXTRK)      ! Parameter diff (fit-input)
      REAL        FMCDIF(MAXMCN)        ! Mass constraint residuals
      REAL        PCON  (MAXVTX,2)      ! Pointing constraint tests
      REAL        SANG  (MAXVTX,2)      ! Pointing constraint tests
      REAL        XYZVRT(3,0:MAXVTX)    ! Primary, Second vertices
      REAL        TRKP4 (MAXTRK,6)      ! Track Px,Py,Pz,E,Pt P
      REAL        VTXP4 (4,MAXVTX)      ! Vertex 4-momentum sum
      REAL        MCNP4 (4,MAXMCN)      ! 4-mom sum for mass constr.
      REAL        DDA   (MAXTRK,8)      ! d(Px,Py)/d(Crv,Phi,Xs,Ys)
      REAL        DXYZPV(3)             ! Primary displacement in fit
      INTEGER     VOFF  (MAXVTX)        ! Vertex      offset pointers
      INTEGER     POFF  (MAXVTX)        ! Point Phi   offset pointers
      INTEGER     COFF  (MAXVTX)        ! Point Cotan offset pointers
      INTEGER     TOFF  (MAXTRK)        ! Track       offset pointers
      INTEGER     MOFF                  ! Mass        offset pointers

C Memory

      INTEGER     TKERR (MAXTRK)        ! Trk error flag, previous fit
      INTEGER     IJKERR(3)             ! Error code reporting
      REAL        XYZPV (3)
      INTEGER     UVWXYZ(UDIM)
      EQUIVALENCE(XYZPV,XYZVRT(1,0))    ! Danger: equiv. statements
      EQUIVALENCE(UVWXYZ,ITER)

C  CTC first approximation and track-vertex step parameters

      REAL    DRMAX,RVMAX,DZMAX,TRNMAX,DSMIN

C Useful conversion parameter

      REAL    PSCALE                    ! Mom.scale conversion
                                        ! [cm^-1 -> Gev/c]

      COMMON /CTVMFR/
     >        RUNNUM,TRGNUM,ITER,NTSCUT
     >       ,NVERTX,NMASSC,NTRACK
     >       ,TRKVTX,TRKMCN,VTXPNT,CMASS,CVTX
     >       ,VTXVTX
     >       ,TKBANK,LIST,TMASS
     >       ,MATDIM
     >       ,TKERR
     >       ,NDOF,CHISQR,CHIT,CHIV,CHIM
     >       ,XYZPV0,EXYZPV
     >       ,XZSLOPE,YZSLOPE                       
     >       ,XYZVRT,DXYZPV
     >       ,PAR,G
     >       ,TRKP4,VTXP4,MCNP4,DDA
     >       ,VOFF,TOFF,POFF,COFF,MOFF
     >       ,PAR0,PARDIF
     >       ,FMCDIF,PCON,SANG
     >       ,DRMAX,RVMAX,DZMAX,TRNMAX,DSMIN
     >       ,IJKERR
     >       ,PSCALE

C ----------------------------------------------------------------------
C Double precision stuff

C Covariance matrix, fit parameters
      REAL*8      VMAT(MAXDIM,MAXDIM+1)
      COMMON /CTVMFD/
     >        VMAT
