C --------- Include file to define dimensions --------------------------
C
C The CTVMFT_ variables are picked up by defines in the preprocessing
C
C Maximum matrix dimension (MAXDIM);
C             = 5 parameters per (vtx + primary) (x,y,z + Phi,Ctg)
C                 (3 positions + 2 pointing variables)
C             + 3 parameters per track (Crv, Phi, Ctg)
C             + number of mass constraints

      INTEGER    MAXVTX,MAXTRK,MAXMCN,MAXDIM
      PARAMETER (MAXVTX = CTVMFT_MAXVTX)  ! Maximum number of vertices
      PARAMETER (MAXTRK = CTVMFT_MAXTRK)  ! Maximum number of tracks
      PARAMETER (MAXMCN = CTVMFT_MAXMCN)  ! Maximum number of mass constraints
      PARAMETER (MAXDIM = 5*(MAXVTX+1) + 3*MAXTRK + MAXMCN)

C Parameters important to the procedure
      INTEGER    MAXITR
      PARAMETER (MAXITR = CTVMFT_MAXITR)  ! Maximum number of iterations
C
      INTEGER    UDIM
      Parameter (UDIM   = MAXTRK*(MAXVTX+MAXMCN+60)
     &                  + MAXVTX*(MAXVTX+18)
     &                  + MAXMCN*7
     &                  + MAXITR
     &                  + 31)
