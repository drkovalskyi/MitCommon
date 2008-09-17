C --------- Include file to define the CTVMFT communciations -----------
C  Communication back to CTVMFT for odds and ends, especially PARERR ---

      INTEGER FLIP, DBGPRT
      REAL    XSVI(2),YSVI(2),DZSVI(2),SS(2,2),ZZ(2,2),PARERR(5,MAXTRK)
      COMMON /CTVMUV/ XSVI,YSVI,DZSVI,SS,ZZ, FLIP, PARERR, DBGPRT
