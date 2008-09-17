c $Id:$

* Inital version of CMS port of CTVMFT (CDF vertex fitter).

      SUBROUTINE VZERO (A,N)
C
C CERN PROGLIB# F121    VZERO           .VERSION KERNFOR  4.40  940929
C ORIG. 01/07/71, modif. 24/05/87 to set integer zero
C                 modif. 25/05/94 to depend on QINTZERO
C
      INTEGER   A(*)
C
      IF (N.LE.0)  RETURN
      DO 9 I= 1,N
    9 A(I)= 0
C
      RETURN
      END
