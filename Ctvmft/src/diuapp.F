c $Id:$

C VAX/DEC CMS REPLACEMENT HISTORY, Element DIUAPP.CDF
C *1    25-MAR-1992 13:06:59 SANSONI "calculates crude secondary vtx position"
C VAX/DEC CMS REPLACEMENT HISTORY, Element DIUAPP.CDF
	SUBROUTINE DIUAPP(PAR0,XSV,YSV,ZSV,DLZ,IERR)
C	============================================
C
C	Initial approximation to the secondary vertex position
C	This subroutine intersects the two circles that are
C	the projections of the helix trajectories of the
C	tracks with the given input parameters.  If the circles intersect,
C	the program choses the one with the smaller z difference
C	between the 2 helices.  If the circles do not intersect
C	it simply exits.
C
C	Author:  John Marriner
C	Original creation: Feb 22, 1991
C
C===============================================================================

C $$IMPLICIT

	INTEGER  IERR,I
	DOUBLE PRECISION R1, X1, Y1, R2, X2, Y2, DX, DY, S, Z1, Z2
	DOUBLE PRECISION A, B, C, RAD, SRAD, XA, XB
	DOUBLE PRECISION XSVI(2), YSVI(2), ZSVI(2), DZSVI(2)
	REAL CKA, CKB, DZMAX, PI, TANG, XHAT, YHAT, PAR0(5,2)
	REAL XSV, YSV, ZSV, DLZ

	IERR = -1
	XSV  =0
	YSV  =0
	ZSV  =0
	DLZ  =0

	DZMAX = 1.E+4	! Maximum acceptable z difference for a solution
	PI = 4.0*ATAN(1.0)	    ! Pi(the constant)
	R1 = 1.0/(2.0*PAR0(2,1))    ! Radius of first track (with sign)
	X1 = -(R1+PAR0(4,1))*SIN(PAR0(5,1)) ! Center (x) of the first track
	Y1 =  (R1+PAR0(4,1))*COS(PAR0(5,1)) ! Center (y) of the first track
	R2 = 1.0/(2.0*PAR0(2,2))    ! Radius of the second track
	X2 = -(R2+PAR0(4,2))*SIN(PAR0(5,2)) ! Center (x) of the second track
	Y2 =  (R2+PAR0(4,2))*COS(PAR0(5,2)) ! Center (y) of the second track
	DX = X2 - X1		    ! Center to center difference in x
	DY = Y2 - Y1		    ! Center to center difference in y
C
C  Vertex y is the solution of A x y**2 + B x y + C = 0
C
	A = 4.0*(DX**2+DY**2)
	B = 4.0*DY*(R2*R2-R1*R1-DX*DX-DY*DY)
	C = (R2*R2-R1*R1)**2 + (DX*DX+DY*DY)**2 -
     -	    2.0*DX*DX*(R2*R2+R1*R1) - 2.0*DY*DY*(R2*R2-R1*R1)
	RAD = B*B - 4.0*A*C	! Argument of square root

C  RAD<0 means circles don't intersect.It occurs for parallel trajectories only
C  in which case the whole strategy of the approximation may not make sense.

	IF (RAD .LT. 0.) RETURN
	SRAD = DSQRT(RAD)

C  Get y at circle crossings
	YSVI(1) = (-B+SRAD)/(2.0*A)
	YSVI(2) = (-B-SRAD)/(2.0*A)

C  Loop over 2 potential vertices
	DO 100 I=1,2

C  The solution for X involves a SQRT with a sign ambiguity.
C  The following code resolves the ambiguity

C	   XA = x if plus sign is correct
	   XA =  DSQRT(R1*R1-YSVI(I)*YSVI(I)) + X1
C	   XB = x if minus sign is correct
	   XB = -DSQRT(R1*R1-YSVI(I)*YSVI(I)) + X1
C	   Translate Y to standard coordinates
	   YSVI(I) = YSVI(I) + Y1
C	   CKA = 0, if XA is correct
	   CKA = R2**2 - (XA-X2)**2 - (YSVI(I)-Y2)**2
C	   CKB = 0, if XB is correct
	   CKB = R2**2 - (XB-X2)**2 - (YSVI(I)-Y2)**2
C	   Choose vertex x according to whether CKA or CKB is smaller
	   IF (ABS(CKA).LT.ABS(CKB)) THEN
		XSVI(I) = XA
	   ELSE
		XSVI(I) = XB
	   ENDIF

  100	CONTINUE    ! Continue loop on Y solutions
C
C  The following choses which of the two circle crossings in the "best"
C  based on the z difference of the trajectories at that point.
C
C  Loop over crossings

	DO 200 I=1,2

	   DZSVI(I) = DZMAX	! Set z difference = big as a flag

C  Xhat, Yhat = unit vector from circle center to vertex
	   XHAT = (XSVI(I) - X1)/R1
	   YHAT = (Y1 - YSVI(I))/R1

C  Turning angle.  Angular change from distance of closest
C  approach (x,y=0) and vertex
	   TANG = ATAN2(XHAT,YHAT) - PAR0(5,1)
	   IF (TANG .GT.  PI) TANG = TANG - 2.*PI   ! Put -Pi < TANG < Pi
	   IF (TANG .LT. -PI) TANG = TANG + 2.*PI
C  Don't allow turning angles beyond 81 degrees
	   IF (ABS(TANG) .GT. 0.45*PI) GO TO 200

	   S = R1*TANG			! Projected distance along the track
	   Z1 = PAR0(3,1) + S*PAR0(1,1)	! Track z at this vertex
C
C  Xhat, Yhat = unit vector from circle center to vertex
	   XHAT = (XSVI(I) - X2)/R2
	   YHAT = (Y2 - YSVI(I))/R2

C  Turning angle.  Angular change from distance of closest
C  approach (x,y=0) and vertex
	   TANG = ATAN2(XHAT,YHAT) - PAR0(5,2)
	   IF (TANG .GT.  PI) TANG = TANG - 2.*PI   ! Put -Pi < TANG < Pi
	   IF (TANG .LT. -PI) TANG = TANG + 2.*PI
C  Don't allow turning angles beyond 81 degrees
	   IF (ABS(TANG) .GT. 0.45*PI) GO TO 200

	   S = R2*TANG			! Projected distance along the track
	   Z2 = PAR0(3,2) + S*PAR0(1,2)	! Track z at this vertex

C  DZSVI = z difference between the two trajectories figure of merit
	   DZSVI(I) = Z2 - Z1
	   ZSVI(I) = 0.5*(Z2+Z1)	! Initial z = halfway between tracks

  200	CONTINUE	! Continue loop on crossings

C  Set I=best solution (smaller z difference)
	IF (ABS(DZSVI(2)) .LT. ABS(DZSVI(1))) THEN
	    I = 2
	ELSE
	    I = 1
	ENDIF

C  Done, if chosen solution is acceptable set IERR to 0

	IF (ABS(DZSVI(I)) .LT. DZMAX) THEN
	    IERR=0
  	    XSV = XSVI(I)
	    YSV = YSVI(I)
	    ZSV = ZSVI(I)
	    DLZ = DZSVI(I)
	ELSE
	    IERR=-1
	    XSV = 0
	    YSV = 0
	    ZSV = 0
	    DLZ = 0
	END IF

	RETURN

  	END
