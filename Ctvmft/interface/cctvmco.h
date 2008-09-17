C --------- Include file to define the most commn constants ------------
C  Definitions of pi pi/ radian, degree etc.

C---- the following include should be used for all definitions of pi

      REAL      PI,TWOPI,HALFPI,RADDEG,DEGRAD
      PARAMETER(PI     = 3.141592653)  
      PARAMETER(TWOPI  = 2.0*PI)    
      PARAMETER(HALFPI = 0.5*PI)   
      PARAMETER(RADDEG = 180./PI) 
      PARAMETER(DEGRAD = PI/180.) 
