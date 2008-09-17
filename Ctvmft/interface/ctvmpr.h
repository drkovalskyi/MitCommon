C ========================= Include ctvmpr.inc =========================
C ============ include file for the profiling of the code ==============
C
      INTEGER NCAL00,NCALFA,NCALVF,NCALCF,NCALPF,NCALMF,NCALPR,NCAL01,
     +        NCALMC,NCALDC
      INTEGER TIME00,TIMEFA,TIMEVF,TIMECF,TIMEPF,TIMEMF,TIMEPR,TIME01,
     +        TIMEMC,TIMEDC

      INTEGER NCALTO
      INTEGER TIMETO

      COMMON /CTVMPR/
     +        NCAL00,NCALFA,NCALVF,NCALCF,NCALPF,NCALMF,NCALPR,NCAL01,
     +        NCALMC,NCALDC,NCALTO,
     +        TIME00,TIMEFA,TIMEVF,TIMECF,TIMEPF,TIMEMF,TIMEPR,TIME01,
     +        TIMEMC,TIMEDC,TIMETO
