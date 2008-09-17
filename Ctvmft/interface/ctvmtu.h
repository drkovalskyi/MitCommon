C ========================= Include ctvmtu.inc =========================
C ============ include file for fine tuning convergence parameters =====
C
      REAL    C20MAX ! Chi2Mx cut after 0. iteration (above is failure)
      REAL    C21MAX ! Chi2Mx cut after 1. iteration (above is failure)
      REAL    LXY0MN ! LxyMin cut after 0. iteration (below is failure)
      REAL    LXY1MN ! LxyMin cut after 1. iteration (below is failure)     
      COMMON /CTVMTU/ C20MAX,C21MAX,LXY0MN,LXY1MN
