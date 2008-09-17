C ========================= Include ctvmtr.inc =========================
C ============ include file for the exchange of track parameters =======
C
      REAL    TRHELIX  (5,MAXTRK) ! Track helix parameters
      REAL    TREM   (5,5,MAXTRK) ! Track Helix parameter errors matrix
      COMMON /CTVMTR/ TRHELIX,TREM
