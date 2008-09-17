c $Id:$
#include "MitCommon/Ctvmft/interface/cdimensions.hh"

c-----------------------------------------------------------------------
      integer function cctvmfr_address()
c
c    Get address of CTVMFR common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"

      integer address
      
      cctvmfr_address = address(vmat(1,1))

      return
      end

c-----------------------------------------------------------------------
      integer function cctvmq_address()
c
c    Get address of CTVMQ common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmft.h"

C  CTC first approximation and track-vertex step parameters for CTVMFT
      DATA  DRMAX /  2.0 /  ! max XY separ., non-intersecting circles
      DATA  DZMAX / 20.0 /  ! max allowed track Z separation at the vtx
      DATA  RVMAX / 70.0 /  ! max allowed vertex radius
      DATA TRNMAX /  0.5 /  ! max permitted track turning angle to vtx
C      DATA  DSMIN / -2.0 /  ! min permitted arc length, track to vtx
      DATA  DSMIN / -101.0 /! min permitted arc length, track to vtx
      DATA   CVTX /MAXVTX*0/
      
      integer address
      
      cctvmq_address=address(runnum)
      
      return
      end

c-----------------------------------------------------------------------
      integer function cctvmtu_address()
c
c    Get address of CTVMTU common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/cctvmtu.h"

      integer address

      DATA    C20MAX / -1.0  /
      DATA    C21MAX / -1.0  /
      DATA    LXY0MN / -1e10 /
      DATA    LXY1MN / -1e10 /
      
      cctvmtu_address = address(C20MAX)

      return
      end

c-----------------------------------------------------------------------
      integer function cfiddle_address()
c
c    Get address of FIDDLE common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/cctvmfi.h"      
      integer address

      cfiddle_address = address(excuse)

      return
      end

c-----------------------------------------------------------------------
      integer function ctrkprm_address()
c
c    Get address of TRKPRM common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/cctvmdi.h"
#include "MitCommon/Ctvmft/interface/cctvmtr.h"

      integer address
      
      ctrkprm_address = address(trhelix(1,1))

      return
      end
