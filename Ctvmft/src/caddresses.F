cc $Id:$
c#include "MitCommon/Ctvmft/interface/cdimensions.hh"
c
cc-----------------------------------------------------------------------
c      integer function cctvmfr_address()
cc
cc    Get address of CTVMFR common
cc
cc-----------------------------------------------------------------------
c      implicit none
c#include "MitCommon/Ctvmft/interface/cctvmdi.h"
c#include "MitCommon/Ctvmft/interface/cctvmft.h"
c
c      integer address
c      
c      cctvmfr_address = address(vmat(1,1))
c
c      return
c      end
c
cc-----------------------------------------------------------------------
c      integer function cctvmq_address()
cc
cc    Get address of CTVMQ common
cc
cc-----------------------------------------------------------------------
c      implicit none
c#include "MitCommon/Ctvmft/interface/cctvmdi.h"
c#include "MitCommon/Ctvmft/interface/cctvmft.h"
c
cC  CTC first approximation and track-vertex step parameters for CTVMFT
c      DATA  DRMAX /  2.0 /  ! max XY separ., non-intersecting circles
c      DATA  DZMAX / 20.0 /  ! max allowed track Z separation at the vtx
c      DATA  RVMAX / 70.0 /  ! max allowed vertex radius
c      DATA TRNMAX /  0.5 /  ! max permitted track turning angle to vtx
cC      DATA  DSMIN / -2.0 /  ! min permitted arc length, track to vtx
c      DATA  DSMIN / -101.0 /! min permitted arc length, track to vtx
c      DATA   CVTX /MAXVTX*0/
c      
c      integer address
c      
c      cctvmq_address=address(runnum)
c      
c      return
c      end
c
cc-----------------------------------------------------------------------
c      integer function cctvmtu_address()
cc
cc    Get address of CTVMTU common
cc
cc-----------------------------------------------------------------------
c      implicit none
c#include "MitCommon/Ctvmft/interface/cctvmtu.h"
c
c      integer address
c
c      DATA    C20MAX / -1.0  /
c      DATA    C21MAX / -1.0  /
c      DATA    LXY0MN / -1e10 /
c      DATA    LXY1MN / -1e10 /
c      
c      cctvmtu_address = address(C20MAX)
c
c      return
c      end
c
cc-----------------------------------------------------------------------
c      integer function cfiddle_address()
cc
cc    Get address of FIDDLE common
cc
cc-----------------------------------------------------------------------
c      implicit none
c#include "MitCommon/Ctvmft/interface/cctvmfi.h"      
c      integer address
c
c      cfiddle_address = address(excuse)
c
c      return
c      end
c
cc-----------------------------------------------------------------------
c      integer function ctrkprm_address()
cc
cc    Get address of TRKPRM common
cc
cc-----------------------------------------------------------------------
c      implicit none
c#include "MitCommon/Ctvmft/interface/cctvmdi.h"
c#include "MitCommon/Ctvmft/interface/cctvmtr.h"
c
c      integer address
c      
c      ctrkprm_address = address(trhelix(1,1))
c
c      return
c      end
