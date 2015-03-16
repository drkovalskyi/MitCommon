cc $Id:$
c#include "MitCommon/Ctvmft/interface/dimensions.hh"
c
cc-----------------------------------------------------------------------
c      integer function ctvmfr_address()
cc
cc    Get address of CTVMFR common
cc
cc-----------------------------------------------------------------------
c      implicit none
c#include "MitCommon/Ctvmft/interface/ctvmdi.h"
c#include "MitCommon/Ctvmft/interface/ctvmft.h"
c
c      integer address
c      
c      ctvmfr_address = address(vmat(1,1))
c
c      return
c      end
c
cc-----------------------------------------------------------------------
c      integer function ctvmq_address()
cc
cc    Get address of CTVMQ common
cc
cc-----------------------------------------------------------------------
c      implicit none
c#include "MitCommon/Ctvmft/interface/ctvmdi.h"
c#include "MitCommon/Ctvmft/interface/ctvmft.h"
c
cC  CTC first approximation and track-vertex step parameters for CTVMFT
c      DATA  DRMAX /  2.0 /  ! max XY separ., non-intersecting circles
c      DATA  DZMAX / 20.0 /  ! max allowed track Z separation at the vtx
c      DATA  RVMAX / 70.0 /  ! max allowed vertex radius
c      DATA TRNMAX /  0.5 /  ! max permitted track turning angle to vtx
cC      DATA  DSMIN / -2.0 /  ! min permitted arc length, track to vtx
c      DATA  DSMIN / -99.0 / ! min permitted arc length, track to vtx,CMS
c      DATA   CVTX /MAXVTX*0/
c      
c      integer address
c      
c      ctvmq_address = address(runnum)
c      
c      return
c      end
c
cc-----------------------------------------------------------------------
c      integer function ctvmtu_address()
cc
cc    Get address of CTVMTU common
cc
cc-----------------------------------------------------------------------
c      implicit none
c#include "MitCommon/Ctvmft/interface/ctvmtu.h"
c
c      integer address
c
c      DATA    C20MAX / -1.0  /
c      DATA    C21MAX / -1.0  /
c      DATA    LXY0MN / -1e10 /
c      DATA    LXY1MN / -1e10 /
c      
c      ctvmtu_address = address(C20MAX)
c
c      return
c      end
c
cc-----------------------------------------------------------------------
c      integer function fiddle_address()
cc
cc    Get address of FIDDLE common
cc
cc-----------------------------------------------------------------------
c      implicit none
c#include "MitCommon/Ctvmft/interface/ctvmfi.h"      
c      integer address
c
c      fiddle_address = address(excuse)
c
c      return
c      end
c
cc-----------------------------------------------------------------------
c      integer function trkprm_address()
cc
cc    Get address of TRKPRM common
cc
cc-----------------------------------------------------------------------
c      implicit none
c#include "MitCommon/Ctvmft/interface/ctvmdi.h"
c#include "MitCommon/Ctvmft/interface/ctvmtr.h"
c
c      integer address
c      
c      trkprm_address = address(trhelix(1,1))
c
c      return
c      end
