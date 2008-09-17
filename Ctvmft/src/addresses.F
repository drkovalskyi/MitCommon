c $Id:$
#include "MitCommon/Ctvmft/interface/dimensions.hh"

c-----------------------------------------------------------------------
      integer function ctvmfr_address()
c
c    Get address of CTVMFR common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/ctvmdi.h"
#include "MitCommon/Ctvmft/interface/ctvmft.h"

      integer address
      
      ctvmfr_address = address(vmat(1,1))

      return
      end

c-----------------------------------------------------------------------
      integer function ctvmq_address()
c
c    Get address of CTVMQ common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/ctvmdi.h"
#include "MitCommon/Ctvmft/interface/ctvmft.h"

C  CTC first approximation and track-vertex step parameters for CTVMFT
      DATA  DRMAX /  2.0 /  ! max XY separ., non-intersecting circles
      DATA  DZMAX / 20.0 /  ! max allowed track Z separation at the vtx
      DATA  RVMAX / 70.0 /  ! max allowed vertex radius
      DATA TRNMAX /  0.5 /  ! max permitted track turning angle to vtx
C      DATA  DSMIN / -2.0 /  ! min permitted arc length, track to vtx
      DATA  DSMIN / -99.0 / ! min permitted arc length, track to vtx,CMS
      DATA   CVTX /MAXVTX*0/
      
      integer address
      
      ctvmq_address = address(runnum)
      
      return
      end

c-----------------------------------------------------------------------
      integer function ctvmtu_address()
c
c    Get address of CTVMTU common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/ctvmtu.h"

      integer address

      DATA    C20MAX / -1.0  /
      DATA    C21MAX / -1.0  /
      DATA    LXY0MN / -1e10 /
      DATA    LXY1MN / -1e10 /
      
      ctvmtu_address = address(C20MAX)

      return
      end

c-----------------------------------------------------------------------
      integer function fiddle_address()
c
c    Get address of FIDDLE common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/ctvmfi.h"      
      integer address

      fiddle_address = address(excuse)

      return
      end

c-----------------------------------------------------------------------
      integer function trkprm_address()
c
c    Get address of TRKPRM common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/ctvmdi.h"
#include "MitCommon/Ctvmft/interface/ctvmtr.h"

      integer address
      
      trkprm_address = address(trhelix(1,1))

      return
      end
