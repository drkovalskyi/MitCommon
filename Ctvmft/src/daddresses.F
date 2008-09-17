c $Id:$
#include "MitCommon/Ctvmft/interface/ddimensions.hh"

c-----------------------------------------------------------------------
      integer function dctvmfr_address()
c
c    Get address of CTVMFR common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/dctvmdi.h"
#include "MitCommon/Ctvmft/interface/dctvmft.h"

      integer address
      
      dctvmfr_address = address(vmat(1,1))

      return
      end

c-----------------------------------------------------------------------
      integer function dctvmq_address()
c
c    Get address of CTVMQ common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/dctvmdi.h"
#include "MitCommon/Ctvmft/interface/dctvmft.h"

C  CTC first approximation and track-vertex step parameters for CTVMFT
      DATA  DRMAX /  2.0 /  ! max XY separ., non-intersecting circles
      DATA  DZMAX / 20.0 /  ! max allowed track Z separation at the vtx
      DATA  RVMAX / 70.0 /  ! max allowed vertex radius
      DATA TRNMAX /  0.5 /  ! max permitted track turning angle to vtx
C      DATA  DSMIN / -2.0 /  ! min permitted arc length, track to vtx
      DATA  DSMIN / -50.0 / ! min permitted arc length, track to vtx
      DATA   CVTX /MAXVTX*0/
      
      integer address
      
      dctvmq_address=address(runnum)
      
      return
      end

c-----------------------------------------------------------------------
      integer function dctvmtu_address()
c
c    Get address of CTVMTU common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/dctvmtu.h"

      integer address

      DATA    C20MAX / -1.0  /
      DATA    C21MAX / -1.0  /
      DATA    LXY0MN / -1e10 /
      DATA    LXY1MN / -1e10 /
      
      dctvmtu_address = address(C20MAX)

      return
      end

c-----------------------------------------------------------------------
      integer function dfiddle_address()
c
c    Get address of FIDDLE common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/dctvmfi.h"      
      integer address

      dfiddle_address = address(excuse)

      return
      end

c-----------------------------------------------------------------------
      integer function dtrkprm_address()
c
c    Get address of TRKPRM common
c
c-----------------------------------------------------------------------
      implicit none
#include "MitCommon/Ctvmft/interface/dctvmdi.h"
#include "MitCommon/Ctvmft/interface/dctvmtr.h"

      integer address
      
      dtrkprm_address = address(trhelix(1,1))

      return
      end
