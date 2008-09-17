c $Id:$
#include "MitCommon/Ctvmft/interface/dimensions.hh"

       subroutine gettrk(jtrk,bnk,mass,radius,beam,vrtx,helix,em,ierr)
c
c    Routine to get track parameters stored in common.
c  
      implicit none
#include "MitCommon/Ctvmft/interface/ctvmdi.h"
#include "MitCommon/Ctvmft/interface/ctvmft.h"
#include "MitCommon/Ctvmft/interface/ctvmtr.h"
c
      integer     jtrk,beam,vrtx,ierr
      real        mass,radius,helix(5),em(5,5)
      character*4 bnk
c      
      integer j,k,jt
c
      ierr=0
      do jt=1,ntrack
         if(jtrk.eq.list(jt)) then
	    ierr=0
	    do j=1,5
	       helix(j)=trhelix(j,jt)
	       do k=1,5
	          em(j,k)=trem(j,k,jt)
	       enddo
	    enddo
ccc	    write(*,'('' Helix: '',5(1x,e15.6))') helix
ccc	    write(*,'(5(6x,5(1x,e15.5),/))') em
	    return
	 endif
      enddo
      ierr=1
c
      return
      end
