      subroutine zeroav(ifl)
      implicit none
      include 'pimc.par'
      integer i,ifl,inc,l
      include 'pimc.cm'
c  initialization
      if (ifl.eq.0) then
       ikin=2
       iev=ikin+2
       iepot=iev+2
       ietot=iepot+14
       ievar=ietot+4
       irg=ievar+4
       ix0xb2=irg+4
       ieff=ix0xb2+2
       irho1=ieff+2*nslices+1
       irho2=irho1+2*nhist+1
       idel=irho2+2*nhist+1
       iinst=idel+2
       idelhist=iinst+2
       iinsthist=idelhist+2*nhist+1
       iepcor=iinsthist+2*nhist+1
       iepcor1=iepcor+2*nslices+1
       nav=iepcor1+2*nslices+1
       print *,'zeroav: iepcor,iepcor1,iekcor,iekcor1'
     &        ,iepcor,iepcor1,iekcor,iekcor1
!      ipr(1)=irho1+2*nhist+1
!      do l=2,nhist
!       ipr(l)=ipr(l-1)+2*nhist+1
!      enddo
!      ipi(1)=ipr(nhist)+2*nhist+1
!      do l=2,nhist
!       ipi(l)=ipi(l-1)+2*nhist+1
!      enddo
!      if (iwigner.eq.0) then
!       nav=ipi(nhist)+2*nhist+1
!      else
!       inc=(nwstep/iwaver)+1
!       iwkin(1)=ipi(nhist)+2*nhist+1
!       iwkin(2)=iwkin(1)+2*inc
!       iwpot(1)=iwkin(2)+2*inc
!       iwpot(2)=iwpot(1)+2*inc
!       iwetot(1)=iwpot(2)+2*inc
!       iwetot(2)=iwetot(2)+2*inc
!       iwtemp(1)=iwetot(2)+2*inc
!       iwtemp(2)=iwtemp(1)+2*inc
!       iwcpp=iwtemp(2)+2*inc
!       nav=iwcpp+2*inc+1
!      endif
       write (*,*) 'dimension of the cumulator vector: nav=',nav,mav
       if (nav.ge.mav) then
        write (*,*) 'nav.ge.mav',nav,mav
        stop
       endif
       do i=1,nav
        av(i)=0.d0
        anorm(i)=0.d0
       enddo
      else
       do i=1,nav
        avp(i)=0.d0
        anormp(i)=0.d0
       enddo
      endif
      return
      end
