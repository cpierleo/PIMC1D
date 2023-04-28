      subroutine wdynamics
c performe a single trajectory for the doubled phase space
      include 'pimc.par'
      integer js,ist,inc,k,ndyn,nstat,icall
      real*8 t0,t1,t2,tdyn,tstat,tp,td,tlocal
      save tdyn,tstat,nstat,inc
      include 'pimc.cm'
      data icall/0/

      if (icall.eq.0) then
       icall=1
       inc=nwstep/iwaver+1
       tdyn=0.0
       tstat=0.0
       nstat=0
      endif

      nstat=nstat+1
      tlocal=0.
!     call second(t0)
      do js=1,2
       call forcew(js)
      enddo
      call quanti(psign,0)
      do ist=1,nwstep
       do js=1,2
        call velverl(js)
       enddo
       if (mod(ist,iwaver).eq.0) then
!       call second(tp)
        call quanti(psign,ist/iwaver)
!       call second(td)
        tlocal=tlocal+(td-tp)
       endif
      enddo
!     call second(t1)
      do js=1,2
       call cumul1(ekin(0,js),avp(iwkin(js)),anormp(iwkin(js)),inc)
       call cumul1(epot(0,js),avp(iwpot(js)),anormp(iwpot(js)),inc)
       call cumul1(et(0,js),avp(iwetot(js)),anormp(iwetot(js)),inc)
       call cumul1(wtemp(0,js),avp(iwtemp(js)),anormp(iwtemp(js)),inc)
      enddo
!     do k=1,nlamb+jaddk
!      sktw0(k)=sktw(0,k)
!     enddo
!     call cumul(sktw0,avwig(ifk0),anwig(ifk0),nlamb+jaddk)
!     do k=1,nlamb+jaddk
!      call cumul(sktw(0,k),avwig(ifktw(k)),anwig(ifktw(k)),inc)
!     enddo
      call cumul(cpp(0),avp(iwcpp),anormp(iwcpp),inc)
!     call second(t2)
      call wsprint(inc,nstat)
c
      tdyn=tdyn+(t1-t0)-tlocal
      tstat=tstat+(t2-t1)+tlocal
      write (69,121) nstat,tdyn,tstat,(tdyn+tstat)
121   format('finished trajectory',i5,/
     &      ,'evolution time (sec.)   = ',f15.7,/
     &      ,'stat. aver. time (sec.) = ',f15.7,/
     &      ,'total CPU time (sec.)   = ',f15.7)
      rewind(69)

      return
      end
