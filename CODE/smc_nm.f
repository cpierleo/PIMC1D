       subroutine smc_nm
       implicit none
       include 'pimc.par'
       integer il,l,js,jsm,i,k,ist,ip1,im1,js1,nsl2
       integer lm,jm,jk,step,jl,lvl,jdn,jup,jmid,kdiv
       integer mover(mmovers),itspan,itar,iorg,kspan,ndiv
       real*8 ht,ek,ekn,eold,enew,weight1
       real*8 fp(mdim,mslices)
       real*8 fpold(mdim,mslices),ft2(mslices,2),ftot(2)
       real*8 phiK(mdim,mslices,2)
       real*8 phiU(mdim,mslices)
       real*8 eps,sg,xold(mdim,mmovers),xm
       real*8 deltac,cold,cnew,dd2(2),dq(mdim),weight
       real*8 qr(mdim,mslices),qrold(mdim,mslices)
       real*8 ronf,ddq(2)
       real*8 gaussn,harmonic
       real*4 etime,tim(2)
       parameter (eps=1.d-32)
       external gaussn,harmonic
       character*80 messagge
       include 'pimc.cm'
 
c do one block of SMC
!      print 'x:'
!      write(*,'(10f15.8)') (x(1,js),js=1,nslices)
       call reduce(x,qr,xtoq,mdim,nslices,mslices)
!      print *,'qr:'
!      write(*,'(10f15.8)') (qr(1,js),js=1,nslices)
!      call reduce(qr,x,qtox,mdim,nslices,mslices)
!      print *,'x:'
!      write(*,'(10f15.8)') (x(1,js),js=1,nslices)

       do step=1,nstep

        dd2(1)=0.d0
        ddq(1)=0.d0
        deltac=0.d0
        do js=1,nslices
         do l=1,ndim
          xold(l,js)=x(l,js)
          qrold(l,js)=qr(l,js)
         enddo
        enddo
        do l=1,ndim
         do js=1,nslices-1
          phiK(l,js,1)=-qr(l,js)/sigma2q(js)     ! real part of the kinetic force in the nm space
         enddo
         phiK(l,nslices,1)=0.d0
         do js=1,nslices-1
          ip1=iwrap(js+1)
          dd2(1)=dd2(1)+(x(l,js)-x(l,ip1))**2     ! kinetic action in coordinate space
          ddq(1)=ddq(1)+qr(l,js)**2/sigma2q(js)
         enddo
         dd2(1)=dd2(1)+(x(l,nslices)-x(l,1))**2
        enddo
        dd2(1)=dd2(1)*cke
        ddq(1)=0.5d0*ddq(1)
!       print *,'ddk,ddq',dd2,ddq
!       stop
        if (omega.ne. 0.d0) then
         cold=0.d0
         do js=1,nslices    !+1
          cold=cold+harmonic(xold(1,js),ndim,fpold(1,js))
         enddo
         call reduceT(fpold,phiU,qtox,mdim,nslices,mslices)  
         do js=1,nslices
          do l=1,ndim
           phiK(l,js,1)=phiK(l,js,1)+phiU(l,js)    ! total force
          enddo
         enddo
        endif

! dynamics x'=x+Df+gaussn()*sqrt(2D)
        do js=1,nslices-1
         do l=1,ndim
          qr(l,js)=qr(l,js)+gaussn()*sqrt(2.d0*hop) 
     &            + hop*phiK(l,js,1)
!         qr(l,js)=qr(l,js)+gaussn()*sqrt(4.d0*sigma2q(js)*hop) 
!    &            + 2.d0*sigma2q(js)*hop*phiK(l,js,1)
         enddo
        enddo
! CoM dynamics
        do l=1,ndim
         qr(l,nslices)=qr(l,nslices)+gaussn()*sqrt(2.d0*hop_cm)
     &                +hop_cm*phiK(l,nslices,1)
        enddo

c compute the reverse forces and the ideal action
        dd2(2)=0.d0
        do l=1,ndim
         do js=1,nslices-1
          phiK(l,js,2)=-qr(l,js)/sigma2q(js)
         enddo
         phiK(l,nslices,2)=0.d0
        enddo
        call reduce(qr,x,qtox,mdim,nslices,mslices)
        ddq(2)=0.d0
        do l=1,ndim
         do js=1,nslices-1
          ip1=iwrap(js+1)
          dd2(2)=dd2(2)+(x(l,js)-x(l,ip1))**2
          ddq(2)=ddq(2)+qr(l,js)**2/sigma2q(js)
         enddo
         dd2(2)=dd2(2)+(x(l,nslices)-x(l,1))**2
        enddo
        dd2(2)=dd2(2)*cke
        ddq(2)=0.5d0*ddq(2)
!       print *,'ddk,ddq',dd2,ddq
!       stop
        if (omega.ne.0.d0) then
         cnew=0.d0
         do js=1,nslices    !+1
          cnew=cnew+harmonic(x(1,js),ndim,fp(1,js))
         enddo
         call reduceT(fp,phiU,qtox,mdim,nslices,mslices)  ! real part of the ext force
         do js=1,nslices
          do l=1,ndim
           phiK(l,js,2)=phiK(l,js,2)+phiU(l,js)
          enddo
         enddo
         deltac=-cnew+cold
!        print *,'deltac=',deltac
        endif
!       pause

        ntry=ntry+1
        call rmar(ronf)
        ronf=dlog(max(eps,ronf))
        weight=-dd2(2)+dd2(1)+deltac
        weight1=0.d0
!       print *,'weight=',weight
        do js=1,nslices-1
         ft2(js,1)=0.d0
         ft2(js,2)=0.d0
         do l=1,ndim
          dq(l)=qr(l,js)-qrold(l,js)
          ft2(js,1)=ft2(js,1)+phiK(l,js,1)**2
          ft2(js,2)=ft2(js,2)+phiK(l,js,2)**2
         enddo
         weight1=weight1-(qr(1,js)**2-qrold(1,js)**2)/sigma2q(js)/2.d0
     &                *hop/2.d0/sigma2q(js)
         weight=weight+0.25d0*hop*(ft2(js,1)-ft2(js,2))
!        weight=weight+0.25d0*2.d0*sigma2q(js)*hop*(ft2(js,1)-ft2(js,2))
         do l=1,ndim
          weight=weight-0.5d0*dq(l)*(phiK(l,js,1)+phiK(l,js,2))
         enddo
        enddo
! add CoM now
        ft2(nslices,1)=0.d0
        ft2(nslices,2)=0.d0
        do l=1,ndim
         dq(l)=qr(l,nslices)-qrold(l,nslices)
         ft2(nslices,1)=ft2(nslices,1)+phiK(l,nslices,1)**2
         ft2(nslices,2)=ft2(nslices,2)+phiK(l,nslices,2)**2
        enddo
        weight=weight+0.25d0*hop_cm*(ft2(nslices,1)-ft2(nslices,2))
        do l=1,ndim
         weight=weight-0.5d0*dq(l)*(phiK(l,nslices,1)+phiK(l,nslices,2))
        enddo
!       if (abs(weight-weight1).gt.0.001) then
!        print *,'weight=',weight,weight1
!        pause
!       endif

        if ( weight.ge.ronf) then
c    ** .. accept **
         nacc = nacc + 1
c update the action
         cpot=cpot-deltac
         action=cpot!-dd2(1)+dd2(2)
        else
c    **  .. reject and replace **
         do js=1,nslices
          do l=1,ndim
           x(l,js) = xold(l,js)
           qr(l,js) = qrold(l,js)
          enddo
         enddo
!        call reduce(qr,x,qtox,mdim,nslices,mslices)
        endif  ! accept or reject

c attempt to displace the entire path
!       call rmar(ronf)
!       if (ronf.le.gammadsp) call displace
!       if (mod(step,iwigner).eq.0) call wigner
c computing averages
        if (mod(step,nanal).eq.0) then
         tt=etime(tim)
         call analyze(1,step)
         ttanal=ttanal+etime(tim)-tt
        endif
        if (mod(step,nspill).eq.0)
     &      write (31) ((x(l,js),l=1,ndim),js=1,nslices)

       enddo         ! over steps

       return
       end
