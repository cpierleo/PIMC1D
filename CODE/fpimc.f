       subroutine fpimc
       implicit none
       include 'pimc.par'
       integer il,l,js,jsm,i,k,ist
       integer lm,jm,jk,step,iflag,jl,lvl,jdn,jup,jmid,kdiv
       integer mover(mmovers),itspan,itar,iorg,kspan,ndiv
       real*8 ht,vcm(mdim),veloc(mdim,mslices),ek,ekn,eold,enew
       real*8 qr(mdim,mslices),qrold(mdim,mslices)
       real*8 eps,sg,xold(mdim,mmovers),xm
       real*8 deltac,cold,cnew,phi,qm,dum(mdim)
       real*8 ronf,ddq,ddx
       real*8 gaussn,harmonic
       real*4 etime,tim(2)
       parameter (eps=1.d-32)
       external gaussn,harmonic
       character*80 messagge
       include 'pimc.cm'
       data iflag/0/
 
c do one block of FPIMC

!      print *,'x:'
!      write(*,'(10f15.8)') (x(1,js),js=1,nslices)
       call reduce(x,qr,xtoq,mdim,nslices,mslices)
!      print *,'qr:'
!      write(*,'(10f15.8)') (qr(1,js),js=1,nslices)
!      call reduce(qr,x,qtox,mdim,nslices,mslices)
!      print *,'x:'
!      write(*,'(10f15.8)') (x(1,js),js=1,nslices)
      
       do step=1,nstep

        ddx=0.d0
        ddq=0.d0
        do l=1,ndim
         do js=1,nslices
          xold(l,js)=x(l,js)
          qrold(l,js)=qr(l,js)
         enddo
        enddo
!       do l=1,ndim
!        do js=1,nslices-1
!         jl=iwrap(js+1)
!         ddx=ddx+(x(l,js)-x(l,jl))**2
!         ddq=ddq+qr(l,js)**2/sigma2q(js)
!        enddo
!        ddx=ddx+(x(l,1)-x(l,nslices))**2
!       enddo
!       ddx=ddx*cke
!       ddq=ddq/2.d0
!       print *,'ddx,ddq',ddx,ddq
! sample the normal modes
! qr(nslices)=sqrt(nslices)*X_cm
        do l=1,ndim
         call rmar(ronf)
         qr(l,nslices)=qrold(l,nslices)+hop*(ronf-0.5d0)
        enddo
! now other NMs
        do l=1,ndim
         do i=1,nslices-1
          qr(l,i)=gaussn()*sigmaq(i)  !/sqrt(2.d0)
         enddo
        enddo
!      print *,'qr:'
!      write(*,'(10f15.8)') (qr(1,js),js=1,nslices)
! back to real space
        call reduce(qr,x,qtox,mdim,nslices,mslices)
!      print *,'x:'
!      write(*,'(10f15.8)') (x(1,js),js=1,nslices)
!       ddx=0.d0
!       ddq=0.d0
!       do l=1,ndim
!        do js=1,nslices-1
!         jl=iwrap(js+1)
!         ddx=ddx+(x(l,js)-x(l,jl))**2
!         ddq=ddq+qr(l,js)**2/sigma2q(js)
!        enddo
!        ddx=ddx+(x(l,1)-x(l,nslices))**2
!       enddo
!       ddx=ddx*cke
!       ddq=ddq/2.d0
!       print *,'ddx,ddq',ddx,ddq
!      call reduce(x,qr,xtoq,mdim,nslices,mslices)
!      print *,'qr:'
!      write(*,'(10f15.8)') (qr(1,js),js=1,nslices)
!      print *,'qi:'
!      write(*,'(10f15.8)') (qi(1,js),js=1,nslices)

!       pause
        deltac=0.d0
c computes the trial action
        if (omega.ne.0.d0) then
         cold=0.d0
         cnew=0.d0
         do js=1,nslices    !+1
          cold=cold+harmonic(xold(1,js),ndim,dum)
          cnew=cnew+harmonic(x(1,js),ndim,dum)
         enddo
         deltac=-cnew+cold
        endif
        ntry=ntry+1
        call rmar(ronf)
        ronf=dlog(max(eps,ronf))
        if ( deltac.ge.ronf) then
c    ** .. accept **
         nacc = nacc + 1
c update the action
         cpot=cpot-deltac
         action=cpot
        else
c    **  .. reject and replace **
         do js=1,nslices
          do l=1,ndim
           x(l,js) = xold(l,js)
           qr(l,js) = qrold(l,js)
          enddo
         enddo
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
