       subroutine hybrid_mc
       implicit none
       include 'pimc.par'
       integer il,l,js,jsm,i,k,ist
       integer lm,jm,jk,step,iflag,jl,lvl,jdn,jup,jmid,kdiv
       integer mover(mmovers),itspan,itar,iorg,kspan,ndiv
       real*8 ht,ht2,vcm(mdim),veloc(mdim,mslices),ek,ekn,eold,enew
       real*8 eps,sg,xold(mdim,mmovers),xm
       real*8 deltac,cold,cnew
       real*8 ronf
       real*8 gaussn,harmonic
       real*4 etime,tim(2)
       parameter (eps=1.d-32)
       external gaussn,harmonic
       save ht,ht2,iflag
       character*80 messagge
       include 'pimc.cm'
       data iflag/0/
 
c do one block of Hybrib MC dynamics
       if (iflag.eq.0) then
        iflag=1
        ht=hop_hmc/nstdyn
        ht2=0.5d0*ht/amass
!       write (*,*) 'hybrid_mc: ht,ht2',ht,ht2
       endif

       do step=1,nstep
        call hforce     ! compute forces and energy
!       write (*,*)((x(l,js),force(l,js),l=1,ndim),js=1,nslices)
!       write(*,*) vcl
! sample the velocities
        ek=0.d0
        do l=1,ndim
         vcm(l)=0.d0
         do i=1,nslices
          veloc(l,i)=gaussn()*sqrt(temp)
          vcm(l)=vcm(l)+veloc(l,i)
         enddo
         vcm(l)=vcm(l)/nslices
         do i=1,nslices
!         veloc(l,i)=veloc(l,i)-vcm(l)
          ek=ek+veloc(l,i)**2
         enddo
        enddo
        ek=0.5d0*amass*ek
        eold=ek+vcl
        rewind(91)
!     write (*,777) 0,vcm(1),eold,ek,vcl
777   format(5g15.7)

        do l=1,ndim
         do js=1,nslices
          xold(l,js)=x(l,js)
         enddo
        enddo
c dynamics
        do ist=1,nstdyn
         do i=1,nslices
          do l=1,ndim
           veloc(l,i)=veloc(l,i)+ht2*force(l,i)
           x(l,i)=x(l,i)+ht*veloc(l,i)
          enddo
         enddo
         call hforce
         ekn=0.d0
         do l=1,ndim
          vcm(l)=0.d0
          do i=1,nslices
           veloc(l,i)=veloc(l,i)+ht2*force(l,i)
           vcm(l)=vcm(l)+veloc(l,i)
           ekn=ekn+veloc(l,i)**2
          enddo
          vcm(l)=vcm(l)/nslices
         enddo
         ekn=0.5d0*amass*ekn
         enew=ekn+vcl
!        if(mod(ist,1).eq.0) write (*,777) 
!    &     ist*ht,vcm(1),enew,ekn,vcl
        enddo       !  over MD steps

!       call flush(91)
        deltac=beta*(eold-enew)
        ntry=ntry+1
        call rmar(ronf)
        ronf=dlog(max(eps,ronf))
        if ( deltac.ge.ronf) then
c    ** .. accept **
         nacc = nacc + 1
c update the action
         cpot=vcl
         action=vcl
        else
c    **  .. reject and replace **
         do js=1,nslices
          do l=1,ndim
           x(l,js) = xold(l,js)
          enddo
         enddo
        endif  ! accept or reject
c attempt to displace the entire path
        call rmar(ronf)
        if (ronf.le.gammadsp) call displace
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
