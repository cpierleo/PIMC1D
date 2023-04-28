       subroutine reptate
       implicit none
       include 'pimc.par'
       integer il,l,js,jsm,i,k
       integer lm,jm,jk,step,iflag,jl,lvl,jdn,jup,jmid,kdiv
       integer mover(mmovers),itspan,itar,iorg,kspan,ndiv,itar1
       real*8 eps,sg,xold(mdim,mmovers),xm,fc
       real*8 deltac,cold,cnew,std,dum(mdim)
       real*8 ronf
       real*8 gaussn,harmonic,lpsit
       real*4 etime,tim(2)
       parameter (eps=1.d-32)
       external gaussn,harmonic,lpsit
       save itspan,iflag,lvl,std
       character*80 messagge
       include 'pimc.cm'
       data iflag/0/
 
c do one block of MC dynamics
c standard Reptation with single link sampling
c set the number of movers to a power of 2
       if (iflag.eq.0) then
        iflag=1
        nmovers=1
        itspan=nmovers+1
        write(*,*) 'effective number of movers=',nmovers
        write(2,*)
        write(2,*) 'effective number of movers=',nmovers
        std=dsqrt(2.d0)*sigma
       endif
       do step=1,nstep
        if (nmovers.gt.0) then
c select the sampling direction at random
         call rmar(ronf)
         if (ronf.le.0.5d0) then
          iorg=1              
          itar=nslices
          itar1=nslices-1
         else
          iorg=nslices
          itar=1
          itar1=2
         endif
         do l=1,ndim
          do js=1,nslices
           xold(l,js)=x(l,js)
          enddo
          call rmar(ronf)
          x(l,itar)=x(l,iorg)+std*gaussn()
         enddo
!        write (*,*) 'iorg,itar',iorg,itar
!        write (*,'(i5,2g20.10)') (js,xold(1,js),x(1,js),js=1,nslices)
!        pause
         
         deltac=0.d0
         if (omega.ne.0.d0) then
          if (nslices>1) then
           cold=lpsit(xold(1,iorg),model,tau)
     .         +lpsit(xold(1,itar),model,tau) 
           cnew=lpsit(x(1,itar),model,tau)
     .         +lpsit(xold(1,itar1),model,tau)
           cold=cold+0.5d0*( harmonic(xold(1,itar),ndim,dum) 
     .                      +harmonic(xold(1,itar1),ndim,dum) )
           cnew=cnew+0.5d0*( harmonic(x(1,itar),ndim,dum)
     .                      +harmonic(xold(1,iorg),ndim,dum) )
          else
           cold=2.d0*lpsit(xold(1,1),model,tau)
           cnew=2.d0*lpsit(x(1,1),model,tau)
          endif
          deltac=cold-cnew
         endif
         ntry=ntry+1
         call rmar(ronf)
         ronf=dlog(max(eps,ronf))
         if ( deltac.ge.ronf) then
c    ** .. accept **
!         write (*,*) 'accepted'
          nacc = nacc + 1
          do l=1,ndim
           xold(l,itar)=x(l,itar)
           do js=1,nslices
            if(iorg.eq.1) then
             jm=iwrap(js-1)
            elseif(iorg.eq.nslices) then
             jm=iwrap(js+1)
            else
             stop 'iorg.ne.1.and.iorg.ne.nslices'
            endif
            x(l,js)=xold(l,jm)
           enddo
          enddo
c update the action
          cpot=cpot-deltac
!         cpot=0.d0
!         do js=1,nslices
!          fc=1.d0
!          if(js.eq.1.or.js.eq.nslices) fc=0.5d0
!          cpot=cpot+fc*harmonic(x(1,js),ndim)
!         enddo
          action=cpot
         else
c    **  .. reject and replace **
!         write (*,*) 'rejected'
          do l=1,ndim
           x(l,itar) = xold(l,itar)
          enddo
         endif  ! accept or reject
        endif   ! (nmovers>0)
!       write (*,'(i5,g20.10)') (js,x(1,js),js=1,nslices)
!       pause
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
