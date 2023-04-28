       subroutine levy
       implicit none
       include 'pimc.par'
       integer il,l,js,jsm,jsp,i,is,mi,ik,jkmin,k!,idbg
       integer lm,jm,jk,step,iflag,jl,lvl,jdn,jup,jmid,kdiv,ifl,ifr
       integer mover(mmovers),itspan,itar,iorg,kspan,ndiv,js1,is1
       real*8 eps,sg,xold(mdim,0:mmovers),pold,pnew,xm,ppo,ppn,pwrap
       real*8 deltac,cold,cnew,dum(mdim)
       real*8 dio,din,ysold,ysnew,pold1,pnew1,ronf,deltaew,ewold,ewnew
       real*8 dbgact,dbgepact,dbgepact1,dbgep12,dbgcpot,epold,epnew
       real*8 gaussn,harmonic
       real*4 etime,tim(2)
       character*5 flag
       parameter (eps=1.d-32)
       external gaussn,coulomb
       save itspan,iflag,lvl,flag
       character*80 messagge
       include 'pimc.cm'
       data iflag/0/
c do one block of MC dynamics
c Levy path reconstruction by bisection
c set the number of movers to a power of 2
       if (iflag.eq.0) then
        iflag=1
        lvl=dint(dlog(dfloat(nmovers))/dlog(2.d0))
        nmovers=2**lvl-1
        itspan=nmovers+1
        write(*,*) 'effective number of movers=',nmovers
        write(2,*)
        write(2,*) 'effective number of movers=',nmovers
        if (nmovers > mmovers) stop 'nmovers > mmovers'
        nrap=nslices/max(1,nmovers)
        flag='#conf'
       endif
       do step=1,nstep
        if (nmovers.gt.0) then
c select the origin at random
         call rmar(ronf)
         iorg=dint(nslices*ronf)+1
c loops recursively over time origins
          do lm=1,nrap
           iorg=iwrap(iorg+nmovers)
           itar=iwrap(iorg+itspan)
           js=iorg
!          do l=1,ndim
!           xold(l,0)=x(l,js)
!          enddo
c loops over levels
           do jm=1,nmovers    !+1
            mover(jm)=iwrap(js+jm)
            do l=1,ndim
             xold(l,jm)=x(l,mover(jm))
            enddo
c           write (*,*) js,iwrap(js+itspan),mover(jm)
           enddo
           kspan=itspan
           ndiv=1
           do jl=lvl,1,-1
            jdn=js
            kspan=kspan/2
            sg=dsqrt(dfloat(kspan))*sigma
            do kdiv=1,ndiv
             jmid=iwrap(jdn+kspan)
             jup=iwrap(jmid+kspan)
             do l=1,ndim
c compute the trial free particle term
              xm=(x(l,jdn)+x(l,jup))/2.d0
              x(l,jmid)=xm+sg*gaussn()
             enddo
             jdn=jup
            enddo
            ndiv=2*ndiv
           enddo
           deltac=0.d0
c computes the trial action
           if (omega.ne.0.d0) then
            cold=0.d0
            cnew=0.d0
            do jm=1,nmovers    !+1
             js=mover(jm)
!            js1=iwrap(js-1)
             cold=cold+harmonic(xold(1,jm),ndim,dum)
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
            do jm=1,nmovers
             js=mover(jm)
             do l=1,ndim
              x(l,js) = xold(l,jm)
             enddo
            enddo
           endif  ! accept or reject
          enddo !  over nrap
        endif   ! (nmovers>0)
c attempt to displace the entire path
        call rmar(ronf)
        if (ronf.le.gammadsp) call displace
c computing averages
        if (mod(step,nanal).eq.0) then
         tt=etime(tim)
!        if (idbg.ge.2) then
!         messagge='levy before anal: step'
!         call compactdbg(messagge,step)
!        endif
         call analyze(1,step)
         if(model==2) call reaction_coordinate
         ttanal=ttanal+etime(tim)-tt
        endif
        if (mod(step,nspill).eq.0)
     &      write (31) ((x(l,js),l=1,ndim),js=1,nslices)
!       ifl=0
!       ifr=0
!       do js=1,nslices
!        if (x(1,js) < -1.0d0) ifl=1
!        if (x(1,js) > 1.0d0) ifr=1
!       enddo
!       if (ifl+ifr > 1) then
!           write (31) flag
!           write (31) ((x(l,js),l=1,ndim),js=1,nslices)
!       endif

       enddo         ! over steps

       return
       end
