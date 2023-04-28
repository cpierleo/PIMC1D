       subroutine displace
       implicit none
       include 'pimc.par'
       integer il,l,js,i
       real*8 eps,xold(mdim,mslices),dx(mdim),ronf
       real*8 deltac,cold,cnew,fc,dum(mdim)
       real*8 harmonic,lpsit
       parameter (eps=1.d-32)
       external harmonic,lpsit
       include 'pimc.cm'
 
       if (omega.eq.0.d0) return
c attempt to displace the entire path as a rigid body
       do l=1,ndim
        call rmar(ronf)
        dx(l)=hop*(ronf-0.5d0)
        do js=1,nslices
         xold(l,js)=x(l,js)
         x(l,js)=x(l,js)+dx(l)   !  hop*(ronf-0.5)
        enddo
       enddo
       if(rqmc) then
        cnew=lpsit(x(1,1),model,tau)+lpsit(x(1,nslices),model,tau)
       else
        cnew=0.d0
       endif
       if (nslices>1) then
        do js=1,nslices
         fc=1.d0
         if (rqmc) then 
          if(js.eq.1.or.js.eq.nslices) fc=0.5d0
         endif
         cnew=cnew+fc*harmonic(x(1,js),ndim,dum)
        enddo
       endif
!      write (*,*) cpot,cold,cnew
!       write (*,'(i5,2g20.10)') (js,xold(1,js),x(1,js),js=1,nslices)
!      pause
       deltac=-cnew+cpot   !cold
       ntrydsp=ntrydsp+1
       call rmar(ronf)
       ronf=dlog(max(eps,ronf))
c      write (*,*) delta,ronf
       if ( deltac .ge. ronf ) then
c    ** .. accept **
        naccdsp = naccdsp + 1
c update the action
        cpot=cnew
       else
c    ** .. reject and replace **
        do js=1,nslices
         do l=1,ndim
          x(l,js) = xold(l,js)
         enddo
        enddo
       endif
       action=cpot

       return
       end
