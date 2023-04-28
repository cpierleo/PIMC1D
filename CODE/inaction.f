        subroutine inaction
        implicit none
        include 'pimc.par'
        integer il,nb,l,js,jsp,i,is,k,ik,nsl
        integer n,m,js1,is1,it,ismin,mover(mmovers)
        real*8 harmonic,fc,ha,lpsit
        external harmonic
        include 'pimc.cm'
 
        cpot=0.0
        nsl=nslices
        if(rqmc) nsl=nslices-1
c   Initial action
        do js=1,nsl
         js1=iwrap(js+1)
         spring(js)=0.d0
         do l=1,ndim
          spring(js) = spring(js)+(x(l,js)-x(l,js1))**2
         enddo
         spring(js)=spring(js)*cke
        enddo
        if (omega.ne.0.d0) then
         if (rqmc) action=lpsit(x(1,1),model,tau)
     &                   +lpsit(x(1,nslices),model,tau)
         if (nslices>1) then
          do js=1,nslices
           fc=1.d0
           if (rqmc) then
            if(js.eq.1.or.js.eq.nslices) fc=0.5d0
           endif
           cpot=cpot+fc*harmonic(x(1,js),ndim,force(1,js))  !harmonic=tau*V
          enddo
          action=action+cpot
         endif
         cpot=action
         write (2,*) 'initial action = ',action
         write (*,*) 'initial action = ',action
        endif

        return
        end
