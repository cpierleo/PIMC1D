      subroutine analyze(lvl,step)
      implicit none
      include 'pimc.par'
      integer k,j,l,js,is,ik,ncount(mslices),lvl,im,nl
      integer intv,nsl,nt,jm,js1,it
      integer ip,kp,step,ix,jsp1,jsm1,jmid,jin
      real*8 rg(2),etot(2),ekine,epote(mslices),cm(mslices)
      real*8 evir,dd2,evirc,ahist(nhist,2),phist(nhist,nhist,2)
      real*8 harmvir,x12,ex12,ppp,el1,el2,eloc,xcm(mdim),evar(2)
      real*8 epcor(mslices),epmix,harmonic,x0xb2
      real*8 epint(2,2),eploc(6),dum(mdim)
      external harmvir,eloc,harmonic
      include 'pimc.cm'
c
c  compute block average
c
      intv=1
      nsl=nslices  !/intv
      if (rqmc) nsl=nslices-1
      jmid=nsl/2+1
!     do nt=1,intv
c
c  observables
c
       etot(1)=0.0
       etot(2)=0.0
       rg(1)=0.d0
       rg(2)=0.d0
       evirc=0.d0
       ekine=0.0
       x0xb2=0.d0
       do k=1,nsl
        epote(k)=0.0
        cm(k)=0.d0
        ncount(k)=0.d0
       enddo
       do k=1,nhist
        ahist(k,1)=0.d0
        if (rqmc) ahist(k,2)=0.d0
       enddo
!      do js=1,2
!       do k=1,nhist
!        do l=1,nhist
!         phist(k,l,js)=0.d0
!        enddo
!       enddo
!      enddo

c computes the centroid position
       do l=1,ndim
        xcm(l)=0.d0
        do js=1,nslices
         xcm(l)=xcm(l)+x(l,js)
        enddo
        xcm(l)=xcm(l)/dble(nslices)
       enddo
       
       if (rqmc) then
        el1=eloc(x(1,1),model,tau,alamb)
        el2=eloc(x(1,nslices),model,tau,alamb)
        etot(1)=0.5d0*(el1+el2)
        evar(1)=el1*el2
        evar(2)=0.5d0*(el1**2+el2**2)
        do js=1,nslices
         epote(js)=harmonic(x(1,js),ndim,dum)/tau
        enddo
        epmix=0.5d0*(epote(1)+epote(nslices))
        do js=1,nslices
         epcor(js)=0.5d0*(el1*epote(js)+el2*epote(nslices-js+1))
        enddo
        do k=1,2   ! k=1 int_0^beta; k=2 int_beta/2^beta
         epint(k,1)=0.d0
         epint(k,2)=0.d0
         jin=1+(k-1)*nslices/2
         do js=jin+1,nslices-1
          epint(k,1)=epint(k,1)+epcor(js)
          epint(k,2)=epint(k,2)+epote(js)
         enddo
         epint(k,1)=tau*(epint(k,1)+0.5d0*(epcor(jin)+epcor(nslices)))
         epint(k,2)=tau*(epint(k,2)+0.5d0*(epote(jin)+epote(nslices)))
        enddo
        eploc(1)=epote(jmid)
        eploc(2)=epmix
        eploc(3)=epote(jmid)-2.d0*epint(2,1)
        eploc(4)=epmix-epint(1,1)
        eploc(5)=2.d0*epint(2,2)
        eploc(6)=epint(1,2)
       else    ! pimc
        do js=1,nsl
         js1=iwrap(js+1)
         dd2=0.d0
         do l=1,ndim
          dd2=dd2+(x(l,js1)-x(l,js))**2
         enddo
         ekine=ekine+dd2
         if (omega.ne.0.0d0) 
     &     evirc=evirc+harmvir(xcm(1),x(1,js),ndim)
        enddo  ! over js
        ekine=(dble(ndim)/2.d0-cke*ekine/nsl)/tau  !(tau*intv)
        epote(1)=cpot/beta
        evir=0.5d0*(dble(ndim)+evirc)/beta
        etot(1)=ekine+epote(1)
        etot(2)=evir +epote(1)
!       etot(1)=etot(1)-dble(ndim)/2.d0/beta
!       etot(2)=etot(2)-dble(ndim)/2.d0/beta
       endif
       do js=1,nslices    !nt,nslices-1+nt,intv
        js1=iwrap(js+nslices/2)
        do l=1,ndim
         rg(1)=rg(1)+(x(l,js)-xcm(l))**2
         x0xb2=x0xb2+x(l,js)*x(l,js1)
        enddo
       enddo
       rg(1)=rg(1)/dfloat(nsl)
       rg(2)=xcm(1)
       x0xb2=x0xb2/float(nslices)
c imaginary time diffusion
       if(rqmc) then
        if (nslices>1) then
         do js=nbeg+1,nsl-nbeg-1
          do ik=1,js-nbeg-2
           is=js-ik
           do l=1,ndim
            cm(ik)=cm(ik)+(x(l,js)-x(l,is))**2
           enddo
           ncount(ik)=ncount(ik)+1
          enddo
         enddo
         do k=1,nsl-2*(nbeg+1)
          cm(k)=cm(k)  !/(2.d0*ndim*alamb*tau*intv*k)
     &                /dfloat(max(ncount(k),1))
         enddo
        endif
       else   !pimc
        js=1  !nt
        do jm=1,nsl-1
         do ik=1,jm-1
          is=iwrap(js-ik*intv)
          k=ik
          if (k.gt.nsl/2) k=nsl/2-(k-nsl/2)
c          write (*,*) js,is,k
          do l=1,ndim
           cm(k)=cm(k)+(x(l,js)-x(l,is))**2
          enddo
          ncount(k)=ncount(k)+1
         enddo
         js=js+intv
        enddo 
        do k=1,nsl/2
         cm(k)=cm(k)  !/(2.d0*ndim*alamb*tau*intv*k)
     &               /dfloat(max(ncount(k),1))
        enddo
       endif
c rho(r)
!      do l=1,ndim
       if (rqmc) then
         js=nslices/2+1
         ppp=x(1,js)
         ix=min(nhist,max(1,int((ppp-rphist(1))*rphist(3)+rphist(4))+1))
         ahist(ix,1)=ahist(ix,1)+1.d0
         ppp=x(1,1)
         ix=min(nhist,max(1,int((ppp-rphist(1))*rphist(3)+rphist(4))+1))
         ahist(ix,2)=ahist(ix,2)+0.5d0
         ppp=x(1,nslices)
         ix=min(nhist,max(1,int((ppp-rphist(1))*rphist(3)+rphist(4))+1))
         ahist(ix,2)=ahist(ix,2)+0.5d0
       else
        do js=1,nslices
         jsp1=iwrap(js+1)
         jsm1=iwrap(js-1)
         ppp=x(1,js)    !(x(1,1)-x(1,nslices))
         ix=min(nhist,max(1,int((ppp-rphist(1))*rphist(3)+rphist(4))+1))
         ahist(ix,1)=ahist(ix,1)+1.d0/nslices
!        x12=0.5d0*cke*(x(1,jsp1)-x(1,jsm1))**2
!        ex12=exp(x12)
!        do k=1,nhist
!         ppp=exp(-2.d0*alamb*tau*pimp(k)**2)
!         phist(k,ix,1)=ppp*ex12*cos(x12*pimp(k))
!         phist(k,ix,2)=ppp*ex12*sin(x12*pimp(k))
!        enddo
        enddo
       endif
!      enddo
       
         
c update local averages
       call cumul1(etot(1),avp(ietot),anormp(ietot),2)
       call cumul1(rg(1),avp(irg),anormp(irg),2)
       if(rqmc) then
        call cumul1(evar(1),avp(ievar),anormp(ievar),2)
        call cumul1(cm(1),avp(ieff),anormp(ieff),nsl-2*(nbeg+1))
        call cumul1(ahist(1,1),avp(irho1),anormp(irho1),nhist)
        call cumul1(ahist(1,2),avp(irho2),anormp(irho2),nhist)
        call cumul1(eploc(1),avp(iepot),anormp(iepot),6)
        call cumul1(epcor(1),avp(iepcor),anormp(iepcor),nslices)
        call cumul1(epote(1),avp(iepcor1),anormp(iepcor1),nslices)
       else
        call cumul1(ekine,avp(ikin),anormp(ikin),1)
        call cumul1(cm(1),avp(ieff),anormp(ieff),nsl/2)
        call cumul1(evir,avp(iev),anormp(iev),1)
        call cumul1(epote(1),avp(iepot),anormp(iepot),2)
        call cumul1(ahist(1,1),avp(irho1),anormp(irho1),nhist)
        call cumul1(x0xb2,avp(ix0xb2),anormp(ix0xb2),1)
       endif
!      do l=1,nhist
!       call cumul1(phist(1,l,1),avp(ipr(l)),anormp(ipr(l)),nhist)
!       call cumul1(phist(1,l,2),avp(ipi(l)),anormp(ipi(l)),nhist)
!      enddo
 
!     enddo     ! loop over nt

      return
      end
