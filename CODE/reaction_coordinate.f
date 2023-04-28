      subroutine reaction_coordinate
      include 'pimc.par'
      integer  js,ic,iorg(mslices),js1,delay(mslices),iflag
      integer instb(mslices),insta(mslices),ix,io,ip1,im1
      real*8 xa,xb,delhist(nhist),insthist(nhist)
      real*8 z(mslices),dz
      integer instt(nslices),nmax
      include 'pimc.cm'

      xa=-2.d0
      xb=2.d0
      do js=1,nhist
       delhist(js)=0.d0
       insthist(js)=0.d0
      enddo

      do js=1,nslices   ! mappa la posizione in modo da avere i minimi in -0.5 e 0.5
       z(js)=(x(1,js)-xa)/(xb-xa)-0.5d0
      enddo

      ic=0
      js=0
      nmax=nslices
      do while (js < nmax)
       js=js+1
       if (abs(z(js)).lt.0.1) then
        dz=0.d0
        do k=1,150
         dz=max(dz,abs(z(iwrap(js+k))-z(iwrap(js-k))))
         if(dz > 0.8d0) goto 10
        enddo
        goto 11
10      ic=ic+1
        iorg(ic)=js
        instb(ic)=0
        insta(ic)=0
        io=iorg(ic)
!       io1=iwrap(io-1)
        ip1=iwrap(io+1)
        do while (abs(z(ip1)) < 0.5d0)
         instb(ic)=instb(ic)+1
         ip1=iwrap(ip1+1)
        enddo
        im1=iwrap(io-1)
        do while(abs(z(im1)) < 0.5d0)
         insta(ic)=insta(ic)+1
         im1=iwrap(im1-1)
        enddo
        if (z(ip1)*z(im1) < 0.d0) then     ! new instanton found
         instt(ic)=insta(ic)+instb(ic)
         print *,' new instanton found:',ic,js,dz,instt(ic)
         js=js-1+instb(ic)
         if (im1 > iorg(ic)) nmax=im1
!        pause
        else                               ! fake instanton
         ic=ic-1
        endif
11      continue
       endif
      enddo
!     print *,'# of instantons : ',ic
      if (mod(ic,2) .ne. 0) then
       write(*,*) 'WARNING odd numbers of instantons: ',ic
       write(2,*) 'WARNING odd numbers of instantons: ',ic
      endif
      do js=1,ic
       js1=js+1
       if(js1>ic) js1=1
       delay(js)=iwrap( iwrap(iorg(js1)-insta(js1))
     &                 -iwrap(iorg(js)+instb(js)) )
       ix=min(nhist,max(1,int((tau*delay(js)-rdelhist(1))*rdelhist(3)
     &                    +rdelhist(4))+1))
       delhist(ix)=delhist(ix)+1.d0
       ix=min(nhist,max(1,int((tau*instt(js)-rinsthist(1))*rinsthist(3)
     &                    +rinsthist(4))+1))
       insthist(ix)=insthist(ix)+1.d0
      enddo
      do js=1,ic
       call cumul1(dble(delay(js)),avp(idel),anormp(idel),1)
       call cumul1(dble(instt(js)),avp(iinst),anormp(iinst),1)
      enddo
      if (ic>0) then
       call cumul1(delhist(1),avp(idelhist),anormp(idelhist),nhist)
       call cumul1(insthist(1),avp(iinsthist),anormp(iinsthist),nhist)
      endif

      return
      end
