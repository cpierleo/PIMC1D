      subroutine wigner
      include 'pimc.par'
      integer ndyn
      parameter (ndyn=1) 
c driver for classical dynamics of wigner function
      integer js1,js2,l,i,iorg,ii
      real*8 gaussn,pcm1,pcm2,psmall,pbig,weight,ronf,err(2,2)
      external gaussn
      include 'pimc.cm'

      call rmar(ronf)
      iorg=dint(nslices*ronf)+1
      do ii=0,ndyn-1
       if (nslices.eq.1) then
        if (ndyn.ne.1) stop 'wigner: classical system with ndyn>1'
        do l=1,ndim
         pcm1=0.0
!        do i=1,nparts
          qw(l,1)=x(l,1)
          qw(l,2)=x(l,1)
          pw(l,1)=sigmapb*gaussn()
!         pcm1=pcm1+pw(l,1)
!        enddo
!        pcm1=pcm1/nparts
!        do i=1,nparts
          pw(l,1)=pw(l,1)  !-pcm1
          pw(l,2)=-pw(l,1)
!        enddo
        enddo
       else
        iorg=iwrap(iorg+ii*nslices/4)
        js1=iorg
        js2=iwrap(iorg+nslices/2)
        do l=1,ndim
         pcm1=0.0
         pcm2=0.0
!        do i=1,nparts
          qw(l,1)=x(l,js1)
          qw(l,2)=x(l,js2)
          pbig=sigmapb*gaussn()
          psmall=sigmaps*gaussn()
          pw(l,1)=pbig+0.5*psmall
          pw(l,2)=pbig-0.5*psmall
!         pcm1=pcm1+pw(l,1)
!         pcm2=pcm2+pw(l,2)
!        enddo
!        pcm1=pcm1/nparts
!        pcm2=pcm2/nparts
!        do i=1,nparts
          pw(l,1)=pw(l,1)   !-pcm1
          pw(l,2)=-(pw(l,2))   !-pcm2)      !backward in time
          do kw=1,2
           do ky=1,2
            call polin2(qgrid,pgrid,phip(1,1,ky,1),nhist,nhist
     .              ,abs(qw(l,kw)),abs(pw(l,kw)),wt(ky,kw),err(ky,kw))
           enddo
           wt(2,kw)=sign(pw(l,kw),wt(2,kw))  ! odd function of pw
          enddo
!        enddo
        enddo
       endif
       call wdynamics
      enddo
      return
      end
