        subroutine normal_mode
        implicit none
        include 'pimc.par'
        integer nb,is,js,k,ms,ls,is1
        real*8 fc,fci,fp,ff
!        real*8 xtoq(mslices,mslices)
!        real*8 qtox(mslices,mslices)
        real*8 One(mslices,mslices)
        include 'pimc.cm'
c setting up the normal modes transformations:
c
c    qr_i=\sum_{j=1}^nslices xtoqr(i,j) x_j     real part
c    qi_i=\sum_{j=1}^nslices xtoqi(i,j) x_j     imaginary part
c
c    x_i=\sum_j=1^nslices [xtoqr(i,j) qr_j + xtoqi(i,j) qi_j]
c
        fc=1.d0/dsqrt(dble(nslices))
        fp=pi/nslices
!       write(*,*) 'XTOQ'
        do js=1,nslices
         do is=1,nslices-2
          if (mod(is,2)==0) then
           xtoq(is,js)=dcos(fp*js*is)*fc
          else
           xtoq(is,js)=dsin(fp*js*(is+1))*fc
          endif
         enddo
         xtoq(nslices-1,js)=(-1.d0)**js*fc
         xtoq(nslices,js)=fc
!        write(*,'(10f15.8)') (xtoq(is,js),is=1,nslices)
        enddo

!       write(*,*) 'QTOX'
        do js=1,nslices
         do is=1,nslices-2
          if (mod(is,2) == 0) then
           qtox(js,is)=2.d0*dcos(fp*js*is)*fc
          else
           qtox(js,is)=2.d0*dsin(fp*js*(is+1))*fc
          endif
         enddo
         qtox(js,nslices-1)=(-1.d0)**js*fc
         qtox(js,nslices)=fc
!        write(*,'(10f15.8)') (qtox(is,js),is=1,nslices)
        enddo
!       do js=1,nslices
!        do is=1,nslices
!         One(js,is)=0.d0
!         do ls=1,nslices
!          One(js,is)=One(js,is)+xtoq(js,ls)*qtox(ls,is)
!         enddo
!        enddo
!       enddo
!       write(*,*) 'One'
!       do js=1,nslices
!        write(*,'(10f15.8)') (One(is,js),is=1,nslices)
!       enddo
c eigenvalues
        do js=1,nslices-1
         if (mod(js,2) == 0) then
          sigma2q(js)= alamb*tau/(1.d0-dcos(fp*js))/2.d0
          sigma2q(js-1)= sigma2q(js)
          sigmaq(js)=dsqrt( sigma2q(js) )
          sigmaq(js-1)=sigmaq(js)
          print *,js-1,sigma2q(js-1),sigmaq(js-1)
          print *,js,sigma2q(js),sigmaq(js)
         endif
        enddo
        sigma2q(nslices-1) = alamb*tau/2.d0
        sigmaq(nslices-1)=dsqrt( sigma2q(nslices-1) )
        print *,nslices-1,sigma2q(nslices-1),sigmaq(nslices-1)
        return
        end

