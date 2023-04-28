      subroutine wsprint(inc,ntraj)
      include 'pimc.par'
      integer inc,i,j,k,k1,k2,k3,k4,m,ncc,ntraj,kk
      real*8 an,a0
      character name*30
      include 'pimc.cm'

      an=real(max(1,ntraj-1))
!     ln=index(filen,' ')-1
      name=filen(1:ln)//'.wther'
      open(71,file=name,form='formatted',status='unknown')
      do i=0,inc-1
       j=2*i
       k=j+1
       write(71,'(f10.5,8g15.6)') 2.*i*iaver*hstep
     &      ,(av(iwkin(1)+j)+av(iwkin(2)+j))/2.
     &      ,sqrt((av(iwkin(1)+k)+av(iwkin(2)+k))/2./an)
     &      ,(av(iwpot(1)+j)+av(iwpot(2)+j))/2.
     &      ,sqrt((av(iwpot(1)+k)+av(iwpot(2)+k))/2./an)
     &      ,(av(iwetot(1)+j)+av(iwetot(2)+j))/2.
     &      ,sqrt((av(iwetot(1)+k)+av(iwetot(2)+k))/2./an)
     &      ,(av(iwtemp(1)+j)+av(iwtemp(2)+j))/2.
     &      ,sqrt((av(iwtemp(1)+k)+av(iwtemp(2)+k))/2./an)
      enddo
      close(71)
!     name=filen(1:ln)//'.w0'
!     open(71,file=name,form='formatted',status='unknown')
!     j=-2
!     do k=1,nlamb+jaddk
!      do kk=kmult(k-1),kmult(k)
!      j=j+2
!      m=j+1
!      write (71,'(f10.5,10g15.6)') 
!    &       rknorm(k),av(ifk0+j),sqrt(av(ifk0+m))/an
!      enddo
!     enddo
!     close(71)
!     ncc=0
!     do k=1,nlamb+jaddk,5
!      ncc=ncc+1
!      name=filen(1:ln)//'.w'
!      call numapp(name,ncc)
!      open(71,file=name,form='formatted',status='unknown')
!      k1=k+1
!      k2=k+2
!      k3=k+3
!      k4=k+4
!      write(71,111) rknorm(k),rknorm(k1),rknorm(k2)
!    &               ,rknorm(k3),rknorm(k4)
111    format(1h#,5f15.7)
!      do i=0,inc-1
!       j=2*i
!       m=j+1
!       write(71,'(f10.5,10g15.6)') 2.*i*iaver*hstep/0.91183189
!    &            ,av(ifktw(k)+j),sqrt(av(ifktw(k)+m)/an)
!    &            ,av(ifktw(k1)+j),sqrt(av(ifktw(k1)+m)/an)
!    &            ,av(ifktw(k2)+j),sqrt(av(ifktw(k2)+m)/an)
!    &            ,av(ifktw(k3)+j),sqrt(av(ifktw(k3)+m)/an)
!    &            ,av(ifktw(k4)+j),sqrt(av(ifktw(k4)+m)/an)
!      enddo
!      close(71)
!     enddo
! velocity-velocity tcf
      name=filen(1:ln)//'.cvv'
      open(71,file=name,form='formatted',status='unknown')
      do i=0,inc-1
       j=2*i
       m=j+1
       write(71,'(f10.5,2g15.6)') 2.*i*iaver*hstep
     &            ,av(iwcpp+j),sqrt(av(iwcpp+m)/an)
      enddo
      close(71)

      return 
      end
