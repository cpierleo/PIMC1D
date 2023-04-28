      subroutine reduce(x,y,t,m,n,nm)
      implicit none
! to go from one space to the other
!
! in input: x(m,n) input vector
!           t(n,n) transformation matrix
!
! in output: y(m,n) transformed vector
!
      integer n,l,m,j,k,nm
      real*8 x(m,nm),y(m,nm),t(nm,nm),d
      do k=1,m
       do l=1,n
        y(k,l)=0.d0
        do j=1,n
         y(k,l)=y(k,l)+t(l,j)*x(k,j)
        enddo
       enddo
      enddo
      return
      end
