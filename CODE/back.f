      subroutine back(x,y,t,m,mn,k,n)
      implicit none
! to go from one space to the other
!
! in input: x(m,n) input vector
!           t(n,n) transformation matrix
!
! in output: y(m,n) transformed vector
!
      integer n,l,m,mn,k,il,in
      real*8 x(m,mn,k),y(m,mn),t(mn,mn,k),d
      do l=1,m
       do il=1,n
        y(l,il)=0.d0
        do in=1,n
         do ik=1,k
          y(l,il)=y(l,il)+t(il,in,k)*x(l,in,ik)
         enddo
        enddo
       enddo
      enddo
      return
      end
