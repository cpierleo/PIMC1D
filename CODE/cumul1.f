      subroutine cumul1(a,b,c,n)
      implicit none
      integer n,k,i
      real*8 a(n),b(n),c(n)
      do i=1,n
       b(i)=b(i)*c(i)+a(i)
       c(i)=c(i)+1.d0
       b(i)=b(i)/c(i)
      enddo
      return
      end
