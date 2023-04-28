      subroutine cumul(a,b,c,n)
      implicit none
      integer n,k,i
      real*8 a(n),b(2*n),c(2*n)
      k=0
      do i=1,n
       k=k+2
       b(k)=(b(k)+b(k-1)**2)*c(k)+a(i)**2
       b(k-1)=b(k-1)*c(k-1)+a(i)
       c(k-1)=c(k-1)+1.d0
       c(k)=c(k)+1.d0
      enddo
!DEC$ NOVECTOR
      k=0
      do i=1,n
       k=k+2
!      print *,i,k,c(k-1),c(k)
       b(k-1)=b(k-1)/c(k-1)
       b(k)=b(k)/c(k)-(b(k-1)**2)
      enddo
      return
      end
