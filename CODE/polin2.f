      subroutine polin2(x1a,x2a,ya,m,n,x1,x2,y,dy)
! given arrays x1a(1:m) and x2a(1:n) of independent variables, and an m by n 
! array of function values ya(1:m,1:n), tabulated at the grid points defined 
! by x1a and x2a; and given values x1 and x2 of the independent variables, 
! this routine returnes an interpolated function value y, and an accuracy 
! indication dy (based only on the interpolation in the x1 direction, however).
!
      implicit none
      integer m,n,nmax,mmax
      real*8 dy,x1,x2,y,x1a(m),x2a(n),ya(m,n)
      parameter (nmax=20,mmax=20)
      integer j,k
      real*8 ymtmp(mmax),yntmp(nmax)
      do j=1,m
       do k=1,n
        yntmp(k)=ya(j,k)
       enddo
       call polint(x2a,yntmp,n,x2,ymtmp(j),dy)
      enddo
      call polint(x1a,ymtmp,m,x1,y,dy)
      return
      end

