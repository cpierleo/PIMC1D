      subroutine forcew(js)
      include 'pimc.par'
      integer js
c     real*8 ss,ee,r2,vcut
      include 'pimc.cm'
!     dist(a,b,l)=a-b-anint((a-b)*elli(l))*ell(l)
c     ss=2.5
c     ee=10.2
c     r2=(ss/8.5)**2
c     vcut=4.*ee*r2**3*(r2**3-1.)

      totpot(js)=0.0
      do l=1,ndim
       fw(l,js) = -8.d0*potpar*x(l,js)*(x(l,js)**2-1.d0)
       totpot(js) = totpot(js) + potpar*2.d0*(x(l,js)**2-1.d0)**2
      enddo
      return
      end
