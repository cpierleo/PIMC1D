      function sang(ndim)
c solid angle in ndim dimensions
      real*8 sang,pi
      integer ndim
      pi=3.1415 92653 58979d0
      sang=2.d0*pi*dfloat(ndim-1)
      if(ndim.eq.1)sang=1.d0
      return
      end

