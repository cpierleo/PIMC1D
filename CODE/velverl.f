      subroutine velverl(js)
c velocity verlet algorithm
      include 'pimc.par'
      integer js,icall,l,i
      real*8 fact,fact1
      save icall,fact,fact1
      include 'pimc.cm'
      data icall/0/

      if(icall.eq.0) then
       icall=1
       fact=hstep/amass
       fact1=hstep/2.
      endif

      do l=1,ndim
       pw(l,js)=pw(l,js)+fact1*fw(l,js)
       qw(l,js)=qw(l,js)+fact*pw(l,js)
      enddo
      
      call forcew(js)

      do l=1,ndim
       pw(l,js)=pw(l,js)+fact1*fw(l,js)
      enddo
      return
      end
