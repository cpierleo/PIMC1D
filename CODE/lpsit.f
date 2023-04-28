      function lpsit(x,model,tau)
      implicit none
      integer model
      real*8 lpsit,x,tau
      real*8 c0,c2,c4,c6
      data c0,c2,c4,c6/1.46486,-0.26217,0.0386943,0d0/
      real*8 var
      common/twf/var

!     var=1.0d0
      if (model==1) then   ! lpsit=-log(psi_t)
       lpsit=var*x**2/2.d0
      else
       lpsit=c0+x**2*(c2+x**2*(c4+x**2*c6))
      endif
      return
      end

      function eloc(x,model,tau,alamb)
      implicit none
c  eloc is H\psi/\psi+v=-lambda(u"-u'^2)+v
c  u1=u'
c  u2=u"
      integer model
      real*8 x,eloc,y,tau,dum,u1,u2,alamb
      real*8 harmonic
      real*8 c0,c2,c4,c6
      data c0,c2,c4,c6/1.46486,-0.26217,0.0386943,0d0/
      real*8 var
      common/twf/var

      if(model==1) then
       u1=var*x
       u2=var
      else
       u1=2.d0*x*(c2+x**2*(2.d0*c4+x**2*3.d0*c6))
       u2=2.d0*(c2+x**2*(6.d0*c4+x**2*15.d0*c6))
      endif
      eloc = alamb*(u2-u1**2) + harmonic(x,1,dum)/tau
  
      return
      end
