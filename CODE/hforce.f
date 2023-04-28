       subroutine hforce
       implicit none
c potential and forces
       include 'pimc.par'
       integer il,l,js,js1,jsm1
       include 'pimc.cm'

       do l=1,ndim
        do js=1,nslices
         force(l,js)=0.d0
        enddo
       enddo
       vcl=0.d0
       do l=1,ndim
        do js=1,nslices
         js1=iwrap(js+1)
         jsm1=iwrap(js-1)
         vcl=vcl+potpar*2.d0*(x(l,js)**2-1.d0)**2
     .       +cke*(x(l,js)-x(l,js1))**2  ! internal spring term
!        force(l,js)=(-2.d0*potpar*x(l,js)
         force(l,js)=(-8.d0*potpar*x(l,js)*(x(l,js)**2-1.d0)
     .                -2.d0*cke*(2.d0*x(l,js)-x(l,js1)-x(l,jsm1)))/beta
        enddo
       enddo
       vcl=vcl/beta
       return
       end
