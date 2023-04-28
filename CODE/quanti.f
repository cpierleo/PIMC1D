      subroutine quanti(ist)
      include 'pimc.par'
      integer ist,js,it,jt,k,kp,l,i,ik
      real*8 pcmw(ndim,2),ecm(2)
      include 'pimc.cm'

c energies
      do js=1,2
       ekin(ist,js)=0.0
       epot(ist,js)=0.0
       ecm(js)=0.0
       do l=1,ndim
!       pcmw(l,js)=0.0
!       do i=1,nparts
!        pcmw(l,js)=pcmw(l,js)+pw(l,i,js)
         ekin(ist,js)=ekin(ist,js)+pw(l,js)*pw(l,js)
!       enddo
!       ecm(js)=ecm(js)+pcmw(l,js)**2
       enddo
       ekin(ist,js)=ekin(ist,js)/(2.*amass*nparts)
       wtemp(ist,js)=ekin(ist,js)/1.5
!      do it=1,ntypes
!       do jt=1,it
         epot(ist,js)=epot(ist,js)+totpot(js)
!       enddo
!      enddo
       epot(ist,js)=epot(ist,js)/nparts
       et(ist,js)=epot(ist,js)+ekin(ist,js)
!      press(ist,js)=(nparts/beta+virtot(js)/3.)/vol
      enddo
      write(69,'(i6,4g20.12,4g15.8)') ist*iaver,et(ist,1),et(ist,2)
     &    ,wtemp(ist,1),wtemp(ist,2),epot(ist,1),epot(ist,2)
c symmetrized intermediate scattering function
c <\rho_k(q_1)\rho_{-k}(q_2)>
!     do js=1,2
!      call rhokall(qw(1,1,js),rhokw(1,1,js),pww(1,1,js),nvects)
!     enddo
!     skcntw=skcntw+wt
!     do k=1,nlamb+jaddk
!      sktw(ist,k)=0.0
!      do kp=2*kmult(k-1)+1,2*kmult(k)   !,2
!       sktw(ist,k)=sktw(ist,k)+wt*rhokw(kp,1,1)*rhokw(kp,1,2) 
c       sktw(ist,k)=sktw(ist,k)+wt*(rhokw(kp,1,1)*rhokw(kp,1,1) 
c    &                             +rhokw(kp,1,2)*rhokw(kp,1,2))/2.
!      enddo
!      sktw(ist,k)=sktw(ist,k)/(kmult(k)-kmult(k-1))/nparts
!     enddo
c symmetrized velocity-velocity c.f.
      cpp(ist)=0.0
      do i=1,nparts
       do l=1,ndim
        cpp(ist)=cpp(ist)+pw(l,1)*pw(l,2)
       enddo
      enddo
      cpp(ist)=cpp(ist)/nparts/ndim

      return
      end
