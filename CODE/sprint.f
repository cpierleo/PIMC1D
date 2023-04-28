      subroutine sprint(nb)
      implicit none
      include 'pimc.par'
      integer k,l,nb,nsl,intv,kp,ki,lp,lk,nend,k1
      real*8 a1,a2
      real*8 cnorm,q
      include 'pimc.cm'
c output file (2)
      write (2,1011) int(anorm(1))
1011  format ('############## block',i5,' completed ################')
c     write (2,*)
      cnorm=max(anorm(1)-1.d0,1.d0)

      if (rqmc) then

       a1=dsqrt(max(av(ietot+1),0.d0)/cnorm)
       write(2,'(3g20.10,18H total   energy   )')avp(ietot),av(ietot),a1
       a1=dsqrt(max(av(ievar+1),0.d0)/cnorm)
       write(2,'(3g20.10,18H variance(pure es))')
     &                                 avp(ievar),av(ievar),a1
       a1=dsqrt(max(av(ievar+3),0.d0)/cnorm)
       write(2,'(3g20.10,18H variance(mix  es))')
     &                                 avp(ievar+1),av(ievar+2),a1
!      print *,'sprint 1',av(iepot+1)
       a1=dsqrt(max(av(iepot+1),0.d0)/cnorm)
       write(2,'(3g20.10,18H poten ener (pure))')
     &                                 avp(iepot),av(iepot),a1
!      print *,'sprint 3',av(iepot+3)
       a1=dsqrt(max(av(iepot+3),0.d0)/cnorm)
       write(2,'(3g20.10,18H poten ener (mixt))')
     &                                 avp(iepot+1),av(iepot+2),a1
!      print *,'sprint 5',av(iepot+5)
       a1=dsqrt(max(av(iepot+5),0.d0)/cnorm)
       write(2,'(3g20.10,18H pot en imp (pure))')
     &                                 avp(iepot+2),av(iepot+4),a1
!      print *,'sprint 7',av(iepot+7)
       a1=dsqrt(max(av(iepot+7),0.d0)/cnorm)
       write(2,'(3g20.10,18H pot en imp (mixt))')
     &                                 avp(iepot+3),av(iepot+6),a1
!      a1=dsqrt(max(av(ikin+1),0.d0)/cnorm)
!      write(2,'(18Hkinet ener (pure):,3g20.10)')
!    &                                 avp(ikin),av(ikin),a1
!      a1=dsqrt(max(av(ikin+3),0.d0)/cnorm)
!      write(2,'(18Hkinet ener (mixt):,3g20.10)')
!    &                                 avp(ikin+1),av(ikin+2),a1
!      a1=dsqrt(max(av(ikin+5),0.d0)/cnorm)
!      write(2,'(18Hkin en imp (pure):,3g20.10)')
!    &                                 avp(ikin+2),av(ikin+4),a1
!      a1=dsqrt(max(av(ikin+7),0.d0)/cnorm)
!      write(2,'(18Hkin en imp (mixt):,3g20.10)')
!    &                                 avp(ikin+3),av(ikin+6),a1

      else

       a1=dsqrt(max(av(ikin+1),0.d0)/cnorm)
       write(2,'(3g20.10,18H kinetic energy   )')avp(ikin),av(ikin),a1
       if (omega.ne.0.0) then
        a1=dsqrt(max(av(iev+1)/cnorm,0.d0))
        write(2,'(3g20.10,18H virial kin. en.  )')avp(iev),av(iev),a1
       endif
       a1=dsqrt(max(av(iepot+1),0.d0)/cnorm)
       write(2,'(3g20.10,18H potential energy )')avp(iepot),av(iepot),a1
       a1=dsqrt(max(av(ietot+1),0.d0)/cnorm)
       write(2,'(3g20.10,18H total   energy   )')avp(ietot),av(ietot),a1
       a1=dsqrt(max(av(ietot+3),0.d0)/cnorm)
       write(2,'(3g20.10,18H vir tot energy   )')
     &                                 avp(ietot+1),av(ietot+2),a1
       a1=dsqrt(max(av(ix0xb2+1),0.d0)/cnorm)
       write(2,'(3g20.10,18H two points cor f )')
     &                                 avp(ix0xb2),av(ix0xb2),a1

      endif

      if (nslices>1) then
       a1=dsqrt(max(av(irg+1)/cnorm,0.d0))
       if(av(irg).ne.0.d0)
     &   write(2,'(3g20.10,18H giration radius  )')dsqrt(avp(irg)),
     &              dsqrt(av(irg)),a1/2/dsqrt(abs(av(irg)))
      endif
      if(model==2) then
       a1=dsqrt(max(av(idel)/cnorm,0.d0))
       write(2,'(3g20.10,18H steady state time )')
     &                                 avp(idel),av(idel),a1
       a1=dsqrt(max(av(iinst)/cnorm,0.d0))
       write(2,'(3g20.10,18H transient time    )')
     &                                 avp(iinst),av(iinst),a1
      endif
      mratio=1.0
      if (ntry.gt.0) mratio = real(nacc) / real(ntry)
      write(2,'(1g20.10,18H Acceptance levy  )') mratio
      if (gammadsp.gt.0.0) then
       if(ntrydsp.ne.0) mratio = real(naccdsp) / real(ntrydsp)
       write(2,'(1g20.10,18H Acceptance displ )') mratio
      endif
      call flush(2)

c effective mass file (10)
      k=ieff-1
      intv=2**(level-1)
      nsl=nslices/intv
      if(rqmc) then
       nend=nsl-2*(nbeg+1)
      else
       nend=nsl/2
      endif
      do l=1,nend
       k=k+2
       a1=dsqrt(max(av(k)/cnorm,0.d0))
       write(10,'(3g20.10)') l*intv*tau,av(k-1),a1
      enddo
      call flush(10)
      rewind(10)

c rho(x)
      k=irho1-1
      k1=irho2-1
      do l=1,nhist
       k=k+2
       a1=dsqrt(max(av(k)/cnorm,0.d0))
       if (.not.rqmc) then
        write(20,'(3g20.10)') 
     &      rphist(2)*(l-1)+rphist(1),av(k-1),a1
       else
        k1=k1+2
        a1=dsqrt(max(av(k1)/cnorm,0.d0))
        write(20,'(5g20.10)') 
     &      rphist(2)*(l-1)+rphist(1),av(k-1),a1,av(k1-1),a1
       endif
      enddo
      call flush(20)
      rewind(20)

      if(rqmc) then
c new estimators
       k=iepcor-1
!      k1=iekcor-1
       do l=1,nslices
        k=k+2
!       k1=k1+2
        a1=dsqrt(max(av(k)/cnorm,0.d0))
!       a2=dsqrt(max(av(k1)/cnorm,0.d0))
        write(24,'(5g20.10)') (l-1)*tau,av(k-1),a1!,av(k1-1),a2
       enddo
       call flush(24)
       rewind(24)
      endif

      if(model==2) then
c plateau time
       k=idelhist-1
       do l=1,nhist
        k=k+2
        a1=dsqrt(max(av(k)/cnorm,0.d0))
        write(22,'(3g20.10)')
     &       rdelhist(2)*(l-1)+rdelhist(1),av(k-1),a1
       enddo
       call flush(22)
       rewind(22)

c transient time
       k=iinsthist-1
       do l=1,nhist
        k=k+2
        a1=dsqrt(max(av(k)/cnorm,0.d0))
        write(23,'(3g20.10)')
     &       rinsthist(2)*(l-1)+rinsthist(1),av(k-1),a1
       enddo
       call flush(23)
       rewind(23)
      endif


c phi(p|q)
!     write (21,'(5g15.7)') (rphist(l),l=1,5)
!     write (21,'(5g15.7)') (pphist(l),l=1,5)
!     lk=0
!     do l=1,nhist
!      q=rphist(2)*(rphist(4)+l-1)-rphist(1)
!      k=ipr(l)-1
!      ki=ipi(l)-1
!      if (mod(l,20).eq.0) lk=lk+1
!      do lp=1,nhist
!       k=k+2
!       ki=ki+2
!       a1=dsqrt(max(av(k)/cnorm,0.d0))
!       a2=dsqrt(max(av(ki)/cnorm,0.d0))
!       write(21,'(6g15.7)') q,pimp(lp),av(k-1),a1,av(ki-1),a2
!       if(mod(l,20).eq.0) then
!        write(21+lk,'(6g20.10)') q,pimp(lp),av(k-1),a1,av(ki-1),a2
!       endif
!      enddo
!      call flush(21+lk)
!      rewind(21+lk)
!     enddo
!     call flush(21)
!     rewind(21)

      return
      end
