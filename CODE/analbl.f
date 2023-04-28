      subroutine analbl
      implicit none
      include 'pimc.par'
      integer l,inc,js,k,k1
      include 'pimc.cm'

c update global averages
      call cumul(avp(ietot),av(ietot),anorm(ietot),2)
      call cumul(avp(irg),av(irg),anorm(irg),2)
      call cumul(avp(irho1),av(irho1),anorm(irho1),nhist)
      if (rqmc) then
       avp(ievar)=avp(ievar)-avp(ietot)**2
       avp(ievar+1)=avp(ievar+1)-avp(ietot)**2
       avp(iepot+2)=avp(iepot+2)+avp(ietot)*avp(iepot+4)
       avp(iepot+3)=avp(iepot+3)+avp(ietot)*avp(iepot+5)
       k=iepcor-1
       k1=iepcor1-1
!      print *,'avp(iepcor)'
       do js=1,nslices
        k=k+1
        k1=k1+1
        avp(k)=avp(k)-avp(ietot)*avp(k1)
!       write(*,'(2i6,g15.5,i6,g15.5)' )js,k,avp(k),k1,avp(k1)
       enddo
       call cumul(avp(iepot),av(iepot),anorm(iepot),4)
       call cumul(avp(ieff),av(ieff),anorm(ieff),nslices-2*(nbeg+1))
       call cumul(avp(irho2),av(irho2),anorm(irho2),nhist)
       call cumul(avp(ievar),av(ievar),anorm(ievar),2)
       call cumul(avp(iepcor),av(iepcor),anorm(iepcor),nslices)
      else
       call cumul(avp(iepot),av(iepot),anorm(iepot),2)
       call cumul(avp(ieff),av(ieff),anorm(ieff),nslices/2)
       call cumul(avp(ikin),av(ikin),anorm(ikin),1)
       call cumul(avp(iev),av(iev),anorm(iev),1)
       call cumul(avp(ix0xb2),av(ix0xb2),anorm(ix0xb2),1)
      endif
      if (anormp(idel) > 0.d0) then
       call cumul(avp(idel),av(idel),anorm(idel),1)
       call cumul(avp(iinst),av(iinst),anorm(iinst),1)
       call cumul(avp(idelhist),av(idelhist),anorm(idelhist),nhist)
       call cumul(avp(iinsthist),av(iinsthist),anorm(iinsthist),nhist)
      endif
!     do l=1,nhist
!      call cumul(avp(ipr(l)),av(ipr(l)),anorm(ipr(l)),nhist)
!      call cumul(avp(ipi(l)),av(ipi(l)),anorm(ipi(l)),nhist)
!     enddo
      if(iwigner.eq.1) then
       inc=nwstep/iwaver+1
       do js=1,2
        call cumul(avp(iwkin(js)),av(iwkin(js)),anorm(iwkin(js)),inc)
        call cumul(avp(iwpot(js)),av(iwpot(js)),anorm(iwpot(js)),inc)
        call cumul(avp(iwetot(js)),av(iwetot(js)),anorm(iwetot(js)),inc)
        call cumul(avp(iwtemp(js)),av(iwtemp(js)),anorm(iwtemp(js)),inc)
        call cumul(avp(iwcpp),av(iwcpp),anorm(iwcpp),inc)
       enddo
      endif
      return
      end
