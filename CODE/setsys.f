        subroutine setsys
        implicit none
        include 'pimc.par'
        integer l,js,i,mslice,k,mstep,mmav,mb,nsl
        integer mmm,mmmm,n,m,nj,mp,js1,it,kp,nndim
        real*8 rronf,ronf
        real*8 rrho,rrs,rksq,rtmin,rtmax,gamma
        logical isex
        real*8 var
        common/twf/var
        include 'pimc.cm'

c    ** Input run parameters **
        write (*,*) 'name of the run'
        read (*,'(a14)') filen
        ln=index(filen,' ')-1
c does this file exist?
        inquire(file=filen(1:ln)//'.in',exist=ifex)
        if(ifex) then
         open(1,file=filen(1:ln)//'.in',status='old',form='formatted')
         open(2,file=filen(1:ln)//'.out',STATUS= 'UNKNOWN',
     &               form='formatted')
         open(10,file=filen(1:ln)//'.eff',STATUS= 'UNKNOWN',
     &         form='formatted')
         open(20,file=filen(1:ln)//'.rho1',STATUS= 'UNKNOWN',
     &         form='formatted')
!         open(21,file=filen(1:ln)//'.phip',STATUS= 'UNKNOWN',
!    &          form='formatted')
        else 
         stop ' file does not exist '
        endif
        read(1,'(a80)') title
        write(*,'(a80)') title
        write(2,'(a80)') title

        read(1,'(a80)') title
        read(1,*) ndim,amass,omega,model,temp,nslices
        write(*,'(a80)') title
        write(*,'(i4,2f15.6,i4,f15.6,i4)') 
     &       ndim,amass,omega,model,temp,nslices
        write(2,'(a80)') title
        write(2,'(i4,2f15.6,i4,f15.6,i4)') 
     &       ndim,amass,omega,model,temp,nslices
        
        if (model /= 1 .and. model /= 2) then
         write (*,*) 'model =',model,' not defined'
         stop
        endif
        
        if(model==2) then
         open(22,file=filen(1:ln)//'.staz',STATUS= 'UNKNOWN',
     &         form='formatted')
         open(23,file=filen(1:ln)//'.inst',STATUS= 'UNKNOWN',
     &         form='formatted')
        endif

        if (mod(nslices,2) /= 0) then
         if (model==1) then
          read(1,'(a80)') title
          read(1,*) var
          write(*,'(a80)') title
          write(*,'(g15.5)') var
          write(2,'(a80)') title
          write(2,'(g15.5)') var
         endif
        endif

        read(1,'(a80)') title
        write(*,'(a80)') title
        write(2,'(a80)') title

        read(1,'(a80)') title
        read(1,*) nblock,nstep,nstart,nmovers,nanal,nbeg
        write(*,'(a80)') title
        write(*,'(6i8)') nblock,nstep,nstart,nmovers,nanal,nbeg
        write(2,'(a80)') title
        write(2,'(6i8)') nblock,nstep,nstart,nmovers,nanal,nbeg

        read(1,'(a80)') title
        read(1,*) nalg,nspill,hop,hop_cm,gammadsp
        write(*,'(a80)') title
        write(*,'(2i8,3f10.5)') nalg,nspill,hop,hop_cm,gammadsp
        write(2,'(a80)') title
        write(2,'(2i8,3f10.5)') nalg,nspill,hop,hop_cm,gammadsp

        read(1,'(a80)') title
        read(1,*) hop_hmc,nstdyn
        write(*,'(a80)') title
        write(*,'(f10.5,i8)') hop_hmc,nstdyn
        write(2,'(a80)') title
        write(2,'(f10.5,i8)') hop_hmc,nstdyn

        read(1,'(a80)') title
        read(1,*) iwigner,hstep,nwstep,iwaver,sigmap,covarp
        write(*,'(a80)') title
        write(*,'(i8,f10.5,2i8,2f10.5)')
     .       iwigner,hstep,nwstep,iwaver,sigmap,covarp
        write(2,'(a80)') title
        write(2,'(i8,f10.5,2i8,2f10.5)')
     .       iwigner,hstep,nwstep,iwaver,sigmap,covarp

        if (nslices.gt.mslices) stop 'nslices.gt.mslices'
        if(mod(nslices,2).eq.0) then
         write(*,*) 'even number of time slices ==> PIMC'
         write(2,*) 'even number of time slices ==> PIMC'
         rqmc=.false.
         nsl=nslices
        else
         write(*,*) 'odd number of time slices ==> RQMC'
         write(2,*) 'odd number of time slices ==> RQMC'
         rqmc=.true.
         nsl=nslices-1
         open(24,file=filen(1:ln)//'.estim',STATUS= 'UNKNOWN',
     &         form='formatted')
        endif
        if (nspill.ge.0.and.nspill.le.nstep) then
         open(31,file=filen(1:ln)//'.pc',form='unformatted',
     &           status='new')
        endif
        pi=4.d0*datan(1.d0)
        spi=dsqrt(pi)
c atomic units: m_electr=1, a_0=hbar2/e2/me=1 (Bohr radius)
c 
        beta = 1.d0 / temp
        tau=beta/max(1,nsl)
c alamb=1/(2 )    
        alamb=0.5d0/amass
        sigma=dsqrt(alamb*tau)
        cke=0.25d0/(alamb*tau)
        potpar=amass*omega**2*tau/2.d0
cpierleo: instanton time analisys
        rdelhist(1)=0
        rdelhist(2)=beta/nhist
        rdelhist(3)=1.d0/rdelhist(2)
        rdelhist(4)=0.d0
        rdelhist(5)=dble(nhist)
        rinsthist(1)=0
        rinsthist(2)=10./nhist
        rinsthist(3)=1.d0/rinsthist(2)
        rinsthist(4)=0.d0
        rinsthist(5)=dble(nhist)
        rmax=6.0d0
        rphist(1)=-rmax
        rphist(2)=2.d0*rmax/nhist
        rphist(3)=1.d0/rphist(2)
        rphist(4)=0.5d0
        rphist(5)=dble(nhist)
        pmax=7.0d0
        pphist(1)=0.d0
        pphist(2)=pmax/nhist
        pphist(3)=1.d0/pphist(2)
        pphist(4)=0.5d0
        pphist(5)=dble(nhist)
        do l=1,nhist
         pimp(l)=pphist(2)*(pphist(4)+l-1)-pphist(1)
        enddo
        write (*,*)
        write (*,*) ' temperature= ',temp
        write (*,*) ' beta*omega = ',beta*omega
        write (*,*) ' tau        = ',tau
        write (*,*) ' lambda*tau = ',sigma**2
        write (2,*) ' temperature= ',temp
        write (2,*) ' beta*omega = ',beta*omega
        write (2,*) ' tau        = ',tau
        write (2,*) ' lambda*tau = ',sigma**2

        if(iwigner.eq.1) then
         sigmapb=(sigmap+covarp)/6.
         sigmaps=(sigmap-covarp)*2./3.
!        write (2,*) 'particle mass: ',amass,'a.u.'
         write (2,*) 'initial squared velocity: ',sigmap
         write (2,*) 'initial covar  velocity: ',covarp
         sigmapb=sqrt(sigmapb)
         sigmaps=sqrt(sigmaps)
        endif

c set the beginning of the random sequence
        inquire(file='rmaseed.dat',exist=ifex)
        if (ifex) then
         call rmaset(1)
        else
         call rmaset(0)
        endif

c    ** Initial position **
        nin=0
        nout=nblock
        inquire(file=filen(1:ln)//'.rs',exist=ifex)
        if(ifex) then
         open(3,file=filen(1:ln)//'.rs',STATUS= 'UNKNOWN',
     &              form='unformatted')
         read(3) nndim,mslice,((x(l,js),l=1,nndim),js=1,mslice)
         write (*,*) 'positions read in qid.rs'
         write (2,*) 'positions read in qid.rs'
         if (nndim.ne.ndim) stop'nndim.ne.ndim'
         if (mslice.ne.nslices) stop'mslice.ne.nslices'
         if (nstart.gt.0) then
          read(3) mb,mstep,mmav,(av(i),anorm(i),i=1,mmav)
          if (mstep.ne.nstep) stop 'mstep.ne.nstep'
          if (mmav.gt.mav) stop 'mmav.ne.mav'
          nin=mb
          nout=nout+mb
         endif
        else
         write (*,*) ' rs file not found '
         write (*,*) 'start from random configuration (sites)'
         write (2,*) ' rs file not found '
         write (2,*) 'start from random configuration (sites)'
         do l=1,ndim
          call rmar(rronf)
          do js=1,nslices
           call rmar(ronf)
           x(l,js) = 0.d0
     &             + (2.0d0*ronf-1.0d0)*sigma/10.d0
          enddo
         enddo
        endif
        if (nspill.ge.0.and.nspill.le.nstep) 
     &  write(31) nslices,ndim,nstep*nblock/nspill

        level=1

        if (nstart.eq.-1) then
         inquire(file=filen(1:ln)//'.pc',exist=ifex)
         if (ifex) then
          open(31,file=filen(1:ln)//'.pc',form='unformatted',
     &         status='old')
          read(31) mslice,mmm,mmmm
          if (mslice.ne.nslices) stop'mslice.ne.nslices in pc file'
          if (mmm.ne.ndim) stop'mmm.ne.ndim in pc file'
          if (mmmm.lt.nstep*nblock) then
           write(*,*) 'number of confs in pc file:',mmmm
           stop 'mmmm.lt.nstep*nblock in pc file'
          endif
          read(1,'(a80)') title
          read(1,'(a80)') title
          read(1,*) level
          write (*,*) 'level has been changed'
          write(*,'(a80)') title
          write(*,'(g15.6,i5)') level
          write (2,*) 'level has been changed'
          write(2,'(a80)') title
          write(2,'(g15.6,i5)') level
          close(1)
         else
          stop 'pc file not found'
         endif
        endif
        return
        end
