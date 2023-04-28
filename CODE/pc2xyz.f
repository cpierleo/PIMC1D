      program pc2xyz
      implicit none
      integer mslices,mdim
      parameter (mslices=4800,mdim=3)
      integer nslices,js,ncs,nskip,ndim,icount,iff,nconfs,l,i,ln,ics
      real*8 x(mdim,mslices)
      save icount,iff
      logical accept,ife
      character typeA*1,typeB*1,typeC*1,qid*30,flag*5
      data iff/0/
      if (iff.eq.0) then
       icount=0
       iff=1
       typeA='O'
       typeB='H'
       typeC='C'
      endif

      
      write(*,*) 'file to process (.pc)'
      read(*,*) qid
      ln=index(qid,' ')-1
      inquire(file=qid(1:ln)//'.pc',exist=ife)
      if(ife) then
       open(1,file=qid(1:ln)//'.pc',form='unformatted',status='old')
       open(2,file=qid(1:ln)//'.xyz',form='formatted',status='unknown')
      else
       write(*,*)'crd file does not found'
       stop
      endif
      read(1) nslices,ndim,nconfs
      write(*,*) 'file content:'
      write(*,*) 'nslices=',nslices
      write(*,*) 'ndim=',ndim
      write(*,*) 'nconfs=',nconfs
     
      if (nslices.gt.mslices) stop 'nslices.gt.mslices'
      if (ndim.gt.mdim) stop'ndim.gt.mdim'

      write (*,*) 'input #confs, #nskip'
      read(*,*) ncs,nskip
      do i=1,99999999
!      read(1,end=200)
       read(1,end=200) 
      enddo
200   ncs=min(i-1-nskip,ncs)
      write(*,*) '#confs in qid.crd:',i-1
      write(*,*) '#confs to skip:',nskip
      write(*,*) '#confs to process:',ncs
      rewind(1)
      read(1) 
      do i=1,nskip
!      read(1) 
       read(1) 
      enddo

      flag='#conf'
      do ics=1,ncs
!      read(1) flag
       read(1) ((x(l,js),l=1,ndim),js=1,nslices)
       write(2,*)
       write(2,*)
       write(2,'(a5,i4)') flag,ics
       do js=1,nslices
        write(2,'(i6,g18.8)') js,x(1,js)
       enddo

       print *,' configuration ',ics
       call react(nslices,mdim,x)

      enddo ! ics
      close (1)
      close(2)
      stop
      end
