       subroutine readpc
       implicit none
       include 'pimc.par'
       integer step,l,js,lp
       integer mover(mmovers)
       real*8 d1,d2,q,p
       real*4 etime,tim(2)
       logical isex
       include 'pimc.cm'
 
       if(iwigner.eq.1) then     ! read the phi(p|q) function
        inquire(file=filen(1:ln)//'.phip',exist=isex)       
        if(isex) then
         open(3,file=filen(1:ln)//'.phip',status='old',form='formatted')
        else
         stop ' file qid.phip does not exist '
        endif
        read(3,*) (rphist(l),l=1,5)
        read(3,*) (pphist(l),l=1,5)
        do l=1,int(rphist(5))
         do lp=1,int(pphist(5))
          read(3,*) q,p,phip(l,lp,1,1),d1,phip(l,lp,2,1),d2
         enddo
        enddo
        do l=1,int(rphist(5))
         qgrid(l)=rphist(2)*(rphist(4)+l-1)-rphist(1)
        enddo
        do l=1,int(pphist(5))
         pgrid(l)=pphist(2)*(pphist(4)+l-1)-pphist(1)
        enddo
        close(3)
       endif
c do one block of MC dynamics
       do step=1,nstep
        read(31) ((x(l,js),l=1,ndim),js=1,nslices)
c       write (*,'(i5,3g20.10)') (js,(x(l,js),l=1,ndim),js=1,nslices)
c compute the action of this path
        call inaction
        write(*,'(i3,g15.6)') step,cpot
c computing averages
        tt=etime(tim)
        call analyze(level,step)
        ttanal=ttanal+etime(tim)-tt
        if(iwigner.eq.1) call wigner
       enddo         ! over steps

       return
       end
