      subroutine spill(nb)
      implicit none
      include 'pimc.par'
      integer iflag,l,js,nb,i,nj
      data iflag/0/
      save iflag
      include 'pimc.cm'

      if (iflag.eq.0) then
       iflag=1
       inquire(file=filen(1:ln)//'.rs',exist=ifex)
       if(ifex) then
        rewind(3)
       else
        open(3,file=filen(1:ln)//'.rs',form='unformatted',status='new')
       endif
      else
       rewind(3)
      endif
      write(3) ndim,nslices,((x(l,js),l=1,ndim),js=1,nslices)
      write(3) nb,nstep,nav,(av(i),anorm(i),i=1,nav)
      return
      end
