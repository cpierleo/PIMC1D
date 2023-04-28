        function gaussn()
c       real *8 r1, r2, dlog, dcos, dsqrt
        real *8 r1, r2
        real*8 gaussn
10      continue
        call rmar(r1)
        call rmar(r2)
         if(abs(r1).lt.0.0000000000001d0) goto 10
        gaussn= dsqrt(-2.0d0*dlog(r1))*dcos(6.2831853d0*r2)
c       print *,r1,r2,gaussn
         if(dabs(gaussn).gt.100.d0) stop
        return
        end

      subroutine rmaset(iunit)
c     initializing routine for rmar, must be called before 
c     beginning to generate pseudorandom numbers by means of the
c     subroutine rmar
c     ranges : 0<=ij<=31328 and 0<=kl<=30081.
c     implicit real*8 (a-h,o-z)
      implicit none
      real*8 u,c,cd,cm,s,t
      integer i,j,k,m,iunit,n,kl,ij,ii,jj

      common/raset1/u(97),c,cd,cm,i,j
      data ij,kl/1802,9373/
c
      if(iunit.eq.0)                            then
       i=mod(ij/177, 177)+2
       j=mod(ij,     177)+2
       k=mod(kl/169, 178)+1
       m=mod(kl,     169)
       print '(a,2i7,4i4)',
     &       ' marsaglia initialized (seed=0): ',ij,kl,i,j,k,m
c
       do ii=1,97
        s=0.0d00
        t=0.5d00
        do jj=1,24
         n=mod(mod(i*j,179)*k, 179)
         i=j
         j=k
         k=n
         m=mod(53*m+1, 169)
         if(mod(m*n,64).ge.32) s=s+t
         t=0.5d00*t
        enddo
        u(ii)=s
       enddo
c
       c =  362436.d00/16777216.d00
       cd= 7654321.d00/16777216.d00
       cm=16777213.d00/16777216.d00
c
       i=97
       j=33
c
      else
      
       print*,'marsaglia initialized, seed=',iunit
       open(41,file='rmaseed.dat',status='old',form='unformatted')
       rewind 41
       read(41) u,c,cd,cm,i,j
       close(41)
      end if
c
      return
      end


      subroutine rmar(xr)          
c     pseudo random number generator 
c     proposed by marsaglia, zaman and tsang
c     implicit real*8 (a-h,o-z)
      implicit none
      real*8 xr,u,c,cd,cm,uni
      integer i,j
      common/raset1/u(97),c,cd,cm,i,j
c
      uni=u(i)-u(j)
      if(uni.lt. 0.0e00) uni=uni+1.0e00
      u(i)=uni           
      i=i-1              
      if(i.eq.0) i=97    
      j=j-1              
      if(j.eq.0) j=97   
      c=c-cd            
      if(c  .lt. 0.0e00)   c=c+cm
      uni=uni-c
      if(uni.lt. 0) uni=uni+1
      xr=uni
c
      return
      end
      
         

      subroutine  rmaget
c     implicit real*8 (a-h,o-z)
      implicit none
      real*8 u,c,cd,cm
      integer i,j
      common/raset1/u(97),c,cd,cm,i,j
      open(41,file='rmaseed.dat',status='unknown',form='unformatted')
      write(41) u,c,cd,cm,i,j
      return
      end
