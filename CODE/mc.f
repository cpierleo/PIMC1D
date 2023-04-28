       subroutine mc
       implicit none
       include 'pimc.par'
       integer il,l,js,jsm,i,k,ist,ip1,im1,js1
       integer lm,jm,jk,step,iflag,jl,lvl,jdn,jup,jmid,kdiv
       integer mover(mmovers),itspan,itar,iorg,kspan,ndiv
       real*8 ht,vcm(mdim),veloc(mdim,mslices),ek,ekn,eold,enew
       real*8 fspring(mdim,mslices,2),fp(mdim,mslices)
       real*8 fpold(mdim,mslices),ft2(mslices,2),ftot(2)
       real*8 eps,sg,xold(mdim,mmovers),xm,x1(mdim,mslices)
       real*8 deltac,cold,cnew,dd2(2),dx(mdim),weight
       real*8 ronf
       real*8 gaussn,harmonic
       real*4 etime,tim(2)
       parameter (eps=1.d-32)
       external gaussn,harmonic
       character*80 messagge
       include 'pimc.cm'
       data iflag/0/
 
c do one block of MC

       do step=1,nstep

        dd2(1)=0.d0
        do i=1,nslices
         ip1=iwrap(i+1)
         do l=1,ndim
          dd2(1)=dd2(1)+(x(l,i)-x(l,ip1))**2
          xold(l,i)=x(l,i)
         enddo
        enddo
        dd2(1)=dd2(1)*cke

! dynamics x'=x+gaussn()
        do i=1,nslices
         do l=1,ndim
          x(l,i)=x(l,i)+sqrt(2.d0*hop)*gaussn()
         enddo
        enddo
c compute the ideal action
        dd2(2)=0.d0
        do i=1,nslices
         ip1=iwrap(i+1)
         do l=1,ndim
          dd2(2)=dd2(2)+(x(l,i)-x(l,ip1))**2
         enddo
        enddo
        dd2(2)=dd2(2)*cke

c computes the potential action
        if (omega.ne.0.d0) then
         cold=0.d0
         cnew=0.d0
         do js=1,nslices    !+1
          cold=cold+harmonic(xold(1,js),ndim,fpold(1,js))
          cnew=cnew+harmonic(x(1,js),ndim,fp(1,js))
         enddo
         deltac=-cnew+cold
        else 
         deltac=0.d0
        endif

        ntry=ntry+1
        call rmar(ronf)
        ronf=dlog(max(eps,ronf))
        weight=-dd2(2)+dd2(1)+deltac
        if ( weight.ge.ronf) then
c    ** .. accept **
         nacc = nacc + 1
c update the action
         cpot=cpot-deltac
         action=cpot-dd2(1)+dd2(2)
        else
c    **  .. reject and replace **
         do js=1,nslices
          do l=1,ndim
           x(l,js) = xold(l,js)
          enddo
         enddo
        endif  ! accept or reject

c attempt to displace the entire path
        call rmar(ronf)
        if (ronf.le.gammadsp) call displace
!       if (mod(step,iwigner).eq.0) call wigner
c computing averages
        if (mod(step,nanal).eq.0) then
         tt=etime(tim)
         call analyze(1,step)
         ttanal=ttanal+etime(tim)-tt
        endif
        if (mod(step,nspill).eq.0)
     &      write (31) ((x(l,js),l=1,ndim),js=1,nslices)

       enddo         ! over steps

       return
       end
