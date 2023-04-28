        program pimc
c C. Pierleoni, Phys. Dep. University of L'Aquila (ITALY)
c Path Integral MC code for a single quantum particle
c in external potential
        implicit none
        include 'pimc.par'
        integer nb
        real*4 etime,tim(2),te
        include 'pimc.cm'

c setting up the system and the run (read the input file)
        call setsys 
c
c    ** Zero average value **
        call zeroav(0)
c
c compute initial action
        call inaction
        if (nalg.ge.3) call normal_mode

        ttanal=0.
        te=etime(tim)
c    ** Carry on the run **
        do nb=nin+1,nout
         ntry = 0
         nacc = 0
         ntrydsp = 0
         naccdsp = 0
         call zeroav(1)
         if (nstart.eq.-1) call readpc
         if (nstart.ge.0) then
          if (.not.rqmc) then
           if (nalg.eq.0) then
            call levy
           elseif (nalg==3) then
            call fpimc
           elseif(nalg.eq.4) then
            call smc
           elseif (nalg.eq.5) then
            call mc
           elseif (nalg.eq.6) then
            call smc_nm
           else
            write (*,*)'nalg not allowed for closed paths'
            write (*,*)'set nalg=0,3,4,5,6 in the input file'
            stop
           endif
          else
           if (nalg.eq.0) then
            call reptate
           elseif (nalg.eq.1) then
            call bounce
           else
            write (*,*) 'nalg /= 0 or 1 not allowed for opened paths'
            write (*,*) 'set nalg=0 or 1 in the input file and restart'
            stop
           endif
          endif
!         if (nalg.eq.5) call noreversal
!         if (nalg.eq.2) call hybrid_mc
!         if (nalg.eq.3) call hmc_prec
         endif
         if (nacc+naccdsp.eq.0.and.nstart.ge.0) then
          write (*,1166) nb
 1166     format ('block #',i3,': no accepted moves !')
          stop
         endif
c update averages
         call analbl
         anorm(1)=anorm(1)+1.d0
c spill out data
         call spill(nb)
c print out data
         call sprint(nb)
        enddo        ! over blocks
        tt=etime(tim)-te
        write(2,*)
        write(2,*)'Run completed:'
        write(2,*)'startup time (sec.) = ',te
        write(2,*)'run time (sec.) = ',tt
        write(2,*)'analysis time (sec.) = ',ttanal
        write(2,*)'total time (sec.) = ',tt+te
        call rmaget(0)
        stop
        end
