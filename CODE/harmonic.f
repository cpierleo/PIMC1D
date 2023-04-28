       function harmonic(xy,ldim,fff)
       implicit none
c harmonic external potential
       include 'pimc.par'
       integer il,ldim,l
       real*8 xy(ldim),harmonic,fff(ldim)
       include 'pimc.cm'

       harmonic=0.d0
       do l=1,ldim
        if (model == 1) then                        ! harmonic oscillator
         harmonic=harmonic+potpar*xy(l)*xy(l)
         fff(l)=-2.d0*potpar*xy(l)
        elseif(model == 2) then                     ! double well potential
         harmonic=harmonic+potpar*2.d0*(xy(l)**2-4.d0)**2
         fff(l)=-potpar*8.d0*xy(l)*(xy(l)-4.d0)
        else
         write(*,*) 'model system not defined',model
         stop
        endif
       enddo
       return
       end

       function harmvir(x1,xy,ldim)
       implicit none
c harmonic external potential
       include 'pimc.par'
       integer il,ldim,l
       real*8 xy(ldim),harmvir,x1(ldim)
       include 'pimc.cm'

       harmvir=0.d0
       do l=1,ldim
        if (model == 1) then
         harmvir=harmvir+2.d0*potpar*xy(l)*(xy(l)-x1(l))
        elseif(model == 2) then
         harmvir=harmvir+8.d0*potpar*xy(l)*(xy(l)**2-4.d0)*(xy(l)-x1(l))
        endif
       enddo
       return
       end
