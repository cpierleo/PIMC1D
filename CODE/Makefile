#####   INTEL
FC = ifort -g -C -r8 -g  -w -fpe0 -traceback -trapuv
FC = ifort -O3 -r8 -g  -w -fpe0 -traceback -trapuv
#COMP = 
XTRALIB=
########   MacBook ##############
#FC = gfortran -C
#FC = gfortran -O3
#COMP = -g
############

FFLAGS = $(COMP) 
LDFLAGS = $(XTRALIB) #-Bstatic

DEPEND	= Makefile pimc.par pimc.cm 
OBJS = pimc.o zeroav.o cumul.o analyze.o analbl.o spill.o \
       sprint.o cumul1.o readpc.o reptate.o levy.o lpsit.o \
       displace.o setsys.o inaction.o random.o harmonic.o hybrid_mc.o hforce.o \
       wigner.o forcew.o wdynamics.o velverl.o quanti.o \
       wsprint.o polin2.o polint.o bounce.o noreversal.o reaction_coordinate.o \
       normal_mode.o fpimc.o reduce.o reduceT.o smc.o mc.o smc_nm.o #wigner.o forcew.o wdynamics.o velverl.o quanti.o

TOUT  =  pimc correl findrlv pc2xyz
all : $(TOUT)


pimc : $(OBJS) $(DEPEND)
	$(FC) $(OBJS) -o pimc $(LDFLAGS)

pc2xyz : pc2xyz.o react.o
	$(FC) $(LDFLAGS) pc2xyz.o react.o -o pc2xyz 

seers : seers.o 
	$(FC) seers.o -o seers

chrs : chrs.o 
	$(FC) chrs.o -o chrs

correl : correl.o $(DEPEND)
	$(FC) $(LDFLAGS) correl.o -o correl

sergio : sergio.o $(DEPEND)
	$(FC) sergio.o -o sergio

findrlv : findrlv.o sites.o shells.o cossin.o functions.o $(DEPEND)
	$(FC) findrlv.o sites.o shells.o cossin.o functions.o -o findrlv

clean : 
	rm -f *.o  $(TOUT)
  
$(OBJS) : $(DEPEND)
correl.o : $(DEPEND)
findrlv.o : $(DEPEND)
pc2xyz.o : $(DEPEND)
react.o : $(DEPEND)
