include ../../make.inc

FSRC =  ACCUMULATECOLUMNROTATIONS.f \
                ACCUMULATEROWROTATIONS.f \
                PDGGHRD.f \
                GRN2LRN.f \
                KRNLCOLUMNANNIHILATE.f \
                KRNLCOLUMNUPDATE.f \
                KRNLROWUPDATE.f \
                KRNLUPDATEANDREDUCECOLUMN.f \
                KRNLACCUMULATEROWROTATIONS.f \
                KRNLACCUMULATECOLUMNROTATIONS.f \
                UPDATEANDREDUCECOLUMN.f \
                UPDATEANDREDUCECOLUMN_ROOT.f \

F90SRC =        SLIVERROWUPDATE.f90     \
                SLIVERHESSCOLUMNUPDATE.f90 \
                SLIVERCOLUMNUPDATE.f90 \
                BLOCKSLIVERROWUPDATE.f90 \
                DUOBLOCKSLIVERCOLUMNUPDATE.f90 \
                BLOCKSLIVERCOLUMNUPDATE.f90

FOBJ = $(FSRC:.f=.o)
F90OBJ = $(F90SRC:.f90=.o)



all: lib
lib: $(FOBJ) $(F90OBJ)

clean:
	rm -f $(FOBJ) $(F90OBJ) 

%.o : %.f
	$(FC) $(FFLAGS) -c $<

%.o : %.f90
	$(FC) $(FFLAGS) -c $<

