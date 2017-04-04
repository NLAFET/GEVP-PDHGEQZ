include ../../make.inc

FSRC =	BDLAGPP.f BDTGEX2.f BDTGEXC.f \
	PDTGSEN.f PDTGORD.f dtgexc.f dtgex2.f \
	dgetc2.f
F90SRC =          \

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


