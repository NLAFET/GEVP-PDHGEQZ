include ../../make.inc

FSRC = 	qzdackit.f qzinf.f kkqz.f qzconf.f \
	qzearly.f qzlap.f kkqzconf.f \
	invhse.f qzdack.f qzfcol.f qzlapit.f

F90SRC =	\

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

