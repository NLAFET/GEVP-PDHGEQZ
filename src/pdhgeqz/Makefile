include ../../make.inc

FSRC =  PDHGEQZ.f \
                PDHGEQZ0.f \
                PDHGEQZ1.f \
                PDHGEQZ2.f \
                PDHGEQZ3.f \
                PDHGEQZ4.f \
                PDHGEQZ5.f \
                PDHGEQZ6.f \
                PDHGEQZ7.f PDHGEQZ8.f PDHGEQZ9.f \
                PDHGEQZA.f PDHGEQZB.f \
                PDLASMSUB.f PDLACP4.f \
                pilaenvx.f piparmq.f PDROT.f DHGEQZ5.f DHGEQZ7.f \

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

