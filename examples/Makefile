include ../make.inc

OBJS =	\
	../tools/PDLAPRNT.o \
	../tools/DLAPRNT.o \
	../tools/pmatgeninc.o \
	../tools/pdmatgen2.o \
	../tools/PQZHLPS.o \
	../tools/PDRDMM.o \
	../tools/mmio.o

LIBS = ../libpdhgeqz.a $(SCALAPACK) $(LAPACK) $(BLAS)

all: $(OBJS) t1 t2 t3

 
DISTURBING = *~

t1: $(OBJS) EXRAND1.o
	$(FC) $(FFLAGS) -o EXRAND1 EXRAND1.o $(OBJS) $(LIBS)
t2: $(OBJS) EXRAND2.o
	$(FC) $(FFLAGS) -o EXRAND2 EXRAND2.o $(OBJS) $(LIBS)
t3: $(OBJS) EXRAND3.o
	$(FC) $(FFLAGS) -o EXRAND3 EXRAND3.o $(OBJS) $(LIBS)
.PHONY : all

clean:
	/bin/rm -f -v EXRAND1 EXRAND2 EXRAND3 *.o
	/bin/rm -f -v $(OBJS)
	/bin/rm -f -v $(DISTURBING)

