include ../make.inc

OBJS =	pdhgeqz/PDHGEQZB.o \
	pdhgeqz/DHGEQZ5.o pdhgeqz/DHGEQZ7.o pdhgeqz/PDHGEQZA.o  \
	pdhgeqz/PDHGEQZ0.o pdhgeqz/PDHGEQZ1.o pdhgeqz/PDHGEQZ2.o \
	pdhgeqz/PDHGEQZ3.o pdhgeqz/PDHGEQZ4.o pdhgeqz/PDLACP4.o \
	pdhgeqz/PDHGEQZ6.o pdhgeqz/PDHGEQZ7.o pdhgeqz/PDHGEQZ8.o \
	pdhgeqz/PDHGEQZ9.o pdhgeqz/PDROT.o pdhgeqz/PDLASMSUB.o \
	pdhgeqz/pilaenvx.o pdhgeqz/piparmq.o pdhgeqz/PDHGEQZ.o \
	pdhgeqz/PDHGEQZ5.o reorder/BDLAGPP.o reorder/BDTGEX2.o \
	reorder/BDTGEXC.o reorder/PDTGSEN.o reorder/PDTGORD.o \
	reorder/dtgexc.o reorder/dtgex2.o reorder/dgetc2.o \
	kkqz/qzdackit.o kkqz/qzinf.o kkqz/kkqz.o kkqz/qzconf.o \
	kkqz/qzearly.o kkqz/qzlap.o kkqz/kkqzconf.o kkqz/invhse.o \
	kkqz/qzdack.o kkqz/qzfcol.o kkqz/qzlapit.o kkqz/qzconf.o \
	ab2ht/ACCUMULATECOLUMNROTATIONS.o ab2ht/ACCUMULATEROWROTATIONS.o \
	ab2ht/PDGGHRD.o ab2ht/GRN2LRN.o ab2ht/KRNLCOLUMNANNIHILATE.o \
	ab2ht/KRNLCOLUMNUPDATE.o ab2ht/KRNLROWUPDATE.o \
	ab2ht/KRNLUPDATEANDREDUCECOLUMN.o ab2ht/KRNLACCUMULATEROWROTATIONS.o \
	ab2ht/KRNLACCUMULATECOLUMNROTATIONS.o ab2ht/UPDATEANDREDUCECOLUMN.o \
	ab2ht/UPDATEANDREDUCECOLUMN_ROOT.o ab2ht/SLIVERHESSCOLUMNUPDATE.o \
	ab2ht/SLIVERHESSCOLUMNUPDATE.o ab2ht/SLIVERCOLUMNUPDATE.o \
	ab2ht/SLIVERROWUPDATE.o ab2ht/BLOCKSLIVERCOLUMNUPDATE.o \
	ab2ht/BLOCKSLIVERROWUPDATE.o ab2ht/DUOBLOCKSLIVERCOLUMNUPDATE.o

LIB = ../libpdhgeqz.a

all: lib
lib: 	
	$(MAKE) -C ./ab2ht 
	$(MAKE) -C ./pdhgeqz
	$(MAKE) -C ./reorder
	$(MAKE) -C ./kkqz
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

DISTURBING = *~
clean:
	rm -rf $(LIB) 
	rm -rf ab2ht/*.o
	rm -rf kkqz/*.o
	rm -rf reorder/*.o
	rm -rf pdhgeqz/*.o

