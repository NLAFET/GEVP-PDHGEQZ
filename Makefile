all: lib ex tests

lib:
	( cd src; make > /dev/null; cd ..)

tests:
	( cd testing; make > /dev/null; ./runmpi.sh; cd .. )
ex:
	(cd examples; make > /dev/null; cd ..)
clean:
	( cd src; make clean; \
        cd ../testing; make clean; cd ../examples; make clean;\
        cd .. )

