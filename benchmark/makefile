#TODO: download datasets
#TODO: decompress datasets
#TODO: fzip compress
#TODO: fzip decompress

.SECONDARY:

all : test.trace.fz

%.fpc :
	wget http://www.cs.txstate.edu/~burtscher/research/datasets/FPdouble/$@

%.bin : %.fpc
	fpc < $< > $@

%.fz : %.bin
	fzip < $< > $@

clean :
	rm *.fz

clean-all :
	rm -rf *.fz *.fpc *.bin
