#TODO: download datasets
#TODO: decompress datasets
#TODO: fzip compress
#TODO: fzip decompress

.SECONDARY:

all : test.trace.fz

%.fpc :
	wget http://www.csl.cornell.edu/~burtscher/research/FPC/$@

%.bin : %.fpc
	fpc < $< > $@

%.fz : %.bin
	fzip < $< > $@

clean :
	rm *.fz

clean-all :
	rm -rf *.fz *.fpc *.bin
