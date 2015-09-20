#TODO: download datasets
#TODO: decompress datasets
#TODO: fzip compress
#TODO: fzip decompress

.SECONDARY:


test : test.trace.fz

obs : obs_error.trace.fz obs_info.trace.fz obs_spitzer.trace.fz obs_temp.trace.fz

num : num_brain.trace.fz num_comet.trace.fz num_control.trace.fz num_plasma.trace.fz

msg : msg_bt.trace.fz msg_lu.trace.fz msg_sp.trace.fz msg_sppm.trace.fz msg_sweep3d.trace.fz

all : test obs

%.fpc :
	wget http://www.cs.txstate.edu/~burtscher/research/datasets/FPdouble/$@

%.bin : %.fpc
	fpc < $< > $@

%.fz : %.bin
	../fzip < $< > $@

vim :
	vim -p makefile createBin.cpp

createBin: createBin.cpp
	g++ -std=c++11 -o createBin createBin.cpp

example.bin: createBin
	createBin < example.dat > example.bin

clean :
	rm *.fz

clean-all :
	rm -rf *.fz *.fpc *.bin
