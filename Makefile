#############################################
#Makefile for MultiWell and MultiWell Tools.#
#              17/04/2019                   #
#                                           #
#############################################
all:    
	echo '';\
	echo 'Multiwell and all of the helper applications will now be compiled';\
	echo '';\
	date;\
	cd src/multiwell;\
	make multiwell;\
	cd ../densum;\
	make densum;\
	cd ../thermo;\
	make thermo;\
	cd ../mominert;\
	make mominert;\
	cd ../X2multi;\
	make gauss2multi;\
	cd ../bdens;\
	make bdens;\
	cd ../ktools;\
	make ktools;\
	cd ../sctst;\
	make sctst;\
	cd ../parsctst;\
	make parsctst;\
	cd ../lamm;\
	make lamm;\
	cd ../paradensum;\
	make paradensum;\
	cd ../;\
	cp gauss2lamm.sh ../bin;

clean:
	cd src/multiwell;\
	make clean;\
	cd ../densum;\
	make clean;\
	cd ../thermo;\
	make clean;\
	cd ../mominert;\
	make clean;\
	cd ../X2multi;\
	make clean;\
	cd ../bdens;\
	make clean;\
	cd ../ktools;\
	make clean;\
	cd ../sctst;\
	make clean;\
	cd ../parsctst;\
	make clean;\
	cd ../lamm;\
	make clean;\
	cd ../paradensum;\
	make clean;\
	cd ../../bin;\
	rm gauss2lamm.sh;
