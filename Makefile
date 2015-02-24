default:	kuramoto_classic kuramoto_swept

COPT = -O3 -std=c++11
# COPT = -O0 -g -std=c++11

kuramoto_classic:	kuramoto_classic.o
	mpic++ $(COPT) -o kuramoto_classic kuramoto_classic.o -L/master/home/qiqi/lib -lpng

kuramoto_swept:    kuramoto_swept.o
	mpic++ $(COPT) -o kuramoto_swept kuramoto_swept.o -L/master/home/qiqi/lib -lpng

kuramoto_classic.o:  kuramoto_classic.cpp kuramoto.h pde_classic.h pde_common.h PngWriter.hpp
	mpic++ $(COPT) -I/master/home/qiqi/include/libpng16 -c $< 

kuramoto_swept.o:  kuramoto_swept.cpp kuramoto.h pde_swept.h pde_common.h PngWriter.hpp
	mpic++ $(COPT) -I/master/home/qiqi/include/libpng16 -c $< 

run:	kuramoto_classic kuramoto_swept
	# mpirun ./kuramoto_classic 8
	# mpirun ./kuramoto_classic 16
	# mpirun ./kuramoto_classic 32
	# mpirun ./kuramoto_classic 64
	# mpirun ./kuramoto_classic 128
	# mpirun ./kuramoto_classic 256
	# mpirun ./kuramoto_classic 512
	# mpirun ./kuramoto_classic 1024
	# mpirun ./kuramoto_classic 2048
	# mpirun ./kuramoto_classic 4096
	# mpirun ./kuramoto_classic 8192
	mpirun ./kuramoto_swept 8
	mpirun ./kuramoto_swept 16
	mpirun ./kuramoto_swept 32
	mpirun ./kuramoto_swept 64
	mpirun ./kuramoto_swept 128
	mpirun ./kuramoto_swept 256
	mpirun ./kuramoto_swept 512
	mpirun ./kuramoto_swept 1024
	mpirun ./kuramoto_swept 2048
	mpirun ./kuramoto_swept 4096
	mpirun ./kuramoto_swept 8192
	mpirun ./kuramoto_swept 16384
	mpirun ./kuramoto_swept 32769
	mpirun ./kuramoto_swept 65536

clean:
	rm -f *.o kuramoto_swept kuramoto_classic
