default:	kuramoto_classic kuramoto_diamond

COPT = -O3 -std=c++11
# COPT = -O0 -g -std=c++11

kuramoto_classic:	kuramoto_classic.o
	mpic++ $(COPT) -o kuramoto_classic kuramoto_classic.o -L/master/home/qiqi/lib -lpng

kuramoto_diamond:    kuramoto_diamond.o
	mpic++ $(COPT) -o kuramoto_diamond kuramoto_diamond.o -L/master/home/qiqi/lib -lpng

kuramoto_classic.o:  kuramoto_classic.cpp kuramoto.h pde_classic.h pde_common.h PngWriter.hpp
	mpic++ $(COPT) -I/master/home/qiqi/include/libpng16 -c $< 

kuramoto_diamond.o:  kuramoto_diamond.cpp kuramoto.h pde_diamond.h pde_common.h PngWriter.hpp
	mpic++ $(COPT) -I/master/home/qiqi/include/libpng16 -c $< 

run:	kuramoto_classic kuramoto_diamond
	# mpirun -np 4 ./kuramoto_classic 8
	# mpirun -np 4 ./kuramoto_classic 16
	# mpirun -np 4 ./kuramoto_classic 32
	# mpirun -np 4 ./kuramoto_classic 64
	# mpirun -np 4 ./kuramoto_classic 128
	# mpirun -np 4 ./kuramoto_classic 256
	# mpirun -np 4 ./kuramoto_classic 512
	# mpirun -np 4 ./kuramoto_classic 1024
	# mpirun -np 4 ./kuramoto_classic 2048
	# mpirun -np 4 ./kuramoto_classic 4096
	# mpirun -np 4 ./kuramoto_classic 8192
	mpirun -np 4 ./kuramoto_diamond 8
	mpirun -np 4 ./kuramoto_diamond 16
	mpirun -np 4 ./kuramoto_diamond 32
	mpirun -np 4 ./kuramoto_diamond 64
	mpirun -np 4 ./kuramoto_diamond 128
	mpirun -np 4 ./kuramoto_diamond 256
	mpirun -np 4 ./kuramoto_diamond 512
	mpirun -np 4 ./kuramoto_diamond 1024
	mpirun -np 4 ./kuramoto_diamond 2048
	mpirun -np 4 ./kuramoto_diamond 4096
	mpirun -np 4 ./kuramoto_diamond 8192

clean:
	rm -f *.o kuramoto_diamond kuramoto_classic
