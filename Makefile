default:	kuramoto_classic kuramoto_swept heat_classic heat_swept

COPT = -O3 -std=c++11 -Wall
# COPT = -O0 -g -std=c++11 -Wall

# --------------------- HEAT --------------------#

heat_classic:	heat_classic.o
	mpic++ $(COPT) -o heat_classic heat_classic.o -L/master/home/qiqi/lib -lpng

heat_swept:    heat_swept.o
	mpic++ $(COPT) -o heat_swept heat_swept.o -L/master/home/qiqi/lib -lpng

heat_classic.o:  heat_classic.cpp heat.h pde_classic.h pde_common.h PngWriter.hpp
	mpic++ $(COPT) -I/master/home/qiqi/include/libpng16 -c $< 

heat_swept.o:  heat_swept.cpp heat.h pde_swept.h pde_common.h PngWriter.hpp
	mpic++ $(COPT) -I/master/home/qiqi/include/libpng16 -c $< 

run_heat:	heat_classic heat_swept
	mpirun -np 4 ./heat_classic 8
	mpirun -np 4 ./heat_classic 16
	mpirun -np 4 ./heat_classic 32
	mpirun -np 4 ./heat_classic 64
	mpirun -np 4 ./heat_classic 128
	mpirun -np 4 ./heat_classic 256
	mpirun -np 4 ./heat_classic 512
	mpirun -np 4 ./heat_classic 1024
	mpirun -np 4 ./heat_classic 2048
	mpirun -np 4 ./heat_classic 4096
	mpirun -np 4 ./heat_classic 8192
	mpirun -np 4 ./heat_classic 8192
	mpirun -np 4 ./heat_classic 16384
	mpirun -np 4 ./heat_classic 32769
	mpirun -np 4 ./heat_classic 65536
	mpirun -np 4 ./heat_swept 8
	mpirun -np 4 ./heat_swept 16
	mpirun -np 4 ./heat_swept 32
	mpirun -np 4 ./heat_swept 64
	mpirun -np 4 ./heat_swept 128
	mpirun -np 4 ./heat_swept 256
	mpirun -np 4 ./heat_swept 512
	mpirun -np 4 ./heat_swept 1024
	mpirun -np 4 ./heat_swept 2048
	mpirun -np 4 ./heat_swept 4096
	mpirun -np 4 ./heat_swept 8192
	mpirun -np 4 ./heat_swept 16384
	mpirun -np 4 ./heat_swept 32769
	mpirun -np 4 ./heat_swept 65536

# ------------------- KURAMOTO ------------------#

kuramoto_classic:	kuramoto_classic.o
	mpic++ $(COPT) -o kuramoto_classic kuramoto_classic.o -L/master/home/qiqi/lib -lpng

kuramoto_swept:    kuramoto_swept.o
	mpic++ $(COPT) -o kuramoto_swept kuramoto_swept.o -L/master/home/qiqi/lib -lpng

kuramoto_classic.o:  kuramoto_classic.cpp kuramoto.h pde_classic.h pde_common.h PngWriter.hpp
	mpic++ $(COPT) -I/master/home/qiqi/include/libpng16 -c $< 

kuramoto_swept.o:  kuramoto_swept.cpp kuramoto.h pde_swept.h pde_common.h PngWriter.hpp
	mpic++ $(COPT) -I/master/home/qiqi/include/libpng16 -c $< 

run_kuramoto:	kuramoto_classic kuramoto_swept
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
	mpirun -np 4 ./kuramoto_swept 8
	mpirun -np 4 ./kuramoto_swept 16
	mpirun -np 4 ./kuramoto_swept 32
	mpirun -np 4 ./kuramoto_swept 64
	mpirun -np 4 ./kuramoto_swept 128
	mpirun -np 4 ./kuramoto_swept 256
	mpirun -np 4 ./kuramoto_swept 512
	mpirun -np 4 ./kuramoto_swept 1024
	mpirun -np 4 ./kuramoto_swept 2048
	mpirun -np 4 ./kuramoto_swept 4096
	# mpirun -np 4 ./kuramoto_swept 8192
	# mpirun -np 4 ./kuramoto_swept 16384
	# mpirun -np 4 ./kuramoto_swept 32769
	# mpirun -np 4 ./kuramoto_swept 65536

clean:
	rm -f *.o kuramoto_swept kuramoto_classic heat_swept heat_classic
