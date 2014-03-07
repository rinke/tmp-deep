BIN = Particles3D.cpp.func_to_vectorize
FMA = fma
CPPFLAGS = -vec-report2 -openmp -restrict
DEBUG = -g -debug -Wl,-rpath,/opt/intel/composer_xe_2013_sp1.1.106/compiler/lib/mic/

all: 
	icc $(CPPFLAGS) -O3 -mmic -o $(FMA).mic $(FMA).c 
	icc $(CPPFLAGS) -O3 -mmic -no-vec -o $(FMA).novec.mic $(FMA).c 
	icc $(CPPFLAGS) -O3 -mmic -o $(BIN).mic $(BIN).c 
	icc $(CPPFLAGS) -O3 -mmic -no-vec -o $(BIN).novec.mic $(BIN).c 
	icc $(CPPFLAGS) -O3 -xhost -o $(BIN).cpu $(BIN).c 
	icc $(CPPFLAGS) -O3 -xhost -no-vec -o $(BIN).novec.cpu $(BIN).c 

mpi:
	mpicc -cc=icc $(CPPFLAGS) -O3 -mmic -o $(BIN).mpi.mic $(BIN).mpi.c
	mpicc -cc=icc $(CPPFLAGS) -O3 -mmic -no-vec -o $(BIN).novec.mpi.mic $(BIN).mpi.c
	mpicc -cc=icc $(CPPFLAGS) -O3 -xhost -o $(BIN).mpi.cpu $(BIN).mpi.c
	mpicc -cc=icc $(CPPFLAGS) -O3 -xhost -no-vec -o $(BIN).novec.mpi.cpu $(BIN).mpi.c

debug:
	icc $(CPPFLAGS) $(DEBUG) -O3 -mmic -o $(FMA).dbg.mic $(FMA).c 
	icc $(CPPFLAGS) $(DEBUG) -O3 -mmic -no-vec -o $(FMA).novec.dbg.mic $(FMA).c 
	icc $(CPPFLAGS) $(DEBUG) -O3 -mmic -o $(BIN).dbg.mic $(BIN).c 
	icc $(CPPFLAGS) $(DEBUG) -O3 -mmic -no-vec -o $(BIN).novec.dbg.mic $(BIN).c 
	icc $(CPPFLAGS) $(DEBUG) -O3 -xhost -o $(BIN).dbg.cpu $(BIN).c 
	icc $(CPPFLAGS) $(DEBUG) -O3 -xhost -no-vec -o $(BIN).novec.dbg.cpu $(BIN).c 


clean:
	rm -rf $(BIN){*.mic,*.cpu,*.s} $(FMA){*.mic,*.cpu,*.s}
