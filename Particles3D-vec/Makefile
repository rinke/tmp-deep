BIN = Particles3D.cpp.func_to_vectorize
FMA = fma
CPPFLAGS = -vec-report2 -openmp

all: 
	icc $(CPPFLAGS) -O3 -mmic -o $(FMA).mic $(FMA).c 
	icc $(CPPFLAGS) -O3 -mmic -no-vec -o $(FMA).novec.mic $(FMA).c 
	icc $(CPPFLAGS) -O3 -mmic -o $(BIN).mic $(BIN).c 
	icc $(CPPFLAGS) -O3 -mmic -no-vec -o $(BIN).novec.mic $(BIN).c 
	icc $(CPPFLAGS) -O3 -xhost -o $(BIN).cpu $(BIN).c 
	icc $(CPPFLAGS) -O3 -xhost -no-vec -o $(BIN).novec.cpu $(BIN).c 

clean:
	rm -rf $(BIN){*.mic,*.cpu} $(FMA){*.mic,*.cpu}
