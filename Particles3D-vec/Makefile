BIN = Particles3D.cpp.func_to_vectorize

all: 
	icc -O3 -mmic -o $(BIN).mic $(BIN).c 
	icc -O3 -mmic -no-vec -o $(BIN).novec.mic $(BIN).c 
	icc -O3 -xhost -o $(BIN).cpu $(BIN).c 
	icc -O3 -xhost -no-vec -o $(BIN).novec.cpu $(BIN).c 

clean:
	rm -rf $(BIN){*.mic,*.cpu}
