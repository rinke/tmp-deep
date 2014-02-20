#!/bin/sh

export LD_LIBRARY_PATH=/opt/intel/composer_xe_2013_sp1.1.106/compiler/lib/mic/
export OMP_NUM_THREADS=1
export KMP_AFFINITY=balanced

exec DEEP/tmp-deep.git/Particles3D-vec/Particles3D.cpp.func_to_vectorize.mic
