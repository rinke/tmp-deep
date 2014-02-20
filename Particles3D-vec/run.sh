#!/bin/sh

export LD_LIBRARY_PATH=/opt/intel/composer_xe_2013_sp1.1.106/compiler/lib/mic/
export OMP_NUM_THREADS=$2
export KMP_AFFINITY=$3

exec DEEP/tmp-deep.git/Particles3D-vec/$1
