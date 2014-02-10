#!/bin/sh

export LD_LIBRARY_PATH=/opt/intel/composer_xe_2013.4.183/compiler/lib/mic/
export OMP_NUM_THREADS=$2
export KMP_AFFINITY=$3

DEEP/tmp-deep.git/Particles3D-vec/$1
