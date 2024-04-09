#!/bin/bash

hip_dir="/scratch/work/vip-amd-cuda/qd-symplectic/common-hip"

cd /home/wang/hudson-20211220/src/common-cuda

export MATLAB_ROOT=/ext3/apps/matlab/2024a

for src in *.[hC] *.cu *.F; do
	echo ${src}
	hip_src="${hip_dir}/${src}"
	hipify-clang --cuda-path=${CUDA_HOME} -I${CUDA_SAMPLES_INC} -I${MATLAB_ROOT}/extern/include ${src} -o ${hip_src}
done

