#!/bin/bash

hip_dir=/scratch/work/vip-amd-cuda/qd-symplectic/cuda/common-hip

cd /scratch/work/vip-amd-cuda/qd-symplectic/cuda/common-cuda

for src in *.[hC] *.cu *.F; do
	echo ${src}
	hip_src="${hip_dir}/${src}"
	hipify-perl ${src} > ${hip_src}
done

