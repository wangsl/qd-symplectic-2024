#!/bin/bash

args=
for i in "$@"; do 
  i="${i//\\/\\\\}"
  args="${args} \"${i//\"/\\\"}\""
done

if [[ "${args}" == "" ]]; then args="/bin/bash"; fi

if [[ -e /dev/nvidia0 ]]; then nv="--nv"; fi

if [[ "${SINGULARITY_CONTAINER}" != "" ]]; then export PATH="/share/apps/apptainer/bin:${PATH}"; fi

singularity exec ${nv} \
--overlay /scratch/work/public/singularity/cuda-11.8.89-cudnn-8.7.sqf:ro \
--overlay /scratch/work/public/singularity/cuda-samples.sqf:ro \
--overlay /scratch/work/public/singularity/matlab-2024a.sqf:ro \
/scratch/work/public/singularity/rocm5.7.1-ubuntu22.04.3.sif \
/bin/bash -c "
unset -f which
export MODULEPATH=/ext3/apps/modulefiles
module load cuda/11.8.89 
export CUDA_SAMPLES_INC=/ext3/cuda/samples/Common
export PATH=/ext3/apps/matlab/2024a/bin:\${PATH}
export MATLAB_ROOT=/ext3/apps/matlab/2024a
export MATLAB_JAVA=/ext3/apps/matlab/openjdk/1.8.0/jre
export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libstdc++.so.6
${args}
"
