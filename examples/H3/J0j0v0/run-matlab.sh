#!/bin/bash

source /etc/profile.d/lmod.sh

module purge
module load matlab/2020a
module load cuda/11.0.2
module load intel/19.1.1

export CUDA_VISIBLE_DEVICES=0

if [ "$SLURM_JOB_TMPDIR" != "" ]; then
    base_dir=$SLURM_JOB_TMPDIR
else
    base_dir=/state/partition1/$USER
fi

mkdir -p $base_dir

export MATLAB_PREFDIR=$(mktemp -d $base_dir/matlab-XXXX)

rm -rf crp/* potwave/*

matlab -nodisplay -r "main(0,0); exit"

exit

export CUDA_VISIBLE_DEVICES=2

screen -d -m \
       bash -c 'taskset -c $(seq -s, 0 2 40) \
       singularity \
       exec --nv --bind $MATLAB_PREFDIR:$HOME/.matlab \
       /home/wang/singularity-images/centos-7.6.1810.simg \
       /share/apps/matlab/2017b/bin/matlab -nodisplay -r \
       "main(0,0); exit" > t2.log 2>&1'

exit

taskset -c $(seq -s, 0 2 27) matlab -nodisplay -r "main(0,0); exit" > t1.log 2>&1 &

wait

exit

n=0
while read -r line; do
    rm -rf crp/* potwave/*
    echo $line | awk '{print $3 " " $4 " " $5 " " $6 " " $7 " " $8}' > a.txt
    n=$((n+1))
    taskset -c $(seq -s, 0 2 27) matlab -nodisplay -r "main(0,0); exit" 2>&1 | tee $n.log
done <t1.txt

exit

rm -rf crp/* potwave/*

export CUDA_VISIBLE_DEVICES=2

#taskset -c $(seq -s, 0 2 27) matlab #-nodesktop

taskset -c $(seq -s, 0 2 27) matlab -nodisplay -r "main(0,0); exit" > t1.log 2>&1 &

wait
