#!/bin/bash

module purge
module load matlab/2017b

export CUDA_VISIBLE_DEVICES=0

#taskset -c $(seq -s, 0 2 27) matlab #-nodesktop

taskset -c $(seq -s, 0 2 27) matlab -nodesktop -r "main; exit" > t1.log 2>&1 &

wait









