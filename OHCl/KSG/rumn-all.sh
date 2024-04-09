#!/bin/bash

module purge
module load matlab/2017b

rm -rf potwave/* crp/*

export CUDA_VISIBLE_DEVICES=1
taskset -c $(seq -s, 0 2 27) matlab -nodisplay -r "main; exit" > t1.log 2>&1 &

wait

