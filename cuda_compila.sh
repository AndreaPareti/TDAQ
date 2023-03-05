#!/bin/bash

if [ $# -eq 0 ]; then
  echo "nccv wrapper"
  echo "usage: $0 file.cu"
  exit 1
fi


EXENAME=`basename $1 .cu`


nvcc --gpu-architecture sm_50 -o $EXENAME $1 $2

if [ $? -eq 0 ]; then
  echo "Created binary: $EXENAME"
fi


