#!/bin/bash

#compile with gfortran
#gfortran -o run test.c vector.c \
#    -L/usr/lib -llapack -lblas \
#    -Wl,-rpath,/usr/lib \
#    -Wl,-rpath,/usr/lib \
#    -std=c99

gcc -I./include HW/PP_LU.c source/*.c\
    -llapack -lblas -std=c99 -o Execute/PP_LU.exe

