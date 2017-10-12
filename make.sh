#!/bin/bash

#compile with gfortran
#gfortran -o run test.c vector.c \
#    -L/usr/lib -llapack -lblas \
#    -Wl,-rpath,/usr/lib \
#    -Wl,-rpath,/usr/lib \
#    -std=c99

gcc -I./include HW4_a.c source/*.c\
    -llapack -lblas -std=c99 -o Execute/HW4_a.exe

