#!/bin/bash

#compile with gfortran
#gfortran -o run test.c vector.c \
#    -L/usr/lib -llapack -lblas \
#    -Wl,-rpath,/usr/lib \
#    -Wl,-rpath,/usr/lib \
#    -std=c99

gcc -I./include GS_QR.c source/*.c\
    -llapack -lblas -std=c99 -o Execute/GS_QR.exe

