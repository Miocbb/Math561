#!/bin/bash

#compile with gfortran
#gfortran -o run test.c vector.c \
#    -L/usr/lib -llapack -lblas \
#    -Wl,-rpath,/usr/lib \
#    -Wl,-rpath,/usr/lib \
#    -std=c99

<<<<<<< HEAD
<<<<<<< HEAD
gcc -I./include HW4_a.c source/*.c\
    -llapack -lblas -std=c99 -o Execute/HW4_a.exe
=======
gcc -I./include GS_QR.c source/*.c\
    -llapack -lblas -std=c99 -o Execute/GS_QR.exe
>>>>>>> 87995adf92c9459c817414c28f46e073c4292b50
=======
gcc -I./include HW4_a.c source/*.c\
    -llapack -lblas -std=c99 -o Execute/HW4_a.exe
>>>>>>> a2810e620ab355e2fedc123bd29479b58a6079dd

