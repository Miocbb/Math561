cmake_minimum_required(VERSION 3.0.2)
project(Math_561)

include_directories(include)
set(SOURCE GS_QR.c)
file(GLOB SOURCE "source/*.c")

set(LIBRARY liblapack.so libm.so)
set(LIB_PATH /usr/lib/lapack /usr/lib/x86_64-linux-gnu)
link_directories(${LIB_PATH})

add_executable(GS_QR.exe ${SOURCE})
target_link_libraries(GS_QR.exe ${LIBRARY})
