cmake_minimum_required(VERSION 3.1)
project(triangle)

add_compile_options(-O -DLINUX -I/usr/X11R6/include -L/usr/X11R6/lib)
add_compile_options(-DTRILIBRARY -DANSI_DECLARATORS)
add_library(triangle
    triangle.c
    triangle.h
)

