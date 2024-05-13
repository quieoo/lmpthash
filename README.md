<!-- [![CodeQL](https://github.com/jermp/pthash/actions/workflows/codeql.yml/badge.svg)](https://github.com/jermp/pthash/actions/workflows/codeql.yml) -->


Introduction
----
Learned Monotonic Perfect Hashing, combining the PTHash and PGM-Index libraries. 

Convert as many consecutive keys as possible into precise segments or PThash based on linear mapping to reduce query overhead.

## Build lmpthash
	mkdir build
	cd build
	cmake ..
	make

Generated files:

`build/utils`: utility programs build on command line

`build/libclmpthash.a`: library for linking with C programs

`include/clmpthash.h`: header file for C programs


## Example of including clmpthash in C source file
	cd example
	gcc -o ctest clmpthash_test.c -L../build -lclmpthash -lstdc++ -lm -lpthread
