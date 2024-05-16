cd ../build
make clean && make
cd ../example
gcc -o ctest clmpthash_test.c -L../build -lclmpthash -lstdc++ -lm -lpthread
