echo "*****************************************************************"
echo "Performance Test examples on trace 'web' with 16 threads"

echo "=====================D-Page====================================="
./example/build/hidpu pagetable 16 ./example/configs/trace_web

echo "=====================D-Learned=================================="
./example/build/hidpu learnedtable 16 ./example/configs/trace_web

echo "=====================HiDPU======================================"
./example/build/hidpu hidpu 16 0 ./example/configs/trace_web

echo "=====================LearnedFTL================================="
./example/build/hidpu learnedftl 16 ./example/configs/trace_web

echo "*****************************************************************"
echo "Scalability Test examples on trace 'web', with scale_factor=4"
./example/build/hidpu scalability 4 ./example/configs/trace_web

echo "*****************************************************************"
echo "Reconstruction Test on trace 'web', which rebuilds the index while reading"
./example/build/hidpu hidpu 16 1 ./example/configs/trace_web

echo "*****************************************************************"
echo "Congratulations! You have finished all the examples!"
