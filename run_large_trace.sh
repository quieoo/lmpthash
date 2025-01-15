cd example
echo "*****************************************************************"
echo "Performance Test examples on trace 'web' with 16 threads"
echo "=====================D-Page====================================="
./build/hidpu pagetable 16 configs/trace_web
echo "=====================D-Learned=================================="
./build/hidpu learnedtable 16 configs/trace_web
echo "=====================HiDPU======================================"
./build/hidpu hidpu 16 0 configs/trace_web
echo "=====================LearnedFTL================================="
./build/hidpu learnedftl 16 configs/trace_web
echo "*****************************************************************"
echo "Scalability Test examples on trace 'web', with scale_factor=4"
./build/hidpu scalability 4 configs/trace_web
echo "*****************************************************************"
echo "Reconstruction Test on trace 'web', which rebuilds the index while reading"
./build/hidpu hidpu 16 1 configs/trace_web
echo "*****************************************************************"
echo "Congratulations! You have finished all the examples!"
