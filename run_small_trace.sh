cd example
echo "*****************************************************************"
echo "Performance Test examples on trace 'wdev' with 16 threads"
echo "=====================D-Page====================================="
./build/hidpu pagetable 16 configs/trace_wdev
echo "=====================D-Learned=================================="
./build/hidpu learnedtable 16 configs/trace_wdev
echo "=====================HiDPU======================================"
./build/hidpu hidpu 16 0 configs/trace_wdev
echo "=====================LearnedFTL================================="
./build/hidpu learnedftl 16 configs/trace_wdev
echo "*****************************************************************"
echo "Scalability Test examples on trace 'wdev', with scale_factor=32"
./build/hidpu scalability 4 configs/trace_wdev
echo "*****************************************************************"
echo "Reconstruction Test on trace 'wdev', which rebuilds the index while reading"
./build/hidpu hidpu 16 1 configs/trace_wdev
echo "*****************************************************************"
echo "Congratulations! You have finished all the examples!"
