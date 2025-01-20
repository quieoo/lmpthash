echo "*****************************************************************"
echo "Performance Test examples on trace 'wdev' with 16 threads"

echo "=====================D-Page====================================="
./example/build/hidpu pagetable 16 ./example/configs/trace_wdev

echo "=====================D-Learned=================================="
./example/build/hidpu learnedtable 16 ./example/configs/trace_wdev

echo "=====================HiDPU======================================"
./example/build/hidpu hidpu 16 0 ./example/configs/trace_wdev

echo "=====================LearnedFTL================================="
./example/build/hidpu learnedftl 16 ./example/configs/trace_wdev

echo "*****************************************************************"
echo "Scalability Test examples on trace 'wdev', with scale_factor=4"
./example/build/hidpu scalability 4 ./example/configs/trace_wdev

echo "*****************************************************************"
echo "Reconstruction Test on trace 'wdev', which rebuilds the index while reading"
./example/build/hidpu hidpu 16 1 ./example/configs/trace_wdev

echo "*****************************************************************"
echo "Congratulations! You have finished all the examples!"
