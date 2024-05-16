#include <iostream>

#include "include/pthash.hpp"
#include "src/util.hpp"  // for functions distinct_keys and check

#include "pgm/pgm_index.hpp"

void test_pgm_float(){
    printf("test_pgm_float\n");
    // Generate some random data
    std::vector<int> data(1000000);
    std::generate(data.begin(), data.end(), std::rand);

    std::vector<int> copyed;
    for(int i=0;i<data.size();i++){
        copyed.push_back(data[i]);
    }
    std::sort(copyed.begin(), copyed.end());

    // Construct the PGM-index
    // const int epsilon = 128; // space-time trade-off parameter
    pgm::PGMIndex<int,64,4,float> index(copyed, 5, 5);

    // query all points and count average query latency
    // get current time
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t steps=0;
    for(int i=0;i<100000;i++){
        auto range=index.search(data[i]);
        auto lo=copyed.begin()+range.lo;
        bool found=false;
        while(lo!=copyed.end()){
            if(*lo==data[i]){
                found=true;
                break;
            }
            lo++;
            steps++;
        }
        
        if(!found){
            printf("%d-th key: %d not found\n", i, data[i]);
        }
    }


    auto end=std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout<<"Average query latency: "<<(double)duration/(data.size())<<" us"<<std::endl;
    std::cout<<"Average steps: "<<(double)steps/(data.size())<<std::endl;
}
void test_pgm_uint(){
    printf("test_pgm_uint\n");
    // Generate some random data
    std::vector<int> data(1000000);
    std::generate(data.begin(), data.end(), std::rand);

    std::vector<int> copyed;
    for(int i=0;i<data.size();i++){
        copyed.push_back(data[i]);
    }
    std::sort(copyed.begin(), copyed.end());

    // Construct the PGM-index
    // const int epsilon = 128; // space-time trade-off parameter
    pgm::PGMIndex<int,64,4,uint32_t> index(copyed, 5, 5);

    // query all points and count average query latency
    // get current time
    auto start = std::chrono::high_resolution_clock::now();
    uint64_t steps=0;
    for(int i=0;i<100000;i++){
        auto range=index.search(data[i]);
        auto lo=copyed.begin()+range.lo;
        bool found=false;
        while(lo!=copyed.end()){
            if(*lo==data[i]){
                found=true;
                break;
            }
            lo++;
            steps++;
        }
        
        if(!found){
            printf("%d-th key: %d not found\n", i, data[i]);
        }
    }


    auto end=std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    std::cout<<"Average query latency: "<<(double)duration/(data.size())<<" us"<<std::endl;
    std::cout<<"Average steps: "<<(double)steps/(data.size())<<std::endl;
}
void test_pgm(){
    test_pgm_float();
    test_pgm_uint();
}



int main() {
    using namespace pthash;

    test_pgm();
    return 0;

    /* Generate 10M random 64-bit keys as input data. */
    static const uint64_t num_keys = 1000;
    static const uint64_t seed = 1234567890;
    std::cout << "generating input data..." << std::endl;
    std::vector<uint64_t> keys = distinct_keys<uint64_t>(num_keys, seed);
    assert(keys.size() == num_keys);

    /* Set up a build configuration. */
    build_configuration config;
    config.c = 6.0;
    config.alpha = 0.94;
    config.minimal_output = true;  // mphf
    config.verbose_output = false;

    config.LinearMapping=true;
    config.max_bucket_size = 100;
    config.pilot_search_threshold = 10000000;

    /* Declare the PTHash function. */
    typedef single_phf<murmurhash2_64,         // base hasher
                       pthash::compact,  // encoder type
                       false                    // minimal
                       >
        pthash_type;

    // config.num_partitions = 50;
    // config.num_threads = 4;
    // typedef partitioned_mphf<murmurhash2_64,        // base hasher
    //                          dictionary_dictionary  // encoder type
    //                          >
    //     pthash_type;

    pthash_type f;

    /* Build the function in internal memory. */
    std::cout << "building the function..." << std::endl;
    auto start = clock_type::now();
    auto timings = f.build_in_internal_memory(keys.begin(), keys.size(), config);
    // auto timings = f.build_in_external_memory(keys.begin(), keys.size(), config);
    double total_seconds = timings.partitioning_seconds + timings.mapping_ordering_seconds +
                           timings.searching_seconds + timings.encoding_seconds;
    std::cout << "function built in " << seconds(clock_type::now() - start) << " seconds"<< std::endl;
    std::cout << "computed: " << total_seconds << " seconds" << std::endl;
    /* Compute and print the number of bits spent per key. */
    double bits_per_key = static_cast<double>(f.num_bits()) / f.num_keys();
    std::cout << "function uses " << bits_per_key << " [bits/key]" << std::endl;

    /* Sanity check! */
    if (check(keys.begin(), f)) std::cout << "EVERYTHING OK!" << std::endl;

    /* Now evaluate f on some keys. */
    for (uint64_t i = 0; i != 10; ++i) {
        std::cout << "f(" << keys[i] << ") = " << f(keys[i]) << '\n';
    }

    /* Serialize the data structure to a file. */
    std::cout << "serializing the function to disk..." << std::endl;
    std::string output_filename("pthash.bin");
    essentials::save(f, output_filename.c_str());

    {
        /* Now reload from disk and query. */
        pthash_type other;
        essentials::load(other, output_filename.c_str());
        for (uint64_t i = 0; i != 10; ++i) {
            std::cout << "f(" << keys[i] << ") = " << other(keys[i]) << '\n';
            assert(f(keys[i]) == other(keys[i]));
        }
    }

    std::remove(output_filename.c_str());
    return 0;
}