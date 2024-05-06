#include <iostream>

#include "include/pthash.hpp"
#include "src/util.hpp"  // for functions distinct_keys and check

#include "pgm/pgm_index.hpp"

void test_pgm(){
    // Generate some random data
    std::vector<int> data(1000000);
    std::generate(data.begin(), data.end(), std::rand);
    data.push_back(42);
    std::sort(data.begin(), data.end());

    // Construct the PGM-index
    const int epsilon = 128; // space-time trade-off parameter
    pgm::PGMIndex<int, epsilon> index(data);

    // Query the PGM-index
    auto q = 42;
    auto range = index.search(q);
    auto lo = data.begin() + range.lo;
    auto hi = data.begin() + range.hi;
    std::cout << *std::lower_bound(lo, hi, q)<<std::endl;
    std::cout<<"pgm index test done"<<std::endl;

    // query all points and count average query latency
    
}

int main() {
    using namespace pthash;

    test_pgm();

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