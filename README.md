<!-- [![CodeQL](https://github.com/jermp/pthash/actions/workflows/codeql.yml/badge.svg)](https://github.com/jermp/pthash/actions/workflows/codeql.yml) -->


Introduction
----
LMPTHASH is an index structure developed by the group for managing PB-level data mapping tables in DPU. It combines learning-based indexing, perfect hashing tables. By leveraging the continuity of logical addresses, it improves the spatial efficiency and query speed of the index.


## LMPTHASH Library 
Run the following commnads to build the LMPTHASH library: 
```
mkdir build
cd build
cmake ..
make
```
Library file `libclmpthash.a` will be outputed in `build`


### Usage Example

The following is an example of how to use the LMPTHASH library, which can also be found at `examples/clmpthash_test.c`:
```
#include "../include/clmpthash.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <stdatomic.h>

int main(){
    char* config_path="random_config";  // config file path
    clmpthash_config cfg;
    clmpthash_lva* uniq_lpns;
    clmpthash_physical_addr* ppns;
    clmpthash_lva* lpns;
    uint64_t num_uniq_lpn, num_total_lpn=100000;    // randomly generate 100000 LPNs

    // generate random LPNs, make sure all lpns are unique and sorted in uniq_lpns
    // first 8 bytes of each ppn is set the value of the corresponding lpn
    clmpthash_parse_configuration(config_path, &cfg, &uniq_lpns, &ppns, &num_uniq_lpn, &lpns, &num_total_lpn);

    // build LMPTHASH index
    void* index = clmpthash_build_index(uniq_lpns, ppns, num_uniq_lpn, &cfg);
    if (index == NULL) {
        printf("error building index\n");
        return -1;
    }

    int all_passed = 1;
    clmpthash_physical_addr pa;
    for (uint64_t i = 0; i < num_total_lpn; ++i) {
        int ret = clmpthash_get_pa(lpns[i], index, &pa);

        // check if the result is correct
        if (ret == 0) {
            clmpthash_lva _lva = 0;
            for (int j = 0; j < 8; j++) { _lva = (_lva << 8) + pa.data[j]; }
            if (_lva != lpns[i]) {
                printf("wrong result: should be 0x%lx, but got 0x%lx\n", lpns[i], _lva);
                all_passed = 0;
                break;
            }
        } else {
            printf("error getting pa\n");
            all_passed = 0;
            break;
        }
    }
    if (all_passed) {
        printf("all tests passed\n");
    }
    clmpthash_clean_index(index);
    clmpthash_clean_bufs(uniq_lpns, ppns, lpns);
    return 0;
}
```

Build the example and link the library:
```
cd example
gcc -o ctest clmpthash_test.c -L../build -lclmpthash -lstdc++ -lm -lpthread
```

## HiDPU Simulator

Based on LMPTHASH index structure, we build `HiDPU` which offload the inner index to the DPU and use the techniques such as parallel memory access and caching to improve index lookup efficiency. 

Since DPU code rely on specific microcode to call the hardware accelerator, we designed and provide the `HiDPU Simulator` to simulate the behavior of HiDPU in a CPU-only environment.

Run `HiDPU` requires the config file, for example: 
```
alpha 0.5
beta 0.5
gamma 0.5
P 1000
hashed_num_bucket_c 6.0
table_size_alpha 0.99
max_bucket_size 10
pilot_search_threshold 65536
dynamic_alpha 1
alpha_limits 0.8
left_epsilon 1
right_epsilon 100
trace_type msr
trace_path /home/quieoo/ftl/trace/MSR-Cambridge/web_2.csv
```
Among those the parameters, `trace_type` and `trace_path` describe the trace file used to create index and run the queries. Trace files can be downloaded from [here](http://iotta.snia.org/traces/block-io/388).

Currently, we support the following trace types:
- `msr`: MSR trace
- `random`: randomly generated trace
- `femu`: the trace captured by running FEMU and recording the LPNs accessed by the guest OS. Each line in the femu trace looks like this: ```"data time lpn"```

Meaning of other parameters can be found in [here](include/clmpthash.h).


Build the simulator:
```
cd example
mkdir build
cmake ..
make
```


HiDPU also integrates two baseline index structures: `Three-level Page Table` and `Learned Table`.
### Three-level Page Table

Run the following command to build and test three-level page table on given trace file:
```
./build/hidpu pagetable <config_path>
```

### Learned Table
Run the following command to build and test three-level page table on given trace file:
```
./build/hidpu learnedtable <config_path>
```

### HiDPU with Reconstruction and Multi-threads
Run the following command to build and test HiDPU on given trace file:
```
./build/hidpu lmpthash <num_threads> <if_reconstruction> <config_path>
```

Among these parameters, `num_threads` is the number of concurrent threads to run the queries, `if_reconstruction` is a flag to indicate whether to do `LocalReconstruction` during the query, if so the value should be `1`, otherwise `0`.

### Scalability Test
In this test, we extend the original trace with a given `scale_factor` while maintaineing the same access pattern. We use the extended trace to build the index, and report the index size.
```
./build/hidpu scalability <scale_factor> <config_path>
```

### Comparison with LearnedFTL
We've also implemented the [LearnedFTL](https://github.com/astlxmu/LearnedFTL) into HiDPU Simulator for comparison. Run the following to build and test LearnedFTL on given trace file:
```
./build/hidpu learnedftl <num_threads> <config_path>
```

For LearnedFTL, the SSD capacity is configured as 1TB with a 4KB page size. Other parameters are set according to the specifications in the paper (e.g., the total index size, including the bitmap, CMT, and learned models, is set to 3% of the mapping table size).

Please note that the original LearnedFTL implementation does not support multi-threading; therefore, the num_threads parameter must be set to 1.
