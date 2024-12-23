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