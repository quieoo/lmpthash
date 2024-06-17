#ifndef clmpthash_h
#define clmpthash_h

#ifdef __cplusplus
extern "C" {
#endif

#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>


typedef uint64_t clmpthash_lva;
typedef struct clmpthash_physical_addr{
    uint8_t data[20];
} clmpthash_physical_addr;

typedef struct clmpthash_lva2pa{
    clmpthash_physical_addr pa;
    clmpthash_lva lva;
}clmpthash_lva2pa;

typedef struct clmpthash_config{
    /*segmentation parameters*/
    // 'alpha' determines the impact of the total number of valid keys within a segment on priority during segmentation
    float alpha;
    // 'beta' determines the impact of the distance between two segments, i.e., the total number of empty keys after merging, on priority during segmentation.
    float beta;
    // 'gamma' determines the preference for maintaining a completely contiguous segment during segmentation, i.e., avoiding its merger with other segments as much as possible, thus creating a Accurate Segment
    float gamma;
    // P determines the number of segments
    uint32_t P;

    /*bucketing parameters*/
    // 'hashed_num_bucket_c' determines the number of buckets when building a hash-mapping based pthash
    double hashed_num_bucket_c;
    // The table size of PThash equals num_keys/'alpha'. When dynamic alpha configuration is adopted, 'alpha' is multiplied by itself after each iteration timeout, thus reducing the search difficulty
    double table_size_alpha;
    // When using linear mapping, the number of buckets is found by binary search between 1 and 'max_bucket_size' to identify the maximum value that meets the requirement. The requirement is to establish PThash after achieving a specified level of difficulty.
    uint32_t max_bucket_size;
    // The complexity threshold for searching the Pilot in each bucket, if exceeded without finding a suitable pilot, leads to an immediate exit and failure.
    uint32_t pilot_search_threshold;
    // If 'dynamic_alpha' is set to 1, the value of 'alpha' is multiplied by itself after each iteration timeout
    bool dynamic_alpha;
    // When employing dynamic alpha, it ensures that alpha does not fall below 'alpha_limit'. If it does, it regresses to hash-mapping
    double alpha_limits;

    /* Learned Index parameters */
    // Binary search to find the minimum epsilon that minimizes the number of layers.
    int left_epsilon;
    int right_epsilon;
}clmpthash_config;

typedef struct clmpthash_pgm_segment{
    uint64_t key;
    uint32_t slope;
    uint32_t intercept;
}clmpthash_pgm_segment;
static uint32_t uintfloat_mask = ((uint32_t)1 << 31);

typedef struct clmpthash_htl_segment{
    uint64_t addr;
    uint64_t meta;
}clmpthash_htl_segment;



void clmpthash_parse_configuration(char* config_path, clmpthash_config* cfg, clmpthash_lva** lvas, clmpthash_physical_addr** pas, uint64_t* num_lva, clmpthash_lva** querys, uint64_t* num_querys);
void clmpthash_clean_bufs(clmpthash_lva* lvas, clmpthash_physical_addr* pas, clmpthash_lva* querys);


void* clmpthash_build_index(clmpthash_lva* lvas, clmpthash_physical_addr* pas, uint64_t num, clmpthash_config* cfg);
int clmpthash_get_pa(clmpthash_lva lva, void* index, clmpthash_physical_addr* pa);
int clmpthash_clean_index(void* index);


void* clmpthash_offload_index(void* index);
int clmpthash_clean_offloaded_index(void* inner_index);

void* clt_build_index(clmpthash_lva* lvas, clmpthash_physical_addr* pas, uint64_t num, clmpthash_config* cfg);

#ifdef __cplusplus
}
#endif


#endif