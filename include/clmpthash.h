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


typedef uint64_t LVA;
typedef struct PhysicalAddr{
    uint8_t data[20];
} PhysicalAddr;

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

void parse_configuration(char* config_path, clmpthash_config* cfg, LVA** lvas, PhysicalAddr** pas, uint64_t* num_lva, LVA** querys, uint64_t* num_querys);

void* build_index(LVA* lvas, PhysicalAddr* pas, uint64_t num, clmpthash_config* cfg);
int get_pa(LVA lva, void* index, PhysicalAddr* pa);

#ifdef __cplusplus
}
#endif


#endif