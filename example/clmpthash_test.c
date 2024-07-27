#include "../include/clmpthash.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <pthread.h>
#include <stdatomic.h>

#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))

int nof_dpu_offload_index(void* index);

int test_host_side_clmpthash(char* config) {
    clmpthash_config cfg;
    clmpthash_lva* lvas;
    clmpthash_physical_addr* pas;
    clmpthash_lva* querys;
    uint64_t num_lva, num_querys;
    clmpthash_parse_configuration(config, &cfg, &lvas, &pas, &num_lva, &querys, &num_querys);

    void* index = clmpthash_build_index(lvas, pas, num_lva, &cfg);
    if (index == NULL) {
        printf("error building index\n");
        return -1;
    }

    clmpthash_physical_addr pa;
    for (uint64_t i = 0; i < num_querys; ++i) {
        if (i % 10000 == 0) {
            printf("\r    %lu / %lu\n", i, num_querys);
            printf("\033[1A");
        }

        int ret = clmpthash_get_pa(querys[i], index, &pa);

        // check if the result is correct: first 8 bytes should be equal to the original lva
        if (ret == 0) {
            clmpthash_lva _lva = 0;
            for (int j = 0; j < 8; j++) { _lva = (_lva << 8) + pa.data[j]; }
            if (_lva != querys[i]) {
                printf("wrong result: should be 0x%lx, but got 0x%lx\n", querys[i], _lva);
                return -1;
            }
        } else {
            printf("error getting pa\n");
            return -1;
        }
        // break;
    }
    printf("\ntest passed\n");
    clmpthash_clean_index(index);
    clmpthash_clean_bufs(lvas, pas, querys);
    return 0;
}

uint64_t MurmurHash2_64(void const* key, size_t len, uint64_t seed) {
    const uint64_t m = 0xc6a4a7935bd1e995ULL;
    const int r = 47;

    uint64_t h = seed ^ (len * m);

#if defined(__arm) || defined(__arm__)
    const size_t ksize = sizeof(uint64_t);
    const unsigned char* data = (const unsigned char*)key;
    const unsigned char* end = data + (std::size_t)(len / 8) * ksize;
#else
    const uint64_t* data = (const uint64_t*)key;
    const uint64_t* end = data + (len / 8);
#endif

    while (data != end) {
#if defined(__arm) || defined(__arm__)
        uint64_t k;
        memcpy(&k, data, ksize);
        data += ksize;
#else
        uint64_t k = *data++;
#endif

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    const unsigned char* data2 = (const unsigned char*)data;

    switch (len & 7) {
        // fall through
        case 7:
            h ^= (uint64_t)(data2[6]) << 48;
        // fall through
        case 6:
            h ^= (uint64_t)(data2[5]) << 40;
        // fall through
        case 5:
            h ^= (uint64_t)(data2[4]) << 32;
        // fall through
        case 4:
            h ^= (uint64_t)(data2[3]) << 24;
        // fall through
        case 3:
            h ^= (uint64_t)(data2[2]) << 16;
        // fall through
        case 2:
            h ^= (uint64_t)(data2[1]) << 8;
        // fall through
        case 1:
            h ^= (uint64_t)(data2[0]);
            h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

typedef struct bucket_cache_entry {
    uint64_t key;
    uint64_t value;
} bucket_cache_entry;
#define CACHE_SIZE (1024 * 1024 / 16)
#define CACHE_BUCKET_SIZE 4
typedef struct bucket_cache {
    uint32_t table_size;
    bucket_cache_entry* data;

    uint64_t total_access;
    uint64_t hit_count;

    atomic_int cnt[CACHE_SIZE/CACHE_BUCKET_SIZE];   // one counter for each four cache entry
} bucket_cache;

bucket_cache* bucket_cache_init(uint32_t table_size) {
    bucket_cache* bc = (bucket_cache*)malloc(sizeof(bucket_cache));
    memset(bc, 0, sizeof(bucket_cache));
    bc->table_size = table_size;
    bc->data = malloc(table_size * sizeof(bucket_cache_entry));
    memset(bc->data, 0xFF, table_size * sizeof(bucket_cache_entry));  // cache initialized to 0xFF, in case of key=0
    for(int i=0;i<CACHE_SIZE/CACHE_BUCKET_SIZE;i++){
        atomic_init(&bc->cnt[i], 0);
    }
    return bc;
}


uint64_t bucket_cache_get(bucket_cache* bc, uint64_t seg_first_key, uint32_t compressed_blk, uint8_t version) {
    if (compressed_blk >= (1 << 24)) {
        printf("error compressed_blk too large: %d\n", compressed_blk);
        return UINT64_MAX;
    }

    uint64_t key = (((uint64_t)seg_first_key << 24) & 0x00ffffffff000000) | (((uint64_t)compressed_blk) & 0xffffff);
    uint32_t idx = MurmurHash2_64(&key, sizeof(key), 0) % (CACHE_SIZE / CACHE_BUCKET_SIZE);

    bc->total_access++;

    for (int i = 0; i < CACHE_BUCKET_SIZE; i++) {
        bucket_cache_entry enty = bc->data[CACHE_BUCKET_SIZE * idx + i];
        uint64_t key2 = enty.key >> 8;
        uint64_t version2 = enty.key & 0xff;

        if (key2 == key && version2 == version) {
            bc->hit_count++;
            return enty.value;
        }else if(key2==key) {
            // printf("version mismatch: %d, %d\n", version, version2);
        }
    }

    return UINT64_MAX;
}

void bucket_cache_put(bucket_cache* bc, uint32_t seg_id, uint32_t compressed_blk, uint64_t value, uint8_t version) {
    uint64_t key = (((uint64_t)seg_id << 24) & 0xffffffff000000) | (((uint64_t)compressed_blk) & 0xffffff);
    uint32_t idx = MurmurHash2_64(&key, sizeof(key), 0) % (CACHE_SIZE / CACHE_BUCKET_SIZE);
    int s_idx = atomic_fetch_add(&bc->cnt[idx], 1) % CACHE_BUCKET_SIZE;
    // printf("idx: %u, s_idx: %u\n", idx, s_idx);

    bucket_cache_entry enty;
    enty.key = ((key << 8) & 0xffffffffffffff00) | ((uint64_t)version & 0xff);
    enty.value = value;
    bc->data[idx * CACHE_BUCKET_SIZE + s_idx] = enty;
}


void bucket_cache_clean(bucket_cache* bc) {
    free(bc->data);
    free(bc);
}

typedef struct Lindex_metadata {
    uint16_t num_levels;
    uint16_t epsilon;
    uint16_t level_offsets[6];
} Lindex_metadata;

typedef struct pgm_segment {
    uint64_t key;
    uint32_t slope;
    uint32_t intercept;
} pgm_segment;

typedef struct pgm_first_key {
    uint64_t key1;
    uint64_t key2;
} pgm_first_key;

typedef struct spram {
    pgm_segment data[4];
} spram;

spram global_spram[64];

#define lmpthash_tmp_spram ((uint8_t*)(&(global_spram[30])))
#define lmpthash_lindex_meta_spram ((uint8_t*)(&(global_spram[32])))

// #define lmpthash_lindex_segment_spram1 ((uint8_t*)(&(global_spram[34])))
// #define lmpthash_lindex_segment_spram2 ((uint8_t*)(&(global_spram[36])))
// #define lmpthash_lindex_segment_spram3 ((uint8_t*)(&(global_spram[38])))
// #define lmpthash_lindex_segment_spram4 ((uint8_t*)(&(global_spram[40])))

// #define lmpthash_htl_first_key_spram1 ((uint8_t*)(&(global_spram[42])))
// #define lmpthash_htl_first_key_spram2 ((uint8_t*)(&(global_spram[44])))
// #define lmpthash_htl_first_key_spram3 ((uint8_t*)(&(global_spram[46])))
// #define lmpthash_htl_first_key_spram4 ((uint8_t*)(&(global_spram[48])))

#define lmpthash_cache_segment_spram ((uint8_t*)(&(global_spram[50])))
#define lmphash_cache_segment_get_spram ((uint8_t*)(&(global_spram[52])))

#define l1_segs(id) ((clmpthash_pgm_segment*)(&(global_spram[34 + ((id) << 1)])))
#define htl_fks(id) ((uint64_t*)(&(global_spram[42 + ((id) & ~1)])) + ((id) & 1))


#define lt_segs(id) ((clmpthash_pgm_segment*)(&(global_spram[34 + ((id) << 1)])))
#define lva2pas_spram ((uint8_t*)(&(global_spram[52])))
#define lva2pa_lva(id)  ((clmpthash_lva*)(lva2pas_spram + id*28))
#define lva2pa_pa(id)  ((clmpthash_physical_addr*)(lva2pas_spram + id*28 + 8))

void load_16B_with_tmp_spram(void* dst, void* inner_index, uint32_t idx){
    memcpy(lmpthash_tmp_spram, inner_index+idx*16, 16);
    uint64_t* d=(uint64_t*)dst;
    uint64_t* s=(uint64_t*)lmpthash_tmp_spram;
    *d=*s;
    *(d+1)=*(s+1);
}

uint64_t get_u64_from_sm(void* sm, uint32_t idx, uint8_t oft) {
    // printf("idx=%d, oft=%d\n", idx, oft);
    clmpthash_htl_segment* sm_16B=(clmpthash_htl_segment*)sm;
    clmpthash_htl_segment tmp = sm_16B[idx];    // load to spram

    uint64_t* p=(uint64_t*)(&tmp);

    return p[oft];
}


int test_lt(char* config) {
    clmpthash_config cfg;
    clmpthash_lva* lvas;
    clmpthash_physical_addr* pas;
    clmpthash_lva* querys;
    uint64_t num_lva;
    uint64_t num_querys;
    clmpthash_parse_configuration(config, &cfg, &lvas, &pas, &num_lva, &querys, &num_querys);

    void* inner_index = clt_build_index(lvas, pas, num_lva, &cfg);
    if (inner_index == NULL) {
        printf("error building index\n");
        return -1;
    }
    printf("build index done\n");

    clmpthash_physical_addr pa1;
    uint8_t cache_last_htl_segment=1;
    uint8_t cache_last_htl_segment_hit=0;

    load_16B_with_tmp_spram(lmpthash_lindex_meta_spram, inner_index, 1);
    Lindex_metadata* lmd_ = (Lindex_metadata*)lmpthash_lindex_meta_spram;
    uint16_t sidx_level0 = lmd_->level_offsets[lmd_->num_levels - 1];
    load_16B_with_tmp_spram(lmpthash_lindex_meta_spram+16, inner_index, sidx_level0);
    load_16B_with_tmp_spram(lmpthash_lindex_meta_spram+32, inner_index, sidx_level0+1);

    // get Lindex_metadata
    uint16_t l = lmd_->num_levels;
    uint32_t epsilon = (lmd_->epsilon&0xff) +2;
    uint32_t er=((lmd_->epsilon)>>8) + 1;
    // printf("l: %u, epsilon: %u, er: %u\n", l, epsilon, er);
    uint16_t* level_offsets = lmd_->level_offsets;

    uint8_t i;
    
    uint64_t spram_cache_hit=0;
    uint64_t average_dma_size=epsilon*2*28;
    uint64_t dma_num=0;
    for (uint64_t q = 0; q < num_querys; ++q) {
        // printf("query %lu\n", i);
        if (q % 10000 == 0) {
            printf("\r    %lu / %lu\n", q, num_querys);
            printf("\033[1A");
        }

        clmpthash_lva lva = querys[q];
        

        // cached last htl segment
        clmpthash_pgm_segment* pgm_seg = (clmpthash_pgm_segment*)lmpthash_cache_segment_spram;
        uint64_t s_first_key;
        cache_last_htl_segment_hit=0;
        uint32_t s_idx1;
        uint32_t s_idx2;
        int64_t pos;


        // printf("1: lva %lx\n", lva);
        if(cache_last_htl_segment){
            // cached_first_key <= lva
            if(pgm_seg[0].key <= lva){
                // printf("2: cached_firsy_key %lx\n", pgm_seg[0].key);
            
                // cached_next_key > lva
                // uint64_t* cached_next_key= htl_fks(cached_firsy_key[1]);
                if(pgm_seg[1].key > lva){
                    // printf("3: cached_next_key %lx\n", pgm_seg[1].key);
                    pos=(int64_t)(pgm_seg->slope) * (lva - pgm_seg->key) / ((uint32_t)1 << 31) +pgm_seg->intercept;
                    cache_last_htl_segment_hit=1;
                    spram_cache_hit++;
                }
            }
        }
        if(!cache_last_htl_segment || !cache_last_htl_segment_hit){
            
            // uint32_t num_htl_segment=level_offsets[l+2]-level_offsets[l+1];
            clmpthash_pgm_segment* level0_segments = (clmpthash_pgm_segment*)(lmpthash_lindex_meta_spram+16);

            // predict with cached level[l] segment
            pos =(int64_t)(level0_segments->slope) * (lva - level0_segments->key) / ((uint32_t)1 << 31) +level0_segments->intercept;
            // printf("key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva,level0_segments->key, level0_segments->slope, level0_segments->intercept, pos);
            if (pos < 0) pos = 0;
            if (pos > level0_segments[1].intercept) pos = level0_segments[1].intercept;
            l=lmd_->num_levels;
            while (l > 1) {
                // search and compare segment first key at level[l-1]
                s_idx2 = level_offsets[l - 2] + PGM_SUB_EPS(pos, er);  // searching start point at next level
                // printf("local search: %u : %u\n", s_idx2, s_idx2+er*2);
                uint32_t s_idx1 = s_idx2;
                while (true) {
                    // printf("active_spram: %u\n", active_spram);
                    
                    // load 4 segments at a time
                    // i=0;
                    // for(;i<4;i++)
                    //     memcpy(l1_segs(i), inner_index + (s_idx2+1+i)*16, 16);    
                    memcpy(lt_segs(0), inner_index+(++s_idx2)*16, 16);
                    memcpy(lt_segs(1), inner_index+(++s_idx2)*16, 16);
                    memcpy(lt_segs(2), inner_index+(++s_idx2)*16, 16);
                    memcpy(lt_segs(3), inner_index+(++s_idx2)*16, 16);
                    memcpy(lt_segs(4), inner_index+(++s_idx2)*16, 16);
                    memcpy(lt_segs(5), inner_index+(++s_idx2)*16, 16);
                    memcpy(lt_segs(6), inner_index+(++s_idx2)*16, 16);
                    memcpy(lt_segs(7), inner_index+(++s_idx2)*16, 16);

                    // lt_segs_0_outstand()
                    if(lt_segs(0)->key > lva){
                        if(s_idx1+8==s_idx2){
                            memcpy(lt_segs(7)+1, inner_index + (s_idx1) * 16, 16);    
                        }
                        pos = (int64_t)(lt_segs(7)[1].slope) * (lva - lt_segs(7)[1].key) / ((uint32_t)1 << 31) +  lt_segs(7)[1].intercept;
                        // printf("[i==0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, lt_segs(3)[1].key, lt_segs(3)[1].slope, lt_segs(3)[1].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > lt_segs(0)[0].intercept)pos = lt_segs(0)[0].intercept;
                        pgm_seg[0].key=lt_segs(7)[1].key;
                        pgm_seg[0].slope=lt_segs(7)[1].slope;
                        pgm_seg[0].intercept=lt_segs(7)[1].intercept;
                        pgm_seg[1].key=lt_segs(0)[0].key;
                        break;
                    }
                    //lt_segs_1_outstand()
                    if(lt_segs(1)->key > lva){
                        pos = (int64_t)(lt_segs(0)[0].slope) * (lva - lt_segs(0)[0].key) / ((uint32_t)1 << 31) +  lt_segs(0)[0].intercept;
                        // printf("[i>0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, lt_segs(i-1)[0].key, lt_segs(i-1)[0].slope, lt_segs(i-1)[0].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > lt_segs(1)[0].intercept)pos = lt_segs(1)[0].intercept;
                        pgm_seg[0].key=lt_segs(0)[0].key;
                        pgm_seg[0].slope=lt_segs(0)[0].slope;
                        pgm_seg[0].intercept=lt_segs(0)[0].intercept;
                        pgm_seg[1].key=lt_segs(1)[0].key;
                        break;
                    }
                    //lt_segs_2_outstand()
                    if(lt_segs(2)->key > lva){
                        pos = (int64_t)(lt_segs(1)[0].slope) * (lva - lt_segs(1)[0].key) / ((uint32_t)1 << 31) +  lt_segs(1)[0].intercept;
                        // printf("[i>0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, lt_segs(i-1)[0].key, lt_segs(i-1)[0].slope, lt_segs(i-1)[0].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > lt_segs(2)[0].intercept)pos = lt_segs(2)[0].intercept;
                        pgm_seg[0].key=lt_segs(1)[0].key;
                        pgm_seg[0].slope=lt_segs(1)[0].slope;
                        pgm_seg[0].intercept=lt_segs(1)[0].intercept;
                        pgm_seg[1].key=lt_segs(2)[0].key;
                        break;
                    }
                    //lt_segs_3_outstand()
                    if(lt_segs(3)->key > lva){
                        pos = (int64_t)(lt_segs(2)[0].slope) * (lva - lt_segs(2)[0].key) / ((uint32_t)1 << 31) +  lt_segs(2)[0].intercept;
                        // printf("[i>0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, lt_segs(i-1)[0].key, lt_segs(i-1)[0].slope, lt_segs(i-1)[0].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > lt_segs(3)[0].intercept)pos = lt_segs(3)[0].intercept;
                        pgm_seg[0].key=lt_segs(2)[0].key;
                        pgm_seg[0].slope=lt_segs(2)[0].slope;
                        pgm_seg[0].intercept=lt_segs(2)[0].intercept;
                        pgm_seg[1].key=lt_segs(3)[0].key;
                        break;
                    }
                    //lt_segs_4_outstand()
                    if(lt_segs(4)->key > lva){
                        pos = (int64_t)(lt_segs(3)[0].slope) * (lva - lt_segs(3)[0].key) / ((uint32_t)1 << 31) +  lt_segs(3)[0].intercept;
                        // printf("[i>0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, lt_segs(i-1)[0].key, lt_segs(i-1)[0].slope, lt_segs(i-1)[0].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > lt_segs(4)[0].intercept)pos = lt_segs(4)[0].intercept;
                        pgm_seg[0].key=lt_segs(3)[0].key;
                        pgm_seg[0].slope=lt_segs(3)[0].slope;
                        pgm_seg[0].intercept=lt_segs(3)[0].intercept;
                        pgm_seg[1].key=lt_segs(4)[0].key;
                        break;
                    }
                    //lt_segs_5_outstand()
                    if(lt_segs(5)->key > lva){
                        pos = (int64_t)(lt_segs(4)[0].slope) * (lva - lt_segs(4)[0].key) / ((uint32_t)1 << 31) +  lt_segs(4)[0].intercept;
                        // printf("[i>0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, lt_segs(i-1)[0].key, lt_segs(i-1)[0].slope, lt_segs(i-1)[0].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > lt_segs(5)[0].intercept)pos = lt_segs(5)[0].intercept;
                        pgm_seg[0].key=lt_segs(4)[0].key;
                        pgm_seg[0].slope=lt_segs(4)[0].slope;
                        pgm_seg[0].intercept=lt_segs(4)[0].intercept;
                        pgm_seg[1].key=lt_segs(5)[0].key;
                        break;
                    }
                    //lt_segs_6_outstand()
                    if(lt_segs(6)->key > lva){
                        pos = (int64_t)(lt_segs(5)[0].slope) * (lva - lt_segs(5)[0].key) / ((uint32_t)1 << 31) +  lt_segs(5)[0].intercept;
                        // printf("[i>0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, lt_segs(i-1)[0].key, lt_segs(i-1)[0].slope, lt_segs(i-1)[0].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > lt_segs(6)[0].intercept)pos = lt_segs(6)[0].intercept;
                        pgm_seg[0].key=lt_segs(5)[0].key;
                        pgm_seg[0].slope=lt_segs(5)[0].slope;
                        pgm_seg[0].intercept=lt_segs(5)[0].intercept;
                        pgm_seg[1].key=lt_segs(6)[0].key;
                        break;
                    }
                    //lt_segs_7_outstand()
                    if(lt_segs(7)->key > lva){
                        pos = (int64_t)(lt_segs(6)[0].slope) * (lva - lt_segs(6)[0].key) / ((uint32_t)1 << 31) +  lt_segs(6)[0].intercept;
                        // printf("[i>0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, lt_segs(i-1)[0].key, lt_segs(i-1)[0].slope, lt_segs(i-1)[0].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > lt_segs(7)[0].intercept)pos = lt_segs(7)[0].intercept;
                        pgm_seg[0].key=lt_segs(6)[0].key;
                        pgm_seg[0].slope=lt_segs(6)[0].slope;
                        pgm_seg[0].intercept=lt_segs(6)[0].intercept;
                        pgm_seg[1].key=lt_segs(7)[0].key;
                        break;
                    }
                    // backup current lt_segs[3][0] to lt_segs[3][1]
                    lt_segs(7)[1]= lt_segs(7)[0];
                }
                l=l-1;
            }
        }
        // reach level-0, search for the right first_key
        
        // printf("pos: %lu\n", pos);
        s_idx1 = (PGM_SUB_EPS(pos, epsilon));
        s_idx2=s_idx1+epsilon*2;
        // get sub_table addr
        uint64_t* sub_table_addr=(uint64_t*)(inner_index+1024*1024+8);

        // uint16_t dma_width=640;
        // printf("bottom level, local search %u : %u \n", s_idx1, s_idx2);
        // if split to 2 bufs, 74897 lva2pas in each dma bufer
        if(s_idx1/74897 == s_idx2/74897){
            // dma with width=(epsilon*2)*28
            uint64_t sb_addr=get_u64_from_sm(sub_table_addr, s_idx1/74897/2, s_idx1/74897%2);
            sb_addr+=(s_idx1%74897)*28;
            uint64_t end_addr=sb_addr+(epsilon*2)*28;
            uint8_t found=0;
            while(sb_addr<end_addr){
                uint16_t dma_width=(uint16_t)(end_addr-sb_addr);
                if(dma_width>640) dma_width=640;
                // dma_width align to 28 bytes
                if(dma_width%28!=0) dma_width-=(dma_width%28);
                memcpy(lva2pas_spram, (void*)sb_addr, dma_width);
                dma_num++;
                sb_addr+=dma_width;
                i=0;
                for(;i<dma_width/28;++i){
                    if(*(lva2pa_lva(i)) == lva){
                        pa1=*(lva2pa_pa(i));
                        found=1;
                        break;
                    }
                }
                if(found)   break;
            }
            // memcpy(lva2pas_spram, (void*)(sb_addr+(s_idx1%74897)*28), epsilon*2*28);
            // dma_num++;
            // i=0;
            // for(;i<epsilon*2;++i){
            //     // printf("i: %u, lva: %lu\n",i, *(lva2pa_lva(i)));
            //     // printf("pa:\n");
            //     // for(int p=0;p<20;++p){
            //     //     printf("%02x ", lva2pa_pa(i)[0].data[p]);
            //     // }
            //     // printf("\n");
            //     if(*(lva2pa_lva(i)) == lva){
            //         pa1=*(lva2pa_pa(i));
            //         break;
            //     }
            // }
        }else{
            // dma 2 times untile find the right pa
            uint32_t num_lva2pa_first_batch=74897-(s_idx1%74897);
            uint64_t sb_addr=get_u64_from_sm(sub_table_addr, s_idx1/74897/2, s_idx1/74897%2);
            sb_addr+=(s_idx1%74897)*28;
            uint64_t end_addr=sb_addr+(num_lva2pa_first_batch)*28;
            uint8_t found=0;
            while(sb_addr<end_addr){
                uint16_t dma_width=(uint16_t)(end_addr-sb_addr);
                if(dma_width>640) dma_width=640;
                // dma_width align to 28 bytes
                if(dma_width%28!=0) dma_width-=(dma_width%28);
                memcpy(lva2pas_spram, (void*)sb_addr, dma_width);
                dma_num++;
                sb_addr+=dma_width;
                i=0;
                for(;i<(dma_width/28);++i){
                    if(*(lva2pa_lva(i)) == lva){
                        pa1=*(lva2pa_pa(i));
                        found=1;
                        break;
                    }
                }
                if(found)   break;
            }
            // memcpy(lva2pas_spram, (void*)(sb_addr+(s_idx1%74897)*28), num_lva2pa_first_batch*28);
            // dma_num++;
            // i=0;
            // for(;i<(num_lva2pa_first_batch);++i){
                
            //     // printf("i: %u, lva: %lu\n",i, *(lva2pa_lva(i)));
            //     // printf("pa:\n");
            //     // for(int p=0;p<20;++p){
            //     //     printf("%02x ", lva2pa_pa(i)[0].data[p]);
            //     // }
            //     // printf("\n");
            //     if(*(lva2pa_lva(i)) == lva){
            //         pa1=*(lva2pa_pa(i));
            //         break;
            //     }
            // }
            if(!found){
                sb_addr=get_u64_from_sm(sub_table_addr, s_idx2/74897/2, s_idx2/74897%2);
                end_addr=sb_addr+(s_idx2%74897)*28;
                while(sb_addr<end_addr){
                    uint16_t dma_width=(uint16_t)(end_addr-sb_addr);
                    if(dma_width>640) dma_width=640;
                    // dma_width align to 28 bytes
                    if(dma_width%28!=0) dma_width-=(dma_width%28);
                    memcpy(lva2pas_spram, (void*)sb_addr, dma_width);
                    dma_num++;
                    sb_addr+=dma_width;
                    i=0;
                    for(;i<(dma_width/28);++i){
                        if(*(lva2pa_lva(i)) == lva){
                            pa1=*(lva2pa_pa(i));
                            found=1;
                            break;
                        }
                    }
                    if(found)   break;
                }

                // memcpy(lva2pas_spram, (void*)sb_addr, (s_idx2%74897)*28);
                // i=0;
                // dma_num++;
                // for(;i<(s_idx2%74897);++i){
                    
                //     // printf("i: %u, lva: %lu\n",i, *(lva2pa_lva(i)));
                //     // printf("pa:\n");
                //     // for(int p=0;p<20;++p){
                //     //     printf("%02x ", lva2pa_pa(i)[0].data[p]);
                //     // }
                //     // printf("\n");
                //     if(*(lva2pa_lva(i)) == lva){
                //         pa1=*(lva2pa_pa(i));
                //         break;
                //     }
                // }
            }
        }
        clmpthash_lva _lva = 0;
        for (int j = 0; j < 8; j++) { _lva = (_lva << 8) + pa1.data[j]; }
        if (_lva != lva) {
            printf("%lu: wrong result, should be %lu, but got %lu\n", q, lva, _lva);
            return -1;
        }
    }
    printf("pass all queries\n");
    printf("spram cache hit ratio: %f\n",(double)spram_cache_hit/num_querys);
    printf("dma size: %lu, average dma num: %f\n", average_dma_size, (double)dma_num/num_querys);
}


#define t_lmpthash_tmp_spram(t_spram) ((uint8_t*)(&((t_spram)[30])))
#define t_lmpthash_lindex_meta_spram(t_spram) ((uint8_t*)(&((t_spram)[32])))
#define t_lmpthash_cache_segment_spram(t_spram) ((uint8_t*)(&((t_spram)[50])))
#define t_lmpthash_cache_segment_get_spram(t_spram) ((uint8_t*)(&((t_spram)[52])))
#define t_l1_segs(t_spram, id) ((clmpthash_pgm_segment*)(&((t_spram)[34 + ((id) << 1)])))
#define t_htl_fks(t_spram, id) ((uint64_t*)(&((t_spram)[42 + ((id) & ~1)])) + ((id) & 1))


typedef struct{
    uint64_t num_querys;
    clmpthash_lva* querys;
    int cache_last_htl_segment;
    void* inner_index;
    bucket_cache* bc;
    spram* thread_spram;
} dlmpht_thread_arg;

void dlmpht_query_func(void* arg){
    dlmpht_thread_arg* _arg=(dlmpht_thread_arg*)arg;
    uint64_t num_querys=_arg->num_querys;
    clmpthash_lva* querys=_arg->querys;
    int cache_last_htl_segment=_arg->cache_last_htl_segment;
    cache_last_htl_segment=0;
    void* inner_index=_arg->inner_index;
    bucket_cache* bc=_arg->bc;

    // initialize Lindex_meta and level0_segment
    load_16B_with_tmp_spram(t_lmpthash_lindex_meta_spram(_arg->thread_spram), inner_index, 1);
    Lindex_metadata* lmd_ = (Lindex_metadata*)t_lmpthash_lindex_meta_spram(_arg->thread_spram);
    uint16_t sidx_level0 = lmd_->level_offsets[lmd_->num_levels - 1];
    load_16B_with_tmp_spram(t_lmpthash_lindex_meta_spram(_arg->thread_spram)+16, inner_index, sidx_level0);
    load_16B_with_tmp_spram(t_lmpthash_lindex_meta_spram(_arg->thread_spram)+32, inner_index, sidx_level0+1);
    
    uint8_t i;
    clmpthash_physical_addr pa1;
    uint64_t cache_hit=0;
    uint64_t cache_total=0;
    // num_querys=1;

    for (uint64_t q = 0; q < num_querys; ++q) {
        // printf("query %lu\n", q);
        if (q % 10000 == 0) {
            // printf("\r    %lu / %lu\n", q, num_querys);
            // printf("\033[1A");

            // printf("    cache hit ratio: %f\n",(double)cache_hit/cache_total);
        }
        clmpthash_lva lva = querys[q];

        // query pgm index

        // cached last htl segment
        clmpthash_htl_segment* htl_seg = (clmpthash_htl_segment*)t_lmpthash_cache_segment_get_spram(_arg->thread_spram);
        uint64_t s_first_key;
        int cache_last_htl_segment_hit=0;
        // printf("1: lva %lx\n", lva);
        if(cache_last_htl_segment){
            // cached_first_key <= lva
            uint64_t* cached_firsy_key=(uint64_t*)(t_lmpthash_cache_segment_spram(_arg->thread_spram));
            if(*cached_firsy_key <= lva){
                // printf("2: cached_firsy_key %lx\n", *cached_firsy_key);
                // printf("fk: %lx, offset: %u, seg_offset: %u\n", cached_firsy_key[0], cached_firsy_key[1], cached_firsy_key[2]);
                // cached_next_key > lva
                // uint64_t* cached_next_key= htl_fks(cached_firsy_key[1]);
                if(*(t_htl_fks(_arg->thread_spram, cached_firsy_key[1])) > lva){
                    // printf("3: cached_next_key %lx\n", *(t_htl_fks(_arg->thread_spram,cached_firsy_key[1])));
                    s_first_key=*cached_firsy_key;
                    memcpy(htl_seg, inner_index + cached_firsy_key[2], 16);
                    cache_last_htl_segment_hit=1;
                }
            }
        }

        if(!cache_last_htl_segment || !cache_last_htl_segment_hit){
            // get Lindex_metadata
            Lindex_metadata* lmd = (Lindex_metadata*)t_lmpthash_lindex_meta_spram(_arg->thread_spram);
            uint16_t l = lmd->num_levels;
            uint32_t epsilon = lmd->epsilon+1;
            // printf("l: %u, epsilon: %u\n", l, epsilon);
            uint16_t* level_offsets = lmd->level_offsets;
            // uint32_t num_htl_segment=level_offsets[l+2]-level_offsets[l+1];
            clmpthash_pgm_segment* level0_segments = (clmpthash_pgm_segment*)(t_lmpthash_lindex_meta_spram(_arg->thread_spram)+16);
            
            // predict with cached level[l] segment
            int64_t pos =(int64_t)(level0_segments->slope) * (lva - level0_segments->key) / ((uint32_t)1 << 31) +level0_segments->intercept;
            // printf("key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva,level0_segments->key, level0_segments->slope, level0_segments->intercept, pos);
            if (pos < 0) pos = 0;
            if (pos > level0_segments[1].intercept) pos = level0_segments[1].intercept;

            // clmpthash_pgm_segment* l1_segs[4];
            // l1_segs[0] = (clmpthash_pgm_segment*)lmpthash_lindex_segment_spram1;
            // l1_segs[1] = (clmpthash_pgm_segment*)lmpthash_lindex_segment_spram2;
            // l1_segs[2] = (clmpthash_pgm_segment*)lmpthash_lindex_segment_spram3;
            // l1_segs[3] = (clmpthash_pgm_segment*)lmpthash_lindex_segment_spram4;
                      
            uint32_t s_idx1;
            uint32_t s_idx2;
            while (l > 1) {
                // search and compare segment first key at level[l-1]
                s_idx2 = level_offsets[l - 2] + PGM_SUB_EPS(pos, epsilon);  // searching start point at next level
                uint32_t s_idx1 = s_idx2;
                while (true) {
                    // printf("s_idx1: %u, s_idx2: %u\n", s_idx1, s_idx2);
                    // printf("active_spram: %u\n", active_spram);
                    
                    // load 4 segments at a time
                    // i=0;
                    // for(;i<4;i++)
                    //     memcpy(l1_segs(i), inner_index + (s_idx2+1+i)*16, 16);
                    
                    memcpy(t_l1_segs(_arg->thread_spram,0), inner_index+(++s_idx2)*16, 16);
                    memcpy(t_l1_segs(_arg->thread_spram,1), inner_index+(++s_idx2)*16, 16);
                    memcpy(t_l1_segs(_arg->thread_spram,2), inner_index+(++s_idx2)*16, 16);
                    memcpy(t_l1_segs(_arg->thread_spram,3), inner_index+(++s_idx2)*16, 16);
                    // for(int seg = 0; seg < 4; seg++){
                    //     printf("l1_segs[%d] key: %lu, slope: %u, intercept: %u\n", seg, l1_segs(seg)->key, l1_segs(seg)->slope, l1_segs(seg)->intercept);
                    // }

                    // l1_segs_0_outstand()
                    if(t_l1_segs(_arg->thread_spram,0)->key > lva){
                        if(s_idx1+4==s_idx2){
                            memcpy(t_l1_segs(_arg->thread_spram,3)+1, inner_index + (s_idx1) * 16, 16);    
                        }
                        pos = (int64_t)(t_l1_segs(_arg->thread_spram,3)[1].slope) * (lva - t_l1_segs(_arg->thread_spram,3)[1].key) / ((uint32_t)1 << 31) +  t_l1_segs(_arg->thread_spram,3)[1].intercept;
                        // printf("[i==0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, l1_segs(3)[1].key, l1_segs(3)[1].slope, l1_segs(3)[1].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > t_l1_segs(_arg->thread_spram,0)[0].intercept)pos = t_l1_segs(_arg->thread_spram,0)[0].intercept;
                        break;
                    }
                    //l1_segs_1_outstand()
                    if(t_l1_segs(_arg->thread_spram,1)->key > lva){
                        pos = (int64_t)(t_l1_segs(_arg->thread_spram,0)[0].slope) * (lva - t_l1_segs(_arg->thread_spram,0)[0].key) / ((uint32_t)1 << 31) +  t_l1_segs(_arg->thread_spram,0)[0].intercept;
                        // printf("[i==1] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, l1_segs(0)[0].key, l1_segs(0)[0].slope, l1_segs(0)[0].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > t_l1_segs(_arg->thread_spram,1)[0].intercept)pos = t_l1_segs(_arg->thread_spram,1)[0].intercept;
                        break;
                    }
                    //l1_segs_2_outstand()
                    if(t_l1_segs(_arg->thread_spram,2)->key > lva){
                        pos = (int64_t)(t_l1_segs(_arg->thread_spram,1)[0].slope) * (lva - t_l1_segs(_arg->thread_spram,1)[0].key) / ((uint32_t)1 << 31) +  t_l1_segs(_arg->thread_spram,1)[0].intercept;
                        // printf("[i==2] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, l1_segs(1)[0].key, l1_segs(1)[0].slope, l1_segs(1)[0].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > t_l1_segs(_arg->thread_spram,2)[0].intercept)pos = t_l1_segs(_arg->thread_spram,2)[0].intercept;
                        break;
                    }
                    //l1_segs_1_outstand()
                    if(t_l1_segs(_arg->thread_spram,3)->key > lva){
                        pos = (int64_t)(t_l1_segs(_arg->thread_spram,2)[0].slope) * (lva - t_l1_segs(_arg->thread_spram,2)[0].key) / ((uint32_t)1 << 31) +  t_l1_segs(_arg->thread_spram,2)[0].intercept;
                        // printf("[i==3] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, l1_segs(2)[0].key, l1_segs(2)[0].slope, l1_segs(2)[0].intercept, pos);
                        if (pos < 0) pos = 0;
                        if (pos > t_l1_segs(_arg->thread_spram,3)[0].intercept)pos = t_l1_segs(_arg->thread_spram,3)[0].intercept;
                        break;
                    }
                    // backup current l1_segs[3][0] to l1_segs[3][1]
                    t_l1_segs(_arg->thread_spram,3)[1]= t_l1_segs(_arg->thread_spram,3)[0];
                }
                l=l-1;
            }

            // reach level-0, search for the right first_key
            
            // align s_idx2 to 2
            s_idx1 = PGM_SUB_EPS(pos, epsilon);
            // printf("bottom level, predict sidx: %u\n", s_idx1);
            if (s_idx1 == 0) {
                s_idx1 = -1;
                s_idx2 = level_offsets[lmd->num_levels];
            } else {
                s_idx1 = (s_idx1 & 1) ? s_idx1 : s_idx1 - 1;
                s_idx2 = level_offsets[lmd->num_levels] + s_idx1 / 2 + 1;
            }
            do{
                memcpy(t_htl_fks(_arg->thread_spram,0), inner_index+(s_idx2++)*16, 16);
                memcpy(t_htl_fks(_arg->thread_spram,2), inner_index+(s_idx2++)*16, 16);
                memcpy(t_htl_fks(_arg->thread_spram,4), inner_index+(s_idx2++)*16, 16);
                memcpy(t_htl_fks(_arg->thread_spram,6), inner_index+(s_idx2++)*16, 16);
                // for(int k=0;k<8;k++) {
                //     printf("first_key: %lu\n", htl_fks(k)[0]);
                // }
                i=0;
                // htl_fks_0_outstand()
                if(t_htl_fks(_arg->thread_spram,0)[0]>lva) break;
                ++i;

                if(t_htl_fks(_arg->thread_spram,1)[0]>lva) break;
                ++i;

                // htl_fks_1_outstand()
                if(t_htl_fks(_arg->thread_spram,2)[0]>lva) break;
                ++i;

                if(t_htl_fks(_arg->thread_spram,3)[0]>lva) break;
                ++i;

                // htl_fks_2_outstand()
                if(t_htl_fks(_arg->thread_spram,4)[0]>lva) break;
                ++i;

                if(t_htl_fks(_arg->thread_spram,5)[0]>lva) break;
                ++i;

                // htl_fks_3_outstand()
                if(t_htl_fks(_arg->thread_spram,6)[0]>lva) break;
                ++i;

                if(t_htl_fks(_arg->thread_spram,7)[0]>lva) break;
                ++i;

                // can not find next_first_key > lva
                t_htl_fks(_arg->thread_spram,7)[1]= t_htl_fks(_arg->thread_spram,7)[0];   // backup current first_key
                s_idx1+=8;
            }while(1);
            // load htl segment
            if(i==0){
                if(s_idx2 == (level_offsets[lmd->num_levels] + s_idx1 / 2 + 1)+4){
                    // first round found the right first_key on htl_fks(0)
                    memcpy(t_htl_fks(_arg->thread_spram,6), inner_index+(s_idx2-5)*16, 16);
                    t_htl_fks(_arg->thread_spram,7)[1]= t_htl_fks(_arg->thread_spram,7)[0];
                }
                s_first_key=t_htl_fks(_arg->thread_spram,7)[1];
            }else{
                s_first_key=t_htl_fks(_arg->thread_spram,i-1)[0];
            }
            s_idx1+=i;

            memcpy(htl_seg, inner_index + (level_offsets[lmd->num_levels + 1] + s_idx1) * 16, 16);
            // save current first_key to cache
            if(cache_last_htl_segment){
                uint64_t* cached_p=(uint64_t*)t_lmpthash_cache_segment_spram(_arg->thread_spram);
                cached_p[0]=s_first_key;
                cached_p[1]=i;
                cached_p[2]=(level_offsets[lmd->num_levels + 1] + s_idx1) * 16;
            }
            
            // printf("sidx_1: %x, s_first_key: %lx\n", s_idx1, s_first_key);
        }
        
        // continue;
        uint8_t seg_type =( htl_seg->meta >> 62)&0x3;
        // printf("seg_type: %lx\n", htl_seg->meta);
        if (seg_type == 0) {
            // printf(" accurate segment\n");
            clmpthash_physical_addr* pas = (clmpthash_physical_addr*)(htl_seg->addr + 8);
            uint64_t pos = lva - s_first_key;
            pa1 = pas[pos];
        } else if (seg_type == 4) {
            printf(" error segment\n");
            exit(0);
        } else {
            // printf(" approximate segment\n");
            uint64_t hash = MurmurHash2_64(&lva, sizeof(clmpthash_lva), 0x123456789);
            uint64_t blk;
            uint32_t table_size = htl_seg->meta & 0xffffff;
            uint64_t table_addr = htl_seg->addr + 8;
            uint64_t bucket_addr = table_addr + table_size * sizeof(clmpthash_physical_addr);
            uint32_t slope = (htl_seg->meta >> 32) & 0x3fffff;
            uint64_t width = (htl_seg->meta >> 24)&0xff;
            uint8_t version=(htl_seg->meta >> 54)&0xff;
            if (seg_type == 1) {
                // linear bucketing
                blk = (lva - s_first_key) / slope;

            } else {
                // hash bucketing
                // printf("hash:%lu, slope: %lu\n", hash, slope);
                blk = hash % slope;
            }
            // printf("blk: %lu\n", blk);
            blk = blk * width;
            width = -(width == 64) | ((((uint64_t)1) << width) - 1);
            /*
            bucket_addr+=blk>>3;
            uint64_t p=*(uint64_t*)bucket_addr;
            */
            uint64_t p;
            uint32_t compressed_blk_id = blk >> 3;
            // uint64_t blk_value=bucket_cache_get(bc, sidx, compressed_blk_id);
            uint64_t blk_value = bucket_cache_get(bc, s_first_key, compressed_blk_id, version);
            // printf("lva: %lx, blk: %lx, blk_value: %lx, seg_type: %u\n", lva, blk, blk_value,
            // seg_type);
            cache_total++;
            // blk_value = UINT64_MAX;
            if (blk_value != UINT64_MAX) {
                // cache hit
                cache_hit++;
                p = blk_value;
            } else {
                bucket_addr += blk >> 3;
                p = *(uint64_t*)bucket_addr;
                // cache miss
                bucket_cache_put(bc, s_first_key, compressed_blk_id, p, version);
            }
            p = p >> (blk & 7) & width;
            // printf("p: %lu\n", p);
            p = MurmurHash2_64(&p, sizeof(uint64_t), 0x123456789);
            // printf("hash: %lu, hash_p: %lu\n", hash, p);

            p = (p ^ hash) % table_size;
            // printf("pos: %lx\n", p);
            clmpthash_physical_addr* pas = (clmpthash_physical_addr*)table_addr;
            pa1 = pas[p];
        }

        // clmpthash_lva _lva = 0;
        // for (int j = 0; j < 8; j++) { _lva = (_lva << 8) + pa1.data[j]; }
        // if (_lva != lva) {
        //     printf("%lu: wrong result, should be 0x%lx, but got 0x%lx\n", q, lva, _lva);
        //     return -1;
        // }
    }

    // output cache counter
    // for (int i = 0; i < CACHE_SIZE / CACHE_BUCKET_SIZE; i++) {
    //     printf("cache[%d]: %u\n", i, bc->cnt[i]);
    // }

    printf("test_host_side_dlmpht_multi_threads passed\n");
    printf("cache hit: %u, total: %u\n", cache_hit, cache_total);
    // break;
}

void dlmpht_local_reconstruct(void* inner_index){
    // printf("dlmpht_local_reconstruct\n");
    
    // initialize Lindex_meta and level0_segment
    load_16B_with_tmp_spram(lmpthash_lindex_meta_spram, inner_index, 1);
    Lindex_metadata* lmd_ = (Lindex_metadata*)lmpthash_lindex_meta_spram;
    uint16_t sidx_level0 = lmd_->level_offsets[lmd_->num_levels - 1];
    load_16B_with_tmp_spram(lmpthash_lindex_meta_spram+16, inner_index, sidx_level0);
    load_16B_with_tmp_spram(lmpthash_lindex_meta_spram+32, inner_index, sidx_level0+1);

    uint16_t l = lmd_->num_levels;
    uint32_t epsilon = lmd_->epsilon+1;
    uint16_t* level_offsets = lmd_->level_offsets;


    for(int offset=level_offsets[l+1]; offset<level_offsets[l+2]; offset++){
        // printf("update %d\n", offset);
        clmpthash_htl_segment* htl_seg = (clmpthash_htl_segment*)lmphash_cache_segment_get_spram;
        memcpy(htl_seg, inner_index+offset*16, 16);
        uint8_t version=(htl_seg->meta >> 54)&0xff;
        // printf("update version: %d ", version);
        version++;
        // printf(" to %d\n", version);
        htl_seg->meta &= 0xc0ffffffffffffff;
        htl_seg->meta |= ((uint64_t)version & 0xff) << 54;
        memcpy(inner_index+offset*16, htl_seg, 16);
    }
    // printf("dlmpht_local_reconstruct done\n");
}

void test_host_side_dlmpht_multi_threads(char* config, int num_threads) {
    // build index
    clmpthash_config cfg;
    clmpthash_lva* lvas;
    clmpthash_physical_addr* pas;
    clmpthash_lva* querys;
    uint64_t num_lva;
    uint64_t num_querys;
    clmpthash_parse_configuration(config, &cfg, &lvas, &pas, &num_lva, &querys, &num_querys);

    void* index = clmpthash_build_index(lvas, pas, num_lva, &cfg);
    if (index == NULL) {
        printf("error building index\n");
        return -1;
    }

    void* inner_index = clmpthash_offload_index(index);
    if (inner_index == NULL) {
        printf("error offloading index\n");
        return -1;
    }

    // nof_dpu_offload_index(inner_index);
    // return;

    printf("query with compacted inner index, for DPU to use\n");
    //
    uint16_t* p16;
    clmpthash_physical_addr pa1;
    bucket_cache* bc = bucket_cache_init(CACHE_SIZE);

    uint64_t total_ls_cnt = 0;
    uint64_t sm_cnt[4] = {0};
    uint64_t total_dma_cnt=0;
    uint64_t bucket_dma_cnt=0;

    uint8_t cache_last_htl_segment=1;
    uint8_t cache_last_htl_segment_hit=0;
    if(cache_last_htl_segment){
        printf("Enabled cache last htl segment\n");
    }


    // 创建多个线程 执行查询
    pthread_t threads[num_threads];
    dlmpht_thread_arg args[num_threads];
    for (int i = 0; i < num_threads; i++) {
        args[i].bc = bc;
        args[i].cache_last_htl_segment = cache_last_htl_segment;
        args[i].inner_index = inner_index;
        args[i].querys = querys;
        args[i].num_querys = num_querys;
        args[i].thread_spram=(spram*)malloc(sizeof(spram)*64);
        memset(args[i].thread_spram, 0xff, sizeof(spram)*64);
        int rc=pthread_create(&threads[i], NULL, dlmpht_query_func, &args[i]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
        sleep(1);
    }

    dlmpht_local_reconstruct(inner_index);

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }

    // printf("cache hit ratio: %lf, total access: %lu\n",(double)(bc->hit_count) / (bc->total_access), bc->total_access);
}

int cmp(const void* a, const void* b) {
    return *(uint64_t*)a - *(uint64_t*)b;
}
typedef struct {
    void* inner_index;
    uint64_t num_query;
    uint64_t* querys;
    int thread_id;
    uint64_t* lats;
} thread_arg_t;

void* query_function(void* arg){
    thread_arg_t* t=(thread_arg_t*)arg;
    for (uint64_t q = 0; q < t->num_query; ++q) {
        // printf("query %lu\n", i);
        // if (q % 10000 == 0) {
        //     printf("\r    %lu / %lu\n", q, num_querys);
        //     printf("\033[1A");
        // }
        
        // 开始时间，初始化
        struct timespec start, end;
        start.tv_sec = 0;
        start.tv_nsec = 0;
        end.tv_sec = 0;
        end.tv_nsec = 0;
        clock_gettime(CLOCK_MONOTONIC, &start);

        clmpthash_lva lva = t->querys[q];
        // printf("query: 0x%lx\n", lva);

        uint64_t l1_addr=L1_SEG_ADDR(lva);
        uint64_t* l1_table=(uint64_t*)(t->inner_index);

        uint64_t l2_table_addr=l1_table[l1_addr];
        uint64_t* l2_table=(uint64_t*)l2_table_addr;
        uint64_t l2_addr=L2_SEG_ADDR(lva);

        uint64_t l3_table_addr=l2_table[l2_addr];
        clmpthash_physical_addr* l3_table=(clmpthash_physical_addr*)l3_table_addr;
        uint64_t l3_addr=L3_SEG_ADDR(lva);

        clmpthash_physical_addr pa1=l3_table[l3_addr];
        clmpthash_lva _lva = 0;
        for (int j = 0; j < 8; j++) { _lva = (_lva << 8) + pa1.data[j]; }
        if (_lva != lva) {
            printf("%lu: wrong result, should be %lu, but got %lu\n", q, lva, _lva);
            pthread_exit(NULL);
        }
        // sleep(1);
        // usleep(1);
        // 结束时间
        clock_gettime(CLOCK_MONOTONIC, &end);
        uint64_t lat=(end.tv_sec-start.tv_sec)*1000000000+(end.tv_nsec-start.tv_nsec);
        t->lats[(t->thread_id)*(t->num_query)+q]=lat;
    }
    pthread_exit(NULL);
}

int test_pt(char* config){
    clmpthash_config cfg;
    clmpthash_lva* lvas;
    clmpthash_physical_addr* pas;
    clmpthash_lva* querys;
    uint64_t num_lva;
    uint64_t num_querys;
    clmpthash_parse_configuration(config, &cfg, &lvas, &pas, &num_lva, &querys, &num_querys);

    void* inner_index = cpt_build_index(lvas, pas, num_lva, &cfg);
    if (inner_index == NULL) {
        printf("error building index\n");
        return -1;
    }
    printf("build index done\n");

    //创建指定数量的线程，同时执行查询
    int num_threads=16;
    // num_querys=10;

    uint64_t* ls=(uint64_t*)malloc(sizeof(uint64_t)*num_querys*num_threads);
    pthread_t threads[num_threads];
    thread_arg_t args[num_threads];

    for (int i = 0; i < num_threads; i++) {
        args[i].inner_index=inner_index;
        args[i].num_query=num_querys;
        args[i].querys=querys;
        args[i].thread_id=i;
        args[i].lats=ls;
        int rc=pthread_create(&threads[i], NULL, query_function, &args[i]);
        if (rc) {
            printf("ERROR; return code from pthread_create() is %d\n", rc);
            exit(-1);
        }
    }

    for (int i = 0; i < num_threads; i++) {
        pthread_join(threads[i], NULL);
    }
    
    qsort(ls, num_querys*num_threads, sizeof(uint64_t), cmp);
    // average latency
    uint64_t total_ls=0;
    for (int i = 0; i < num_querys*num_threads; i++) {
        total_ls+=ls[i];
    }
    printf("average latency: %lf us\n",(double)(total_ls) / (num_querys*num_threads) / 1000);

    // tail latency
    printf("50-th percentile latency: %lf us\n",(double)(ls[(uint64_t)(num_querys*num_threads*0.5)]) / 1000);
    printf("90-th percentile latency: %lf us\n",(double)(ls[(uint64_t)(num_querys*num_threads*0.9)]) / 1000);
    printf("99-th percentile latency: %lf us\n",(double)(ls[(uint64_t)(num_querys*num_threads*0.99)]) / 1000);
    printf("99.9-th percentile latency: %lf us\n",(double)(ls[(uint64_t)(num_querys*num_threads*0.999)]) / 1000);
    printf("99.99-th percentile latency: %lf us\n",(double)(ls[(uint64_t)(num_querys*num_threads*0.9999)]) / 1000);
    free(ls);
}


void build_index_with_scale(clmpthash_lva* lvas, clmpthash_physical_addr* pas, uint64_t num, clmpthash_config* cfg, int scale){
    printf("=========build index with scale %d=========\n", scale);
    // scale the lvas
    clmpthash_lva* new_lvas=NULL;
    clmpthash_physical_addr* new_pas=NULL;
    clmpthash_lva max_lva=lvas[num-1];
    printf("max lva: %lu\n", max_lva);
    uint64_t old_num=num;
    num=num*scale;
    new_lvas=(clmpthash_lva*)malloc(num*sizeof(clmpthash_lva));
    if (new_lvas == NULL) {
        printf("error allocating memory for lvas\n");
        return;
    }
    new_pas=(clmpthash_physical_addr*)malloc(num*sizeof(clmpthash_physical_addr));
    if (new_pas == NULL) {
        printf("error allocating memory for lvas\n");
        return;
    }
    for (int i = 0; i < num; i++) {
        new_lvas[i]=lvas[i%old_num]+(i/old_num)*max_lva;
        new_pas[i]=pas[i%old_num];
    }
    printf("old num: %lu, new num: %lu\n", old_num, num);
    printf("--------build page table--------\n");
    cpt_build_index(new_lvas, new_pas, num, cfg);

    printf("--------build learned table--------\n");
    clt_build_index(new_lvas, new_pas, num, cfg);

    printf("--------build lmpht--------\n");

    void* index = clmpthash_build_index(new_lvas, new_pas, num, cfg);
    if (index == NULL) {
        printf("error building index\n");
        return -1;
    }

    clmpthash_offload_index(index);

    free(new_lvas);
    free(new_pas);
}

void scalability_benchmarks(char* config){
    clmpthash_config cfg;
    clmpthash_lva* lvas;
    clmpthash_physical_addr* pas;
    clmpthash_lva* querys;
    uint64_t num_lva;
    uint64_t num_querys;
    clmpthash_parse_configuration(config, &cfg, &lvas, &pas, &num_lva, &querys, &num_querys);

    build_index_with_scale(lvas, pas, num_lva, &cfg, 1);
    build_index_with_scale(lvas, pas, num_lva, &cfg, 2);
    build_index_with_scale(lvas, pas, num_lva, &cfg, 4);
    build_index_with_scale(lvas, pas, num_lva, &cfg, 8);
    build_index_with_scale(lvas, pas, num_lva, &cfg, 16);
}


int main(int argc, char** argv) {
    // test_host_side_clmpthash(argv[1]);
    // test_host_side_compacted(argv[1]);

    // test_host_side_compacted_v2(argv[1]);
    // test_host_side_compacted_batch_without_align(argv[1]);

    // test_lt(argv[1]);

    // test_pt(argv[1]);

    // 根据参数选择不同的功能
    if (strcmp(argv[1], "scalability") == 0) {
        scalability_benchmarks(argv[2]);
    }else if(strcmp(argv[1], "dlmpht-threads") == 0){
        int num_threads = atoi(argv[2]);
        test_host_side_dlmpht_multi_threads(argv[3], num_threads);
    }else{
        printf("unknown command: %s\n", argv[1]);
        printf("Usage: ./clmpthash_test scalability <config>\n");
        return -1;
    }
    return 0;
}
