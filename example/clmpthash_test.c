#include "../include/clmpthash.h"

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

typedef struct bucket_cache {
    uint32_t table_size;
    bucket_cache_entry* data;

    uint64_t total_access;
    uint64_t hit_count;
} bucket_cache;

bucket_cache* bucket_cache_init(uint32_t table_size) {
    bucket_cache* bc = (bucket_cache*)malloc(sizeof(bucket_cache));
    memset(bc, 0, sizeof(bucket_cache));
    bc->table_size = table_size;
    bc->data = malloc(table_size * sizeof(bucket_cache_entry));
    memset(bc->data, 0xFF,
           table_size * sizeof(bucket_cache_entry));  // cache initialized to 0xFF, in case of key=0
    return bc;
}

uint64_t bucket_cache_get(bucket_cache* bc, uint32_t seg_id, uint32_t compressed_blk) {
    uint64_t key = ((uint64_t)seg_id << 32) + (uint64_t)compressed_blk;
    uint32_t idx = MurmurHash2_64(&key, sizeof(key), 0) % (bc->table_size);
    bc->total_access = bc->total_access + 1;
    if (bc->data[idx].key == key) {
        bc->hit_count = bc->hit_count + 1;
        return bc->data[idx].value;
    } else {
        return UINT64_MAX;
    }
}

uint64_t bucket_cache_get_v2(bucket_cache* bc, uint64_t seg_first_key, uint32_t compressed_blk) {
    uint64_t key = (seg_first_key << 32) + (uint64_t)compressed_blk;
    uint32_t idx = MurmurHash2_64(&key, sizeof(key), 0) % (bc->table_size);
    bc->total_access = bc->total_access + 1;
    if (bc->data[idx].key == key) {
        bc->hit_count = bc->hit_count + 1;
        return bc->data[idx].value;
    } else {
        return UINT64_MAX;
    }
}

void bucket_cache_put(bucket_cache* bc, uint32_t seg_id, uint32_t compressed_blk, uint64_t value) {
    uint64_t key = ((uint64_t)seg_id << 32) + (uint64_t)compressed_blk;
    uint32_t idx = MurmurHash2_64(&key, sizeof(key), 0) % (bc->table_size);
    bc->data[idx].key = key;
    bc->data[idx].value = value;
}

void bucket_cache_put_v2(bucket_cache* bc, uint64_t seg_first_key, uint32_t compressed_blk,
                         uint64_t value) {
    uint64_t key = (seg_first_key << 32) + (uint64_t)compressed_blk;
    uint32_t idx = MurmurHash2_64(&key, sizeof(key), 0) % (bc->table_size);
    bc->data[idx].key = key;
    bc->data[idx].value = value;
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

spram global_spram[32];

#define The_Other_Segment_Spram_ID(id) (id ^ 1)




int test_host_side_compacted_v2(char* config) {
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
    bucket_cache* bc = bucket_cache_init(1024 * 1024 / 16);

    uint64_t total_ls_cnt = 0;
    uint64_t sm_cnt[4] = {0};
    uint8_t cache_last_htl_segment=1;
    uint8_t cache_last_htl_segment_hit=0;
    if(cache_last_htl_segment){
        printf(" ** cache last htl segment\n");
    }
    // num_querys = 10;

    // initialize Lindex_meta and level0_segment
    memcpy(&(global_spram[0].data[0]), inner_index + 16, 16);
    Lindex_metadata* lmd_ = (Lindex_metadata*)(&(global_spram[0].data[0]));
    uint16_t sidx_level0 = lmd_->level_offsets[lmd_->num_levels - 1];
    memcpy((&(global_spram[0].data[1])), inner_index + 32 + sidx_level0 * 16, 16);
    memcpy((&(global_spram[0].data[2])), inner_index + 32 + (sidx_level0 + 1) * 16, 16);

    uint8_t i;
    // num_querys=10;
    for (uint64_t q = 0; q < num_querys; ++q) {
        // printf("query %lu\n", i);
        if (q % 10000 == 0) {
            printf("\r    %lu / %lu\n", q, num_querys);
            printf("\033[1A");
        }
        clmpthash_lva lva = querys[q];

        // query pgm index

        // cached last htl segment
        clmpthash_htl_segment* htl_seg = (clmpthash_htl_segment*)(&(global_spram[11]));
        uint64_t s_first_key;
        cache_last_htl_segment_hit=0;
        // printf("1: lva %lx\n", lva);
        if(cache_last_htl_segment){
            // cached_first_key <= lva
            uint64_t* cached_firsy_key=(uint64_t*)(&(global_spram[11].data[1]));
            // printf("2: cached_firsy_key %lx\n", *cached_firsy_key);
            if(*cached_firsy_key <= lva){
                // cached_next_key > lva
                uint64_t* cached_next_key=(uint64_t*)(&(global_spram[10])) + cached_firsy_key[1];
                // printf("3: cached_next_key %lx\n", *cached_next_key);
                if(*cached_next_key > lva){
                    s_first_key=*cached_firsy_key;
                    cache_last_htl_segment_hit=1;
                }
            }
        }

        if(!cache_last_htl_segment || !cache_last_htl_segment_hit){
            // get Lindex_metadata
            Lindex_metadata* lmd = (Lindex_metadata*)(&(global_spram[0].data[0]));
            uint16_t l = lmd->num_levels;
            uint32_t epsilon = lmd->epsilon;
            uint16_t* level_offsets = lmd->level_offsets;
            // uint32_t num_htl_segment=level_offsets[l+2]-level_offsets[l+1];
            clmpthash_pgm_segment* level0_segments = (clmpthash_pgm_segment*)(&(global_spram[0].data[1]));

            // predict with cached level[l] segment
            int64_t pos =(int64_t)(level0_segments->slope) * (lva - level0_segments->key) / ((uint32_t)1 << 31) +level0_segments->intercept;
            // printf("key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva,level0_segments->key, level0_segments->slope, level0_segments->intercept, pos);
            if (pos < 0) pos = 0;
            if (pos > level0_segments[1].intercept) pos = level0_segments[1].intercept;


            clmpthash_pgm_segment* l1_segs[2];
            l1_segs[0] = (clmpthash_pgm_segment*)(&(global_spram[1]));
            l1_segs[1] = (clmpthash_pgm_segment*)(&(global_spram[2]));  // use two sprams for segments
            uint32_t s_idx1;
            uint32_t s_idx2;
            while (l > 1) {
                // search and compare segment first key at level[l-1]
                s_idx2 = level_offsets[l - 2] + PGM_SUB_EPS(pos, epsilon + 1);  // searching start point at next level
                uint32_t s_idx1 = s_idx2;
                uint8_t active_spram = 1;
                while (true) {
                    active_spram = (active_spram ^ 1);
                    // printf("active_spram: %u\n", active_spram);
                    // load 4 segments at a time
                    memcpy(l1_segs[active_spram], inner_index + 32 + (s_idx2 + 1) * 16, 64);
                    total_ls_cnt++;
                    sm_cnt[0]++;

                    i=0;
                    for (; i < 4; i++) {
                        if (l1_segs[active_spram][i].key > lva) { break; }
                    }
                    s_idx2 += i;
                    if (i < 4) {
                        if (i == 0) {
                            // printf("s_idx1: %u, s_idx2: %u\n", s_idx1, s_idx2);
                            // lva < l1_segs[l1_seg_id][0].key
                            if (s_idx1 == s_idx2) {
                                // start segment happen to be the right segment
                                // load the segment on the other spram
                                memcpy(l1_segs[active_spram ^ 1] + 3, inner_index + 32 + (s_idx1) * 16, 16);total_ls_cnt++;sm_cnt[1]++;
                            }
                            // third segment on the other spram is the right segment in level[l-1]
                            pos = (int64_t)(l1_segs[active_spram ^ 1][3].slope) * (lva - l1_segs[active_spram ^ 1][3].key) / ((uint32_t)1 << 31) +  l1_segs[active_spram ^ 1][3].intercept;
                            // printf("[i==0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n",lva, l1_segs[active_spram ^ 1][3].key,l1_segs[active_spram ^ 1][3].slope,l1_segs[active_spram ^ 1][3].intercept, pos);
                            if (pos < 0) pos = 0;
                            if (pos > l1_segs[active_spram][0].intercept)pos = l1_segs[active_spram][0].intercept;
                        } else {
                            // i=1 or 2 or 3
                            // l1_segs[l1_seg_id][i-1].key <= lva < l1_segs[l1_seg_id][i].key
                            // (i-1)-th segment on this spram is the right segment in level[l-1]
                            pos = (int64_t)(l1_segs[active_spram][i - 1].slope) *      (lva - l1_segs[active_spram][i - 1].key) / ((uint32_t)1 << 31) +  l1_segs[active_spram][i - 1].intercept;
                            // printf("[i>0] key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n",lva, l1_segs[active_spram][i - 1].key,l1_segs[active_spram][i - 1].slope,l1_segs[active_spram][i - 1].intercept, pos);
                            if (pos < 0) pos = 0;

                            if (pos > l1_segs[active_spram][i].intercept)pos = l1_segs[active_spram][i].intercept;
                        }
                        break;
                    }
                }
                l=l-1;
            }

            // reach level-0, search for the right first_key
            
            s_idx1 = PGM_SUB_EPS(pos, epsilon + 1);
            // printf("s_idx1: %u\n", s_idx1);
            uint64_t* fks = (uint64_t*)(&(global_spram[10]));
            memcpy(fks, inner_index + 32 + (level_offsets[lmd->num_levels] + s_idx1 / 2) * 16, 64);
            total_ls_cnt++;
            sm_cnt[2]++;
            s_first_key = fks[0];
            s_idx1 = s_idx1 & (~(uint32_t)1);
            // printf("s_idx1: %u, s_first_key: %lu\n", s_idx1, s_first_key);
            do {
                i=0;
                for (; i < 8; i++) {
                    // printf("i: %u, fks[i]: %lu, lva: %lu\n", i, fks[i], lva);
                    if (fks[i] > lva) { break; }
                    s_first_key = fks[i];
                }
                s_idx1 += i;
                if (i < 8) {
                    break;
                } else {
                    memcpy(fks, inner_index + 32 + (level_offsets[lmd->num_levels] + s_idx1 / 2) * 16, 64);
                    sm_cnt[3]++;
                    total_ls_cnt++;
                }
            } while (1);
            if (s_idx1 > 0) --s_idx1;
            // printf("s_idx1: %u, s_first_key: %lu\n", s_idx1, s_first_key);

            // load htl segment
            
            memcpy(htl_seg, inner_index + 32 + (level_offsets[lmd->num_levels + 1] + s_idx1) * 16, 16);

            // save current first_key to cache
            htl_seg[1].addr = s_first_key;
            htl_seg[1].meta = i;
        }

        



        uint8_t seg_type = htl_seg->meta >> 62;
        // printf("seg_type: %u\n", seg_type);
        if (seg_type == 0) {
            // printf(" accurate segment\n");
            clmpthash_physical_addr* pas = (clmpthash_physical_addr*)(htl_seg->addr + 8);
            uint64_t pos = lva - s_first_key;
            pa1 = pas[pos];
        } else if (seg_type == 4) {
            printf(" error segment\n");
            return -1;
        } else {
            // printf(" approximate segment\n");
            uint64_t hash = MurmurHash2_64(&lva, sizeof(clmpthash_lva), 0x123456789);
            uint64_t blk;
            uint32_t table_size = htl_seg->meta & 0xffffff;
            uint64_t table_addr = htl_seg->addr + 8;
            uint64_t bucket_addr = table_addr + table_size * sizeof(clmpthash_physical_addr);
            uint32_t slope = (htl_seg->meta >> 32) & 0x3fffffff;
            uint64_t width = (htl_seg->meta & 0xff000000) >> 24;
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
            uint64_t blk_value = bucket_cache_get_v2(bc, s_first_key, compressed_blk_id);
            // printf("lva: %lx, blk: %lx, blk_value: %lx, seg_type: %u\n", lva, blk, blk_value,
            // seg_type);

            blk_value = UINT64_MAX;
            if (blk_value != UINT64_MAX) {
                // cache hit
                p = blk_value;
            } else {
                bucket_addr += blk >> 3;
                p = *(uint64_t*)bucket_addr;
                // cache miss
                // bucket_cache_put(bc, sidx, compressed_blk_id, p);
                bucket_cache_put(bc, s_first_key, compressed_blk_id, p);
            }
            p = p >> (blk & 7) & width;
            // printf("p: %lu\n", p);
            p = MurmurHash2_64(&p, sizeof(uint64_t), 0x123456789);
            // printf("hash: %lu, hash_p: %lu\n", hash, p);

            p = (p ^ hash) % table_size;
            // printf("pos: %lu\n", p);
            clmpthash_physical_addr* pas = (clmpthash_physical_addr*)table_addr;
            pa1 = pas[p];
        }

        clmpthash_lva _lva = 0;
        for (int j = 0; j < 8; j++) { _lva = (_lva << 8) + pa1.data[j]; }
        if (_lva != lva) {
            printf("%lu: wrong result, should be 0x%lx, but got 0x%lx\n", q, lva, _lva);
            return -1;
        }
    }
    // break;

    printf("all query passed\n");
    printf("cache hit ratio: %lf, total access: %lu\n",(double)(bc->hit_count) / (bc->total_access), bc->total_access);
    printf("average SM accesses for each query (approximate=(height-1)*epsilon*2/4+epsilon*2/2/4): %f\n", (double)(total_ls_cnt) / num_querys);
    printf("sm_cnt: \n");
    for (i=0; i < 4; i++) { printf("%d: %f\n", i, (double)(sm_cnt[i]) / num_querys); }
    bucket_cache_clean(bc);
    clmpthash_clean_index(index);
    clmpthash_clean_bufs(lvas, pas, querys);
    clmpthash_clean_offloaded_index(inner_index);
}

int test_host_side_compacted(char* config) {
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
    bucket_cache* bc = bucket_cache_init(1024 * 1024 / 16);

    uint64_t total_ls_cnt = 0;

    num_querys = 885;
    for (uint64_t i = 884; i < num_querys; ++i) {
        if (i % 10000 == 0) {
            printf("\r    %lu / %lu\n", i, num_querys);
            printf("\033[1A");
        }

        clmpthash_lva lva = querys[i];
        // query pgm index
        p16 = (uint16_t*)(inner_index + 16);
        uint32_t num_levels = p16[0];
        uint32_t epsilon = p16[1];

        uint16_t* level_offsets = p16 + 2;
        uint32_t num_htl_segment = level_offsets[num_levels + 2] - level_offsets[num_levels + 1];
        // for(int j=0;j<=num_levels+2;j++){
        //     printf("%d ", p16[j+2]);
        // }
        uint16_t l = num_levels;
        clmpthash_pgm_segment* segs = (clmpthash_pgm_segment*)(inner_index + 32);
        uint32_t sidx = level_offsets[num_levels - 1];
        while (l > 0) {
            int64_t pos =
                (int64_t)(segs[sidx].slope) * (lva - segs[sidx].key) / ((uint32_t)1 << 31) +
                segs[sidx].intercept;
            printf("key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva,
                   segs[sidx].key, segs[sidx].slope, segs[sidx].intercept, pos);
            if (pos < 0) pos = 0;
            if (pos > segs[sidx + 1].intercept) pos = segs[sidx + 1].intercept;
            if (l <= 1) {
                sidx = PGM_SUB_EPS(pos, epsilon + 1);
                break;
            }
            uint32_t s_idx2 = level_offsets[l - 2] + PGM_SUB_EPS(pos, epsilon + 1);

            while (segs[s_idx2 + 1].key <= lva) {
                // printf("    s_idx2.key: %lu, lva: %lu\n", segs[s_idx2+1].key, lva);
                total_ls_cnt++;
                ++s_idx2;
            }
            sidx = s_idx2;
            --l;
        }
        // search for bottmom level
        // printf("sidx: %u\n", sidx);
        uint64_t* htl_first_key =
            (uint64_t*)(inner_index + 32 +
                        level_offsets[num_levels] * sizeof(clmpthash_pgm_segment));
        while ((sidx + 1) < num_htl_segment && htl_first_key[sidx + 1] <= lva) {
            ++sidx;
            total_ls_cnt++;
        }
        // printf("sidx: %u, htl_first_key: %lu\n", sidx, htl_first_key[sidx]);
        // printf("lva: %lu, sidx: %d\n", lva, sidx);

        clmpthash_htl_segment* htl_segs =
            (clmpthash_htl_segment*)(inner_index + 32 +         level_offsets[num_levels + 1] * sizeof(clmpthash_pgm_segment));
        clmpthash_htl_segment htl_seg = htl_segs[sidx];
        uint8_t seg_type = htl_seg.meta >> 62;
        // printf("seg_type: %u\n", seg_type);
        if (seg_type == 0) {
            // printf(" accurate segment\n");
            clmpthash_physical_addr* pas = (clmpthash_physical_addr*)(htl_seg.addr + 8);
            uint64_t pos = lva - htl_first_key[sidx];
            pa1 = pas[pos];
        } else if (seg_type == 4) {
            printf(" error segment\n");
            return -1;
        } else {
            // printf(" approximate segment\n");
            uint64_t hash = MurmurHash2_64(&lva, sizeof(clmpthash_lva), 0x123456789);
            uint64_t blk;
            uint32_t table_size = htl_seg.meta & 0xffffff;
            uint64_t table_addr = htl_seg.addr + 8;
            uint64_t bucket_addr = table_addr + table_size * sizeof(clmpthash_physical_addr);
            uint32_t slope = (htl_seg.meta >> 32) & 0x3fffffff;
            uint64_t width = (htl_seg.meta & 0xff000000) >> 24;
            if (seg_type == 1) {
                // linear bucketing
                blk = (lva - htl_first_key[sidx]) / slope;

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
            uint64_t blk_value = bucket_cache_get_v2(bc, htl_first_key[sidx], compressed_blk_id);
            // printf("lva: %lx, blk: %lx, blk_value: %lx, seg_type: %u\n", lva, blk, blk_value,
            // seg_type);

            blk_value = UINT64_MAX;
            if (blk_value != UINT64_MAX) {
                // cache hit
                p = blk_value;
            } else {
                bucket_addr += blk >> 3;
                p = *(uint64_t*)bucket_addr;
                // cache miss
                // bucket_cache_put(bc, sidx, compressed_blk_id, p);
                bucket_cache_put(bc, htl_first_key[sidx], compressed_blk_id, p);
            }
            p = p >> (blk & 7) & width;
            // printf("p: %lu\n", p);
            p = MurmurHash2_64(&p, sizeof(uint64_t), 0x123456789);
            // printf("hash: %lu, hash_p: %lu\n", hash, p);

            p = (p ^ hash) % table_size;
            // printf("pos: %lu\n", p);
            clmpthash_physical_addr* pas = (clmpthash_physical_addr*)table_addr;
            pa1 = pas[p];
        }

        clmpthash_lva _lva = 0;
        for (int j = 0; j < 8; j++) { _lva = (_lva << 8) + pa1.data[j]; }
        if (_lva != lva) {
            printf("%lu: wrong result, should be 0x%lx, but got 0x%lx\n", i, lva, _lva);
            return -1;
        }
        // break;
    }
    printf("all query passed\n");
    printf("cache hit ratio: %lf, total access: %lu\n",
           (double)(bc->hit_count) / (bc->total_access), bc->total_access);
    printf("average local search steps for each query: %f\n", (double)(total_ls_cnt) / num_querys);
    bucket_cache_clean(bc);
    clmpthash_clean_index(index);
    clmpthash_clean_bufs(lvas, pas, querys);
    clmpthash_clean_offloaded_index(inner_index);
}

// int nof_dpu_offload_index(void* index){
//     uint32_t num_levels=*((uint16_t*)(index+16));
//     uint16_t* level_offsets=(uint16_t*)(index+20);

//     // visit and offload table data
//     clmpthash_htl_segment*
//     htl_segs=(clmpthash_htl_segment*)(index+32+level_offsets[num_levels+1]*sizeof(clmpthash_pgm_segment));
//     uint32_t num_htl_segment=level_offsets[num_levels+2]-level_offsets[num_levels+1];
//     for(uint32_t i=0;i<num_htl_segment;i++){
//         uint64_t* table=(uint64_t*)(htl_segs[i].addr);
//         uint8_t seg_type=(htl_segs[i].meta)>>62;
//         uint64_t table_size=table[0];
//         // NOTE: current kernel driver noly support allocate 4MB dma buffer
//         // TODO: alloc huge dma buffer by IOMMU（Input-Output Memory Management Unit） or
//         CMA（Contiguous Memory Allocator), which may need to modify the kernel configuration

//         // PS: Only Accurate Segment will generage such a huge dma buffer,
//         if(seg_type==0){
//             // accurate segment
//             uint32_t offload_size = (table_size-8) > 4*1024*1024 ? 4*1024*1024 : (table_size-8);
//             // offload_to_kernel(table+8, offload_size);
//         }else{
//             uint32_t offload_size=(table_size-8);
//             if(offload_size>4*1024*1024){
//                 printf("error: offload size too large for a approximate segment\n");
//                 return -1;
//             }
//             // offload_to_kernel(table+8, offload_size);
//         }
//     }

//     // offload Lindex to DPU
//     offload_to_DPU(index, level_offsets[num_levels+2]*16);

// }

int main(int argc, char** argv) {
    // test_host_side_clmpthash(argv[1]);
    // test_host_side_compacted(argv[1]);

    test_host_side_compacted_v2(argv[1]);

    return 0;
}
