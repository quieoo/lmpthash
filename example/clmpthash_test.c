#include "../include/clmpthash.h"

#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))

int nof_dpu_offload_index(void* index);

int test_host_side_clmpthash(char* config){
    clmpthash_config cfg;
    clmpthash_lva* lvas;
    clmpthash_physical_addr* pas;
    clmpthash_lva* querys;
    uint64_t num_lva, num_querys;
    clmpthash_parse_configuration(config, &cfg, &lvas, &pas, &num_lva, &querys, &num_querys);

    void* index = clmpthash_build_index(lvas, pas, num_lva, &cfg);
    if(index==NULL){
        printf("error building index\n");
        return -1;
    }

    clmpthash_physical_addr pa;
    for(uint64_t i = 0; i < num_querys; ++i) {
        if(i%10000==0){
            printf("\r    %lu / %lu\n", i, num_querys);
            printf("\033[1A");
        }

        int ret = clmpthash_get_pa(querys[i], index, &pa);
        
        // check if the result is correct: first 8 bytes should be equal to the original lva
        if(ret==0){
            clmpthash_lva _lva=0;
            for(int j=0;j<8;j++){
                _lva=(_lva<<8)+pa.data[j];
            }
            if(_lva!=querys[i]){
                printf("wrong result: should be 0x%lx, but got 0x%lx\n", querys[i], _lva);
                return -1;
            }
        }else{
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

typedef struct bucket_cache_entry{
    uint64_t key;
    uint64_t value;
}bucket_cache_entry;

typedef struct bucket_cache{
    uint32_t table_size;
    bucket_cache_entry* data;

    uint64_t total_access;
    uint64_t hit_count;
}bucket_cache;

bucket_cache* bucket_cache_init(uint32_t table_size){
    bucket_cache* bc=(bucket_cache*)malloc(sizeof(bucket_cache));
    memset(bc, 0, sizeof(bucket_cache));
    bc->table_size=table_size;
    bc->data=malloc(table_size*sizeof(bucket_cache_entry));
    memset(bc->data, 0xFF, table_size*sizeof(bucket_cache_entry));  // cache initialized to 0xFF, in case of key=0
    return bc;
}

uint64_t bucket_cache_get(bucket_cache* bc, uint32_t seg_id, uint32_t compressed_blk){
    uint64_t key=((uint64_t)seg_id<<32)+(uint64_t)compressed_blk;
    uint32_t idx=MurmurHash2_64(&key, sizeof(key), 0)%(bc->table_size);
    bc->total_access=bc->total_access+1;
    if(bc->data[idx].key==key){
        bc->hit_count=bc->hit_count+1;
        return bc->data[idx].value;
    }else{
        return UINT64_MAX;
    }
}

void bucket_cache_put(bucket_cache* bc, uint32_t seg_id, uint32_t compressed_blk, uint64_t value){
    uint64_t key=((uint64_t)seg_id<<32)+(uint64_t)compressed_blk;
    uint32_t idx=MurmurHash2_64(&key, sizeof(key), 0)%(bc->table_size);
    bc->data[idx].key=key;
    bc->data[idx].value=value;
}

void bucket_cache_clean(bucket_cache* bc){
    free(bc->data);
    free(bc);
}


int test_host_side_compacted(char* config){
    clmpthash_config cfg;
    clmpthash_lva* lvas;
    clmpthash_physical_addr* pas;
    clmpthash_lva* querys;
    uint64_t num_lva;
    uint64_t num_querys;
    clmpthash_parse_configuration(config, &cfg, &lvas, &pas, &num_lva, &querys, &num_querys);

    void* index = clmpthash_build_index(lvas, pas, num_lva, &cfg);
    if(index==NULL){
        printf("error building index\n");
        return -1;
    }

    void* inner_index=clmpthash_offload_index(index);
    if(inner_index==NULL){
        printf("error offloading index\n");
        return -1;
    }

    // nof_dpu_offload_index(inner_index);
    // return;

    printf("query with compacted inner index, for DPU to use\n");
    // 
    uint16_t* p16;
    clmpthash_physical_addr pa1;
    bucket_cache* bc=bucket_cache_init(1024*1024/16);

    uint64_t total_ls_cnt=0;

    // num_querys=1;
    for(uint64_t i = 0; i < num_querys; ++i) {
        if(i%10000==0){
            printf("\r    %lu / %lu\n", i, num_querys);
            printf("\033[1A");
        }

        clmpthash_lva lva = querys[i];
        // query pgm index
        p16=(uint16_t*)(inner_index+16);
        uint32_t num_levels=p16[0];
        uint32_t epsilon=p16[1];

        uint16_t* level_offsets=p16+2;
        uint32_t num_htl_segment=level_offsets[num_levels+2]-level_offsets[num_levels+1];
        // for(int j=0;j<=num_levels+2;j++){
        //     printf("%d ", p16[j+2]);
        // }
        uint16_t l=num_levels;
        clmpthash_pgm_segment* segs=(clmpthash_pgm_segment*)(inner_index+32);
        uint32_t sidx=level_offsets[num_levels-1];
        while(l>0){
            int64_t pos=(int64_t)(segs[sidx].slope)*(lva-segs[sidx].key) / ((uint32_t)1<<31) + segs[sidx].intercept;
            // printf("key: %lu, first_key: %lu, slope: %u, intercept: %u, pos: %lu\n", lva, segs[sidx].key, segs[sidx].slope, segs[sidx].intercept, pos);
            if(pos<0) pos=0;
            if(pos>segs[sidx+1].intercept)  pos=segs[sidx+1].intercept;
            if(l<=1){
                sidx=PGM_SUB_EPS(pos, epsilon+1);
                break;
            }
            uint32_t s_idx2=level_offsets[l-2]+PGM_SUB_EPS(pos, epsilon+1);
            
            while(segs[s_idx2+1].key <= lva){
                // printf("    s_idx2.key: %lu, lva: %lu\n", segs[s_idx2+1].key, lva);
                total_ls_cnt++;
                ++s_idx2;
            }
            sidx=s_idx2;
            --l;
        }
        // search for bottmom level
        uint64_t* htl_first_key=(uint64_t*)(inner_index+32+level_offsets[num_levels]*sizeof(clmpthash_pgm_segment));
        while((sidx+1)<num_htl_segment && htl_first_key[sidx+1]<=lva){
            ++sidx;
            total_ls_cnt++;
        }
        // printf("lva: %lu, sidx: %d\n", lva, sidx);
        
        clmpthash_htl_segment* htl_segs=(clmpthash_htl_segment*)(inner_index+32+level_offsets[num_levels+1]*sizeof(clmpthash_pgm_segment));
        clmpthash_htl_segment htl_seg=htl_segs[sidx];
        uint8_t seg_type=htl_seg.meta >> 62;
        // printf("seg_type: %u\n", seg_type);
        if(seg_type==0){
            // printf(" accurate segment\n");
            clmpthash_physical_addr* pas=(clmpthash_physical_addr*)(htl_seg.addr+8);
            uint64_t pos=lva-htl_first_key[sidx];
            pa1=pas[pos];
        }else if(seg_type==4){
            printf(" error segment\n");
            return -1;
        }else{
            // printf(" approximate segment\n");
            uint64_t hash=MurmurHash2_64(&lva, sizeof(clmpthash_lva), 0x123456789);
            uint64_t blk;
            uint32_t table_size=htl_seg.meta & 0xffffff;
            uint64_t table_addr=htl_seg.addr+8;
            uint64_t bucket_addr=table_addr+table_size*sizeof(clmpthash_physical_addr);
            uint32_t slope=(htl_seg.meta >>32) & 0x3fffffff;
            uint64_t width=(htl_seg.meta & 0xff000000)>>24;
            if(seg_type==1){
                // linear bucketing
                blk=(lva-htl_first_key[sidx])/slope;
                
            }else{
                // hash bucketing
                // printf("hash:%lu, slope: %lu\n", hash, slope);
                blk=hash%slope;
            }
            // printf("blk: %lu\n", blk);
            blk=blk*width;
            width=-(width==64) | ((((uint64_t)1)<<width)-1);
            /*
            bucket_addr+=blk>>3;
            uint64_t p=*(uint64_t*)bucket_addr;
            */
            uint64_t p;
            uint32_t compressed_blk_id=blk>>3;
            uint64_t blk_value=bucket_cache_get(bc, sidx, compressed_blk_id);
            blk_value=UINT64_MAX;
            if(blk_value!=UINT64_MAX){
                // cache hit
                p=blk_value;
            }else{
                bucket_addr+=blk>>3;
                p=*(uint64_t*)bucket_addr;
                // cache miss
                bucket_cache_put(bc, sidx, compressed_blk_id, p);
            }
            p=p>>(blk&7)&width;
            // printf("p: %lu\n", p);
            p=MurmurHash2_64(&p, sizeof(uint64_t), 0x123456789);
            // printf("hash: %lu, hash_p: %lu\n", hash, p);
            
            p=(p^hash)%table_size;
            // printf("pos: %lu\n", p);
            clmpthash_physical_addr* pas=(clmpthash_physical_addr*)table_addr;
            pa1=pas[p];
        }
        
        clmpthash_lva _lva=0;
        for(int j=0;j<8;j++){
            _lva=(_lva<<8)+pa1.data[j];
        }
        if(_lva!=lva){
            printf("%lu: wrong result, should be 0x%lx, but got 0x%lx\n", i, lva, _lva);
            return -1;
        }
        // break;
        
    }
    printf("all query passed\n");
    printf("cache hit ratio: %lf, total access: %lu\n", (double)(bc->hit_count)/(bc->total_access), bc->total_access);
    printf("average local search steps for each query: %f\n",(double)(total_ls_cnt)/num_querys);
    bucket_cache_clean(bc);
    clmpthash_clean_index(index);
    clmpthash_clean_bufs(lvas, pas, querys);
    clmpthash_clean_offloaded_index(inner_index);

}


// int nof_dpu_offload_index(void* index){
//     uint32_t num_levels=*((uint16_t*)(index+16));
//     uint16_t* level_offsets=(uint16_t*)(index+20);

//     // visit and offload table data
//     clmpthash_htl_segment* htl_segs=(clmpthash_htl_segment*)(index+32+level_offsets[num_levels+1]*sizeof(clmpthash_pgm_segment));
//     uint32_t num_htl_segment=level_offsets[num_levels+2]-level_offsets[num_levels+1];
//     for(uint32_t i=0;i<num_htl_segment;i++){
//         uint64_t* table=(uint64_t*)(htl_segs[i].addr);
//         uint8_t seg_type=(htl_segs[i].meta)>>62;
//         uint64_t table_size=table[0];
//         // NOTE: current kernel driver noly support allocate 4MB dma buffer
//         // TODO: alloc huge dma buffer by IOMMU（Input-Output Memory Management Unit） or CMA（Contiguous Memory Allocator), which may need to modify the kernel configuration
        
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
    test_host_side_compacted(argv[1]);
    return 0;

}
