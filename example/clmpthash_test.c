#include "../include/clmpthash.h"

#define PGM_SUB_EPS(x, epsilon) ((x) <= (epsilon) ? 0 : ((x) - (epsilon)))

int test_host_side_clmpthash(char* config){
    clmpthash_config cfg;
    LVA* lvas;
    PhysicalAddr* pas;
    LVA* querys;
    uint64_t num_lva, num_querys;
    parse_configuration(config, &cfg, &lvas, &pas, &num_lva, &querys, &num_querys);

    void* index = build_index(lvas, pas, num_lva, &cfg);
    if(index==NULL){
        printf("error building index\n");
        return -1;
    }

    PhysicalAddr pa;
    for(int i = 0; i < num_querys; ++i) {
        if(i%10000==0){
            printf("\r    %d / %d\n", i, num_querys);
            printf("\033[1A");
        }

        int ret = get_pa(querys[i], index, &pa);
        
        // check if the result is correct: first 8 bytes should be equal to the original lva
        if(ret==0){
            LVA _lva=0;
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
    clean_index(index);
    clean_bufs(lvas, pas, querys);
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


int test_host_side_compacted(char* config){
    clmpthash_config cfg;
    LVA* lvas;
    PhysicalAddr* pas;
    LVA* querys;
    uint64_t num_lva, num_querys;
    parse_configuration(config, &cfg, &lvas, &pas, &num_lva, &querys, &num_querys);

    void* index = build_index(lvas, pas, num_lva, &cfg);
    if(index==NULL){
        printf("error building index\n");
        return -1;
    }

    void* inner_index=offload_index(index);
    if(inner_index==NULL){
        printf("error offloading index\n");
        return -1;
    }

    printf("query with compacted inner index, for DPU to use\n");
    // 
    uint16_t* p16;
    PhysicalAddr pa1;

    for(int i = 16734; i < num_querys; ++i) {
        if(i%10000==0){
            printf("\r    %d / %d\n", i, num_querys);
            printf("\033[1A");
        }

        LVA lva = querys[i];
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
        pgm_segment* segs=(pgm_segment*)(inner_index+32);
        uint32_t sidx=level_offsets[num_levels-1];
        while(l>0){
            int64_t pos=(int64_t)(segs[sidx].slope)*(lva-segs[sidx].key) / ((uint32_t)1<<31) + segs[sidx].intercept;

            if(pos<0) pos=0;
            if(pos>segs[sidx+1].intercept)  pos=segs[sidx+1].intercept;
            if(l<=1){
                sidx=PGM_SUB_EPS(pos, epsilon+1);
                break;
            }
            uint32_t s_idx2=level_offsets[l-2]+PGM_SUB_EPS(pos, epsilon+1);
            while(segs[s_idx2].key > lva){
                ++s_idx2;
            }
            sidx=s_idx2;
            --l;
        }
        // search for bottmom level
        uint64_t* htl_first_key=(uint64_t*)(inner_index+32+level_offsets[num_levels]*sizeof(pgm_segment));
        while((sidx+1)<num_htl_segment && htl_first_key[sidx+1]<=lva){
            ++sidx;
        }
        // printf("lva: %lu, sidx: %d\n", lva, sidx);
        
        htl_segment* htl_segs=(htl_segment*)(inner_index+32+level_offsets[num_levels+1]*sizeof(pgm_segment));
        htl_segment htl_seg=htl_segs[sidx];
        uint8_t seg_type=htl_seg.meta >> 62;
        // printf("seg_type: %u\n", seg_type);
        if(seg_type==0){
            // printf(" accurate segment\n");
            PhysicalAddr* pas=(PhysicalAddr*)(htl_seg.addr+8);
            uint64_t pos=lva-htl_first_key[sidx];
            pa1=pas[pos];
        }else if(seg_type==4){
            printf(" error segment\n");
            return -1;
        }else{
            // printf(" approximate segment\n");
            uint64_t hash=MurmurHash2_64(&lva, sizeof(LVA), 0x123456789);
            uint64_t blk;
            uint32_t table_size=htl_seg.meta & 0xffffff;
            uint64_t table_addr=htl_seg.addr+8;
            uint64_t bucket_addr=table_addr+table_size*sizeof(PhysicalAddr);
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
            bucket_addr+=blk>>3;
            
            uint64_t p=*(uint64_t*)bucket_addr;
            p=p>>(blk&7)&width;
            // printf("p: %lu\n", p);
            p=MurmurHash2_64(&p, sizeof(uint64_t), 0x123456789);
            // printf("hash: %lu, hash_p: %lu\n", hash, p);
            
            p=(p^hash)%table_size;
            // printf("pos: %lu\n", p);
            PhysicalAddr* pas=(PhysicalAddr*)table_addr;
            pa1=pas[p];
        }
        
        LVA _lva=0;
        for(int j=0;j<8;j++){
            _lva=(_lva<<8)+pa1.data[j];
        }
        if(_lva!=lva){
            printf("%d: wrong result, should be 0x%lx, but got 0x%lx\n", i, lva, _lva);
            return -1;
        }
        // break;
        
    }
    printf("all query passed\n");
    clean_index(index);
    clean_bufs(lvas, pas, querys);
    clean_offloaded_index(inner_index);

}


int main(int argc, char** argv) {
    test_host_side_clmpthash(argv[1]);
    // test_host_side_compacted(argv[1]);
    return 0;

}