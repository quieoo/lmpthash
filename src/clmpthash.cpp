#include "lmpthash.hpp"
#include "clmpthash.h"
#include "include/utils/logger.hpp"

void cparse_msr_cambridge(clmpthash_lva** lvas, clmpthash_physical_addr** pas, uint64_t* num_lva, clmpthash_lva** querys, uint64_t * num_querys, std::string trace_path){
    std::vector<clmpthash_lva> uniq_lpn;
    std::vector<clmpthash_lva> lpns;
    parse_MSR_Cambridge(uniq_lpn, lpns, trace_path);
    *num_lva=uniq_lpn.size();
    *num_querys=lpns.size();

    // *lvas=(LVA*)malloc(sizeof(LVA)*uniq_lpn.size());
    *lvas=new clmpthash_lva[uniq_lpn.size()];
    // *querys=(LVA*)malloc(sizeof(LVA)*lpns.size());
    *querys=new clmpthash_lva[lpns.size()];
    for(uint64_t i=0;i<uniq_lpn.size();i++){
        (*lvas)[i]=uniq_lpn[i];
    }
    for(uint64_t i=0;i<lpns.size();i++){
        (*querys)[i]=lpns[i];
    }
}

void clmpthash_parse_configuration(char* config_path, clmpthash_config* cfg, clmpthash_lva** lvas, clmpthash_physical_addr** pas, uint64_t* num_lva, clmpthash_lva** querys, uint64_t* num_querys){
    lmpthash_config config;
    config.load_config(std::string(config_path));

    // get trace file
    cparse_msr_cambridge(lvas, pas, num_lva, querys, num_querys, config.trace_path);
    printf("    lva_num: %lu, query_num: %lu\n", *num_lva, *num_querys);
    // *pas=(clmpthash_physical_addr*)malloc(sizeof(clmpthash_physical_addr)*(*num_lva));
    *pas=new clmpthash_physical_addr[*num_lva];
    for(uint64_t k=0;k<*num_lva;k++){
        for(int i=0; i<8; i++){
            (*pas)[k].data[i]=((*lvas)[k])>>((7-i)*8);
        }
    }

    // get config
    cfg->alpha=config.alpha;
    cfg->beta=config.beta;
    cfg->gamma=config.gamma;
    cfg->P=config.P;
    cfg->hashed_num_bucket_c=config.hashed_num_bucket_c;
    cfg->table_size_alpha=config.table_size_alpha;
    cfg->max_bucket_size=config.max_bucket_size;
    cfg->pilot_search_threshold=config.pilot_search_threshold;
    cfg->dynamic_alpha=config.dynamic_alpha;
    cfg->alpha_limits=config.alpha_limits;
    cfg->left_epsilon=config.left_epsilon;
    cfg->right_epsilon=config.right_epsilon;
}

void clmpthash_clean_bufs(clmpthash_lva* lvas, clmpthash_physical_addr* pas, clmpthash_lva* querys){
    printf("    clean bufs\n");
    delete[] lvas;
    delete[] pas;
    delete[] querys;
}

void* clmpthash_build_index(clmpthash_lva* lvas, clmpthash_physical_addr* pas, uint64_t num, clmpthash_config* cfg){
    lmpthash_config config;
    config.alpha=cfg->alpha;
    config.beta=cfg->beta;
    config.gamma=cfg->gamma;
    config.P=cfg->P;
    config.hashed_num_bucket_c=cfg->hashed_num_bucket_c;
    config.table_size_alpha=cfg->table_size_alpha;
    config.max_bucket_size=cfg->max_bucket_size;
    config.pilot_search_threshold=cfg->pilot_search_threshold;
    config.dynamic_alpha=cfg->dynamic_alpha;
    config.alpha_limits=cfg->alpha_limits;
    config.left_epsilon=cfg->left_epsilon;
    config.right_epsilon=cfg->right_epsilon;

    LMPTHashBuilder<clmpthash_lva, clmpthash_physical_addr>* builder = new LMPTHashBuilder<clmpthash_lva, clmpthash_physical_addr>(config);
    std::vector<clmpthash_lva> keys(lvas, lvas+num);
    std::vector<clmpthash_physical_addr> values(pas, pas+num);

    builder->Segmenting(keys);
    builder->Learning();
    builder->Multi_Bucketing();
    if(builder->Tabling(keys, values)){
        printf("Build index table failed\n");
        return NULL;
    }

    printf("Build index done\n");

    return static_cast<void*>(builder);
}

int clmpthash_get_pa(clmpthash_lva lva, void* index, clmpthash_physical_addr* pa){
    LMPTHashBuilder<clmpthash_lva, clmpthash_physical_addr>* builder=static_cast<LMPTHashBuilder<clmpthash_lva, clmpthash_physical_addr>*>(index);
    
    pthash::simple_logger slogger;
    // slogger.allowed_func_ids.push_back(4);

    slogger.func_log(4, "lva: %lu\n", lva);
    auto range=builder->pgm_index.search(lva);
    auto lo=range.lo;
    auto hi=range.hi;
    slogger.func_log(4, "lo: %lu, hi: %lu\n", lo, hi);
    if(hi >= builder->lmpt_segments.size()) hi=builder->lmpt_segments.size()-1;
    int ans=-1;
    while(lo <= hi){
        auto mid=(lo+hi)/2;
        slogger.func_log(4, "mid: %lu, first_key: %lu\n", mid, builder->lmpt_segments[mid].first_key);
        if(builder->lmpt_segments[mid].first_key <= lva){
            ans=mid;
            lo=mid+1;
        }else{
            hi=mid-1;
        }
    }
    if(ans==-1){
        printf("error: segment not found\n");
        return -1;
    }
    slogger.func_log(4, "ans: %lu\n", ans);
    int64_t pos=-1;
    if(builder->lmpt_segments[ans].seg_type==0){
        pos=lva-builder->lmpt_segments[ans].first_key;
    }else if(builder->lmpt_segments[ans].seg_type==1 || builder->lmpt_segments[ans].seg_type==2){
        pos=builder->pthash_map[ans](lva);
    }

    clmpthash_physical_addr* subtable=(clmpthash_physical_addr*)(builder->lmpt_segments[ans].next_addr);
    *pa=subtable[pos];

    // output pa
    for(int i=0;i<20;i++){
        slogger.func_log(4, "%02x ", (*pa).data[i]);
    }
    slogger.func_log(4, "\n");

    return 0;
}

int clmpthash_clean_index(void* index){
    LMPTHashBuilder<clmpthash_lva, clmpthash_physical_addr>* builder=static_cast<LMPTHashBuilder<clmpthash_lva, clmpthash_physical_addr>*>(index);
    builder->Cleaning();
    // printf("finish cleanning buffer\n");
    delete builder;
    return 0;
}


void* clmpthash_offload_index(void* index){

    uint8_t* inner_index=new uint8_t[1024*1024];
    memset(inner_index, 0, 1024*1024);

    LMPTHashBuilder<clmpthash_lva, clmpthash_physical_addr>* builder=static_cast<LMPTHashBuilder<clmpthash_lva, clmpthash_physical_addr>*>(index);
    if(builder->Compacting(inner_index)){
        printf("error: compacting failed\n");
        return NULL;
    }
    uint64_t* ptr=(uint64_t*)inner_index;
    return static_cast<void*>(inner_index);
}


int clmpthash_clean_offloaded_index(void* inner_index){
    uint16_t num_level=*((uint16_t*)(inner_index+16));
    uint16_t* level_offsets=(uint16_t*)(inner_index+20);
    uint16_t segments_start=level_offsets[num_level+1];
    uint16_t segments_ends=level_offsets[num_level+2];
    // printf("    clean offloaded index, segment_offset: %d-%d\n", segments_start, segments_ends);

    
    clmpthash_htl_segment* segments=(clmpthash_htl_segment*)(inner_index+32);

    for(int i=segments_start;i<segments_ends;i++){
        uint8_t* ptr=(uint8_t*)(segments[i].addr);
        if(ptr!=nullptr){
            delete[] ptr;
        }
    }

    delete[] inner_index;
    return 0;
}

pthash::simple_logger gslogger;

void* clt_build_index(clmpthash_lva* lvas, clmpthash_physical_addr* pas, uint64_t num, clmpthash_config* cfg){

    // gslogger.allowed_func_ids.push_back(0);

    lmpthash_config config;
    config.alpha=cfg->alpha;
    config.beta=cfg->beta;
    config.gamma=cfg->gamma;
    config.P=cfg->P;
    config.hashed_num_bucket_c=cfg->hashed_num_bucket_c;
    config.table_size_alpha=cfg->table_size_alpha;
    config.max_bucket_size=cfg->max_bucket_size;
    config.pilot_search_threshold=cfg->pilot_search_threshold;
    config.dynamic_alpha=cfg->dynamic_alpha;
    config.alpha_limits=cfg->alpha_limits;
    config.left_epsilon=cfg->left_epsilon;
    config.right_epsilon=cfg->right_epsilon;

    std::vector<clmpthash_lva> keys(lvas, lvas+num);
    std::vector<clmpthash_physical_addr> values(pas, pas+num);
    pgm::PGMIndex<uint64_t, 64,4,uint32_t> pgm_index;
    // sort
    assert(std::is_sorted(keys.begin(), keys.end()));
    
    //binary serach for a minimum epsilon that fits the inner nodes in limited memory
    gslogger.func_log(0, "    binary search for a minimum epsilon that fits the inner nodes in limited memory\n");
    int l=config.left_epsilon, r=config.right_epsilon;
    int ep=-1;

    while(l<=r){
        int mid=(l+r)/2;
        pgm::PGMIndex<uint64_t, 64,4,uint32_t> pgm(keys, mid, mid);
        uint32_t seg_num=pgm.segments_count();
        double seg_size=(double)seg_num*16/1024.0/1024.0;
        gslogger.func_log(0, "    mid: %d, seg_num: %d, size: %f MB\n", mid, seg_num, seg_size);
        // if(seg_size>1.0){
        if(seg_size>0.3){    
            l=mid+1;
        }else{
            r=mid-1;
            ep=mid;
        }
    }

    printf("    min epsilon: %d, \n", ep);
    if(ep==-1){
        printf("error: cannot find a minimum epsilon that fits the inner nodes in limited memory\n");
        return NULL;
    }
    if(ep>13){
        printf("warning: epsilon is too large that cannot be dma to spram in one time\n");
    }

    pgm::PGMIndex<uint64_t, 64,4,uint32_t> pgm(keys, ep, config.right_epsilon);
    int min_height=pgm.height();
    printf("    min height: %d\n", min_height);
    // binary search for a minimum epsilon that generates the minimum height
    gslogger.func_log(0, "    binary search for a minimum epsilon that generates the minimum height\n");
    l=config.left_epsilon, r=config.right_epsilon;
    int epr=-1;
    while(l<=r){
        int mid=(l+r)/2;
        pgm::PGMIndex<uint64_t, 64, 4, uint32_t> pgm(keys, ep, mid);
        int height=pgm.height();
        gslogger.func_log(0, "    mid: %d, height: %d\n", mid, height);
        if(height==min_height){
            epr=mid;
            r=mid-1;
            pgm_index=pgm;
        }else{
            l=mid+1;
        }
    }
    printf("    min_epsilon_recursive: %d\n", epr);

    uint8_t* inner_index=new uint8_t[2*1024*1024];
    memset(inner_index, 0, 2*1024*1024);


    uint8_t* ptr=inner_index;
    // first_key and last_key takes the first 16 bytes
    uint64_t* p64=(uint64_t*)ptr;
    p64[0]=pgm_index.first_key;
    ptr+=16;

    // num level of pgm_index take the next 2 bytes
    uint16_t* p16=(uint16_t*)ptr;
    p16[0]=pgm_index.height();
    if(p16[0]>5){
        printf("Error while Compacting: pgm height > 3\n");
        return NULL;
    }
    ptr+=2;

    // variable_epsilon_value takes the next 2 bytes
    uint16_t ep_value=0;
    ep_value=(pgm_index.variable_er)&0xFF;
    ep_value=ep_value<<8;
    ep_value |= (pgm_index.variable_epsilon_value)&0xFF;

    p16=(uint16_t*)ptr;
    p16[0]=ep_value;
    ptr+=2;
    // array of level_offsets takes the next 12 bytes, each one takes 2 bytes
    p16=(uint16_t*)ptr;
    std::vector<size_t> lofs;
    pgm_index.get_level_offsets(lofs);
    
    int num_level_offsets=lofs.size();
    for(int i=0;i<num_level_offsets;i++){
        p16[i]=lofs[i]+2;
    }
    
    printf("    level_offsets: ");
    for(int i=0; i<num_level_offsets;i++){
        printf("%d ", p16[i]);
    }
    printf("\n");
    ptr+=12;

    // each inner model takes 16 bytes
    std::vector<clmpthash_pgm_segment> segs;
    for(int i=0;i<pgm_index.segments.size();i++){
        clmpthash_pgm_segment ps;
        ps.key=pgm_index.segments[i].key;
        ps.slope=(uint32_t)(pgm_index.segments[i].slope);
        ps.intercept=pgm_index.segments[i].intercept;
        // printf("pgm segment %d: key: %ld, slope: %d, intercept: %d\n", i, ps.key, ps.slope, ps.intercept);
        gslogger.func_log(0, "pgm segment %d: key: %ld, slope: %d, intercept: %d\n", i, ps.key, ps.slope, ps.intercept);
        segs.push_back(ps);
    }
    memcpy(ptr, segs.data(), segs.size()*sizeof(clmpthash_pgm_segment));

    // split lvas and pas to buffers with 2097136 bytes
    uint64_t lva2pa_size=sizeof(clmpthash_lva)+sizeof(clmpthash_physical_addr);
    uint64_t pas_in_buf=2097136/lva2pa_size;
    uint64_t num_bufs = (num + pas_in_buf - 1) / pas_in_buf;
    printf("    num_bufs: %ld\n", num_bufs);

    uint64_t buf_id=0;
    uint64_t* sub_table_addr=(uint64_t*)(inner_index+1024*1024);
    sub_table_addr[buf_id++]=num_bufs;
    uint64_t offset=0;
    while(offset<num){
        uint64_t last=offset+pas_in_buf;
        if(last>num){
            last=num;
        }
        // printf("    offset: %ld, last: %ld\n", offset, last);
        uint8_t* sub_table=new uint8_t[(last-offset)*lva2pa_size];

        for(uint64_t i=offset;i<last;i++){
            memcpy(sub_table+(i-offset)*lva2pa_size, lvas+i, sizeof(clmpthash_lva));
            memcpy(sub_table+(i-offset)*lva2pa_size+sizeof(clmpthash_lva), pas+i, sizeof(clmpthash_physical_addr));
        }

        sub_table_addr[buf_id++]=(uint64_t)sub_table;
        offset=last;
    }

    return inner_index;
}