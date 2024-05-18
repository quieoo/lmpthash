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