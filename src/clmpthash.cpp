#include "lmpthash.hpp"
#include "clmpthash.h"
#include "include/utils/logger.hpp"

void cparse_msr_cambridge(LVA** lvas, PhysicalAddr** pas, uint64_t* num_lva, LVA** querys, uint64_t * num_querys, std::string trace_path){
    std::vector<LVA> uniq_lpn;
    std::vector<LVA> lpns;
    parse_MSR_Cambridge(uniq_lpn, lpns, trace_path);
    *num_lva=uniq_lpn.size();
    *num_querys=lpns.size();

    *lvas=(LVA*)malloc(sizeof(LVA)*uniq_lpn.size());
    *querys=(LVA*)malloc(sizeof(LVA)*lpns.size());

    for(uint64_t i=0;i<uniq_lpn.size();i++){
        (*lvas)[i]=uniq_lpn[i];
    }
    for(uint64_t i=0;i<lpns.size();i++){
        (*querys)[i]=lpns[i];
    }
}

void parse_configuration(char* config_path, clmpthash_config* cfg, LVA** lvas, PhysicalAddr** pas, uint64_t* num_lva, LVA** querys, uint64_t* num_querys){
    lmpthash_config config;
    config.load_config(std::string(config_path));

    // get trace file
    cparse_msr_cambridge(lvas, pas, num_lva, querys, num_querys, config.trace_path);
    printf("    lva_num: %lu, query_num: %lu\n", *num_lva, *num_querys);
    *pas=(PhysicalAddr*)malloc(sizeof(PhysicalAddr)*(*num_lva));
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

void* build_index(LVA* lvas, PhysicalAddr* pas, uint64_t num, clmpthash_config* cfg){
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

    LMPTHashBuilder<LVA, PhysicalAddr>* builder = new LMPTHashBuilder<LVA, PhysicalAddr>(config);
    std::vector<LVA> keys(lvas, lvas+num);
    std::vector<PhysicalAddr> values(pas, pas+num);

    builder->Segmenting(keys);
    builder->Learning();
    builder->Multi_Bucketing();
    builder->Tabling(keys, values);
    printf("Build index done\n");

    return static_cast<void*>(builder);
}

int get_pa(LVA lva, void* index, PhysicalAddr* pa){
    LMPTHashBuilder<LVA, PhysicalAddr>* builder=static_cast<LMPTHashBuilder<LVA, PhysicalAddr>*>(index);
    
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

    PhysicalAddr* subtable=(PhysicalAddr*)(builder->lmpt_segments[ans].next_addr);
    *pa=subtable[pos];

    // output pa
    for(int i=0;i<20;i++){
        slogger.func_log(4, "%02x ", (*pa).data[i]);
    }
    slogger.func_log(4, "\n");

    return 0;
}

int clean_index(void* index){
    LMPTHashBuilder<LVA, PhysicalAddr>* builder=static_cast<LMPTHashBuilder<LVA, PhysicalAddr>*>(index);

    builder->Cleaning();
    // printf("finish cleanning buffer\n");
    return 0;
}