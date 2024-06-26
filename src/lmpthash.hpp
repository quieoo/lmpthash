#include <vector>
#include <queue>
#include <iostream>
#include <vector>
#include "include/pthash.hpp"
#include "include/utils/logger.hpp"
#include "pgm/pgm_index.hpp"
#include <thread>
#include <mutex>
#include <atomic>
#include <condition_variable>
#include <unordered_set>

#include "include/clmpthash.h"

const int page_size = 4096;
typedef pthash::single_phf<pthash::murmurhash2_64, pthash::compact, false> pthash_type;

#define ANSI_CURSOR_UP(n)    "\033[" #n "A"

template <typename Key>
struct MonoSegment{
    std::vector<Key> first_key;
    std::vector<Key> last_key;

    uint64_t NumKeys(){
        uint64_t ret=0;
        for(size_t i=0;i<first_key.size();i++){
            ret+=last_key[i]-first_key[i]+1;
        }
        return ret;
    }

    uint64_t NumEmpty(){
        uint64_t rt=0;
        for(size_t i=0;i<first_key.size()-1;i++){
            rt+=first_key[i+1]-last_key[i];
        }
        return rt;
    }

    uint64_t DistanceTo(MonoSegment& next){
        int64_t ret;
        ret=next.first_key[0]-last_key.back();
        if(ret<0){
            printf("eror DistanceTo: %ld-%ld\n", next.first_key[0], last_key.back());
            return -1;
        }
        return ret;
    }

    void MergeWith(MonoSegment& next){
        for(size_t i=0;i<next.first_key.size();i++){
            first_key.push_back(next.first_key[i]);
            last_key.push_back(next.last_key[i]);
        }
        // clean next
        next.first_key.resize(0);
        next.last_key.resize(0);
    }
};

template <typename Key>
struct MonoSegmentMerger{
    float alpha;
    float beta;
    float gamma;
    uint32_t P;
    std::vector<MonoSegment<Key>> segments;

    MonoSegmentMerger() = default;
    MonoSegmentMerger(float _alpha, float _beta, float _gamma, uint32_t _P) :
        alpha(_alpha), beta(_beta), gamma(_gamma), P(_P){
            // printf("alpha: %f, beta: %f, gamma: %f, P: %d\n", alpha, beta, gamma, P);
    }
    void LoadKeys(std::vector<Key>& keys){
        // make sure keys are sorted
        assert(std::is_sorted(keys.begin(), keys.end()));
        // get all continuous keys and insert into segments
        auto it = keys.begin();
        while (it != keys.end()) {
            auto start = it;
            while (std::next(it) != keys.end() && *it + 1 == *std::next(it)) {
                ++it;
            }
            segments.emplace_back();
            segments.back().first_key.push_back(*start);
            segments.back().last_key.push_back(*it);
            ++it;
        }
        printf("    load keys: %lu, get %lu monotonical segments\n", keys.size(), segments.size());
    }

    void mockkeys(){
        //generate keys and insert into segments
        std::vector<Key> keys;
        
        keys.push_back(0);
        keys.push_back(100);
        keys.push_back(2);
        keys.push_back(4);
        keys.push_back(5);
        keys.push_back(12);
        keys.push_back(16);
        keys.push_back(17);
        LoadKeys(keys);
    }

    uint64_t GetScore(uint32_t seg_id, uint32_t next_seg_id){
        if(seg_id<0 || seg_id>=segments.size()-1){
            printf("invalid seg_id: %d\n", seg_id);
            return -1;
        }
        int num_sparse_seg=0;
        if(segments[seg_id].first_key.size()>1){
            num_sparse_seg++;
        }
        if(segments[next_seg_id].first_key.size()>1){
            num_sparse_seg++;
        }
        double ret=(double)alpha*(segments[seg_id].NumKeys()+segments[next_seg_id].NumKeys()) + (double)beta*(segments[seg_id].DistanceTo(segments[next_seg_id]));
        /*
        seg_id: 1715, next_seg_id: 1716
        */
        for(int i=0;i<num_sparse_seg;i++){
            ret*=gamma;
        }
        // printf("seg-%d: score: %f, num_sparse_seg: %d\n", seg_id, ret, num_sparse_seg);
        
        return ret;
    }

    struct pq_item{
        uint64_t score;
        uint32_t seg_id;
        uint32_t ver;

        pq_item(uint64_t _score, uint32_t _seg_id, uint32_t _ver) : score(_score), seg_id(_seg_id), ver(_ver){}

        bool operator<(const pq_item& rhs) const {
            return score>rhs.score;
        }
    };

    struct dlist_item{
        uint32_t next;
        uint32_t prev;
        uint32_t ver;

        dlist_item(uint32_t _next, uint32_t _prev, uint32_t _ver) : next(_next), prev(_prev), ver(_ver){}
    };



    void GreedyMerge(){
        const uint32_t N=segments.size();
        printf("    Current number of segments: %u, merge to (P=%u) segments\n", N, P);
        if(N<=P){
            return;
        }
        const uint32_t NumMerge=N-P;
        const uint32_t NoEdge=UINT32_MAX;
        
        // init scores and create a min-heap
        std::vector<pq_item> pq_items;
        for(uint32_t i=0; i<N-1; ++i){
            pq_items.push_back(pq_item(GetScore(i, i+1), i, 0));
        }
        std::priority_queue<pq_item> pq(pq_items.begin(), pq_items.end());

        // create a doubly linked list on segments
        std::vector<dlist_item> dlist(N, dlist_item(NoEdge, NoEdge, 0));
        for(uint32_t i=1; i<N-1; ++i){
            dlist[i].next=i+1;
            dlist[i].prev=i-1;
        }
        dlist[0].next=1;
        dlist[N-1].prev=N-2;
                
        uint32_t cnt=0;
        // start merging
        while(!pq.empty() && cnt<NumMerge){
            auto top = pq.top();
            pq.pop();
            uint32_t seg_id=top.seg_id;
            // printf("%u-%u, get seg from top: seg-%d, score: %lu\n",cnt, NumMerge, seg_id, top.score);
            // output status and next_prev of segments
            // printf("    status: ");
            // for(uint32_t i=0; i<N; ++i){
            //     printf("%d-(%d, %d) ", dlist[i].ver, dlist[i].next, dlist[i].prev);
            // }
            // printf("\n");
            uint32_t next_seg_id=dlist[seg_id].next;
            uint32_t previous_seg_id=dlist[seg_id].prev;

            if(top.ver==dlist[seg_id].ver && next_seg_id!=NoEdge){
                // printf("    merge seg-%d and seg-%d\n", seg_id, next_seg_id);
                // printf("        top: version-%d. Dlist: version-%d, next-%d, prev-%d\n", top.ver, dlist[seg_id].ver, dlist[seg_id].next, dlist[seg_id].prev);
                // merge current segment with next segment
                segments[seg_id].MergeWith(segments[next_seg_id]);

                // update dlist
                dlist[seg_id].ver=dlist[seg_id].ver+1;
                dlist[next_seg_id].ver=dlist[next_seg_id].ver+1;
                dlist[seg_id].next=dlist[next_seg_id].next;
                if(dlist[next_seg_id].next!=NoEdge)
                    dlist[dlist[next_seg_id].next].prev=seg_id;           
                // if current and merged segment have next segment(not the tail), it can still be used to merge so push it back to pq
                if(dlist[seg_id].next!=NoEdge){
                    uint64_t new_score=GetScore(seg_id, dlist[seg_id].next);
                    pq_item new_item(new_score, seg_id, dlist[seg_id].ver);
                    pq.push(new_item);
                }
                
                // the score of pervious segment should be also updated
                if(previous_seg_id!=NoEdge){
                    dlist[previous_seg_id].ver=dlist[previous_seg_id].ver+1;
                    uint64_t new_score=GetScore(previous_seg_id, seg_id);
                    pq_item new_item(new_score, previous_seg_id, dlist[previous_seg_id].ver);
                    pq.push(new_item);
                }
                cnt++;

            }
        }

        // remove invalid segments
        std::vector<MonoSegment<Key>> new_segs;
        for(uint32_t i=0; i<N; ++i){
            if(segments[i].first_key.size()>0){
                // invalid segment
                new_segs.push_back(segments[i]);
            }
        }
        segments=new_segs;
    }

    void PrintSeg(uint32_t seg_id){
        MonoSegment<Key>& seg=segments[seg_id];
        printf("seg-%ld ", seg_id);
        for(size_t i=0;i<seg.first_key.size();i++){
            printf("%ld-%ld ",seg.first_key[i], seg.last_key[i]);
        }
        printf("\n");
    }

    void PrintSegs(){
        for (size_t s=0;s<segments.size();s++){
            printf("seg-%ld ", s);
            MonoSegment<Key>& seg=segments[s];
            for(size_t i=0;i<seg.first_key.size();i++){
                printf("%ld-%ld ",seg.first_key[i], seg.last_key[i]);
            }
            if(s<segments.size()-1){
                printf(" score: %ld", GetScore(s, s+1));
            }
            printf("\n");
        }
    }

    void ScoreSegs(){
        uint64_t num_keys_dense=0;
        uint64_t num_keys_sparse=0;
        uint64_t num_slots_sparse=0;
        uint64_t max_dense_seg=0;
        uint64_t max_sparse_seg=0;
        uint64_t num_denseg_segs=0;
        uint64_t num_sparse_segs=0;
        for (size_t s=0;s<segments.size();s++){
            MonoSegment<Key>& seg=segments[s];
            if(seg.first_key.size()==1){
                num_denseg_segs++;
                num_keys_dense+=seg.NumKeys();
                if(seg.NumKeys()>max_dense_seg) max_dense_seg=seg.NumKeys();
            }
            else{
                num_sparse_segs++;
                num_keys_sparse+=seg.NumKeys();
                num_slots_sparse+=seg.NumEmpty();
                if(seg.NumKeys()>max_sparse_seg) max_sparse_seg=seg.NumKeys();
            }
        }
        
        printf("    num_denseg_segs: %ld, num_sparse_segs: %ld, num_keys_dense: %ld, num_keys_sparse: %ld, num_slots_sparse: %ld, max_dense_seg: %ld, max_sparse_seg: %ld\n", num_denseg_segs, num_sparse_segs , num_keys_dense, num_keys_sparse, num_slots_sparse, max_dense_seg, max_sparse_seg);

    }
};


struct lmpthash_config{
    lmpthash_config():
        alpha(0.5),
        beta(0.5),
        gamma(0.5),
        P(65536*9/10),
        hashed_num_bucket_c(6.0),
        table_size_alpha(0.94),
        max_bucket_size(100),
        pilot_search_threshold(10000000){}

    void load_config(std::string config_file){
        std::ifstream fin(config_file);
        std::string line;
        std::string name;
        std::string value;
        while(getline(fin, line)){
            std::istringstream sin(line);
            sin>>name>>value;
            if(name=="alpha") alpha=std::stof(value);
            else if(name=="beta") beta=std::stof(value);
            else if(name=="gamma") gamma=std::stof(value);
            else if(name=="P") P=std::stoi(value);
            else if(name=="hashed_num_bucket_c") hashed_num_bucket_c=std::stof(value);
            else if(name=="table_size_alpha") table_size_alpha=std::stof(value);
            else if(name=="max_bucket_size") max_bucket_size=std::stoi(value);
            else if(name=="pilot_search_threshold") pilot_search_threshold=std::stoi(value);
            else if(name=="dynamic_alpha") dynamic_alpha=std::stoi(value);
            else if(name=="alpha_limits") alpha_limits=std::stof(value);
            else if(name=="left_epsilon") left_epsilon=std::stoi(value);
            else if(name=="right_epsilon") right_epsilon=std::stoi(value);
            else if(name=="trace_path") trace_path=value;
        }

        std::cout<<"================="<<config_file<<"=================="<<std::endl;
        std::cout<<"alpha: "<<alpha<<std::endl;
        std::cout<<"beta: "<<beta<<std::endl;
        std::cout<<"gamma: "<<gamma<<std::endl;
        std::cout<<"P: "<<P<<std::endl;
        std::cout<<"hashed_num_bucket_c: "<<hashed_num_bucket_c<<std::endl;
        std::cout<<"table_size_alpha: "<<table_size_alpha<<std::endl;
        std::cout<<"max_bucket_size: "<<max_bucket_size<<std::endl;
        std::cout<<"pilot_search_threshold: "<<pilot_search_threshold<<std::endl;
        std::cout<<"dynamic_alpha: "<<dynamic_alpha<<std::endl;
        std::cout<<"alpha_limits: "<<alpha_limits<<std::endl;
        std::cout<<"left_epsilon: "<<left_epsilon<<std::endl;
        std::cout<<"right_epsilon: "<<right_epsilon<<std::endl;
        std::cout<<"trace_path: "<<trace_path<<std::endl;
    }


    // segmentation parameters
    float alpha;
    float beta;
    float gamma;
    uint32_t P;

    // bucketing parameters
    double hashed_num_bucket_c;
    double table_size_alpha;
    uint32_t max_bucket_size;
    uint32_t pilot_search_threshold;
    bool dynamic_alpha;
    double alpha_limits;

    // Learned Index parameters
    int left_epsilon;
    int right_epsilon;

    std::string trace_path;
};

struct LMPTSegment{
    uint64_t first_key;
    uint32_t slope;
    uint8_t seg_type;
    uint64_t next_addr=0;
    uint32_t size=0;
};


template <typename Key, typename Value>
struct LMPTHashBuilder{
    std::vector<LMPTSegment> lmpt_segments;
    MonoSegmentMerger<Key> ms_merger;
    lmpthash_config cfg;
    pthash::simple_logger slogger;
    std::unordered_map<uint32_t, pthash_type> pthash_map;
    pgm::PGMIndex<Key, 64,4,uint32_t> pgm_index;
    int epsilon;

    LMPTHashBuilder(lmpthash_config _cfg){
        cfg=_cfg;
        ms_merger=MonoSegmentMerger<Key>(cfg.alpha, cfg.beta, cfg.gamma, cfg.P);
        slogger.log_level=1;
    }

    int Segmenting(std::vector<Key>& keys){
        printf("## Segmentating ##\n");
        printf("    # alpha: %f, beta: %f, gamma: %f, P: %d\n", cfg.alpha, cfg.beta, cfg.gamma, cfg.P);
        // merge segments
        ms_merger.LoadKeys(keys);
        ms_merger.GreedyMerge();
        ms_merger.ScoreSegs();

        lmpt_segments.resize(ms_merger.segments.size());
        for(uint32_t i=0; i<ms_merger.segments.size(); ++i){
            MonoSegment<Key>& seg=ms_merger.segments[i];
            lmpt_segments[i].first_key=seg.first_key[0];
            // printf("    seg %d, first_key: %llu\n", i, lmpt_segments[i].first_key);
        }
        return 0;
    }


    // 结构体用于保存线程处理结果
    struct SegmentResult {
        uint32_t index;
        uint32_t slope;
        pthash_type pthashMap;
    };

    // 全局变量，用于保存所有线程处理结果
    std::vector<SegmentResult> segmentResults;
    std::atomic<int> currentSegment;

    // 处理单个段的函数
    void processSegment(const std::vector<MonoSegment<Key>>& segments, int index, pthash::build_configuration config) {
        MonoSegment seg = segments[index];
        SegmentResult result;
        result.index = index;
        
        if(seg.first_key.size() == 1) {
            // accurate segment
            result.slope = 0;
        } else {
            std::vector<uint32_t> keys;
            for(uint32_t j = 0; j < seg.first_key.size(); ++j) {
                for(uint32_t k = seg.first_key[j]; k <= seg.last_key[j]; ++k) {
                    keys.push_back(k);
                }
            }
            
            config.LinearMapping = true;
            pthash_type f;
            auto ret = f.build_in_internal_memory(keys.begin(), keys.size(), config);
            if(ret.searching_seconds == -1) {
                // printf("seg_%d build_in_internal_memory hash\n", index);
                config.LinearMapping = false;
                ret = f.build_in_internal_memory(keys.begin(), keys.size(), config);
                if(ret.searching_seconds == -1) {
                    printf("build_in_internal_memory hashed failed\n");
                }
            }
            
            result.slope = f.get_slope();
            result.pthashMap = f;
        }
        
        // 保存处理结果到全局变量
        segmentResults[index] = result;
    }

    int Multi_Bucketing() {
        printf("## Bucketing ##\n");
        printf("    # hashed_num_bucket_c: %f, table_size_alpha: %f, max_bucket_size: %d, pilot_search_threshold: %d, dynamic_alpha: %d, alpha_limits: %f\n", cfg.hashed_num_bucket_c, cfg.table_size_alpha, cfg.max_bucket_size, cfg.pilot_search_threshold, cfg.dynamic_alpha, cfg.alpha_limits);
        
        uint32_t seg_num = ms_merger.segments.size();
        pthash::build_configuration config;
        config.c = cfg.hashed_num_bucket_c;
        config.alpha = cfg.table_size_alpha;
        config.minimal_output = true;
        config.verbose_output = false;
        config.seed=0x123456789;
        config.LinearMapping = true;
        config.max_bucket_size = cfg.max_bucket_size;
        config.pilot_search_threshold = cfg.pilot_search_threshold;
        config.dynamic_alpha = cfg.dynamic_alpha;
        config.alpha_limits = cfg.alpha_limits;
        
        
        segmentResults.resize(seg_num);
        currentSegment.store(0);
        std::vector<std::thread> threads;
        int numThreads = std::thread::hardware_concurrency(); 
        printf("    Progress with %d threads: \n", numThreads);
      
        for(int i = 0; i < numThreads; ++i) {
            threads.emplace_back([&](){
                while (true) {
                    int segmentIndex = currentSegment.fetch_add(1);
                    if (segmentIndex >= seg_num) break; // 所有段都处理完毕
                    processSegment(ms_merger.segments, segmentIndex, config);
                    printf("\r    %d / %d\n", segmentIndex, seg_num);
                    std::cout<<ANSI_CURSOR_UP(1);
                }
            });
        }

        for(auto& thread : threads) {
            thread.join();
        }

        for(uint32_t i = 0; i < seg_num; ++i) {
            // printf("\r    %d / %d\n", i, seg_num);
            // std::cout << ANSI_CURSOR_UP(1);
            lmpt_segments[i].slope=segmentResults[i].slope;
            pthash_map[i]=segmentResults[i].pthashMap;
        }
        
        return 0;
    }

    int Bucketing(){
        printf("## Bucketing ##\n");
        printf("    # hashed_num_bucket_c: %f, table_size_alpha: %f, max_bucket_size: %d, pilot_search_threshold: %d, dynamic_alpha: %d, alpha_limits: %f\n", cfg.hashed_num_bucket_c, cfg.table_size_alpha, cfg.max_bucket_size, cfg.pilot_search_threshold, cfg.dynamic_alpha, cfg.alpha_limits);
        uint32_t seg_num=ms_merger.segments.size();

        pthash::build_configuration config;
        config.c = cfg.hashed_num_bucket_c;
        config.alpha = cfg.table_size_alpha;
        config.minimal_output = true;  // mphf
        config.verbose_output = false;
        config.seed=0x123456789;

        config.LinearMapping=true;
        config.max_bucket_size = cfg.max_bucket_size;
        config.pilot_search_threshold = cfg.pilot_search_threshold;

        config.dynamic_alpha = cfg.dynamic_alpha;
        config.alpha_limits = cfg.alpha_limits;
        printf("    Progress: \n");
        for(uint32_t i=0; i<seg_num; ++i){
            printf("\r    %d / %d\n", i, seg_num);
            std::cout<<ANSI_CURSOR_UP(1);
            MonoSegment<Key>& seg=ms_merger.segments[i];
            if(seg.first_key.size()==1){
                // accurate segment
                lmpt_segments[i].slope=0;
            }else{
                // printf("Segment-%d ", i);
                // approximate segment
                // get keys
                std::vector<Key> keys;
                for(uint32_t j=0; j<seg.first_key.size(); ++j){
                    for(uint32_t k=seg.first_key[j]; k<=seg.last_key[j]; ++k){
                        keys.push_back(k);
                    }
                }
                // try to create a linear mapping pthash
                config.LinearMapping=true;
                pthash_type f;
                auto ret=f.build_in_internal_memory(keys.begin(), keys.size(), config);
                if(ret.searching_seconds==-1){
                    // printf("Failed to create a linear mapping pthash, transfer to hash mapping \n");
                    config.LinearMapping=false;
                    ret=f.build_in_internal_memory(keys.begin(), keys.size(), config);
                }
                uint64_t _slope=f.get_slope();
                if(_slope >= UINT32_MAX){
                    printf("Error Bucketing: slope is larger than UINT32_MAX\n");
                    return -1;
                }
                lmpt_segments[i].slope=_slope;
                pthash_map[i]=f;
                // printf("slope: %ld\n", lmpt_segments[i].slope);
            }
        }

        printf("\n");
        return 0;
    }

    int Tabling(std::vector<Key>& keys, std::vector<Value>& values){
        printf("## Tabling ##\n");
        // make sure vector keys are sorted
        assert(std::is_sorted(keys.begin(), keys.end()));
        // divide keys into segments according to first_key in lmpt_segments
        uint32_t seg_id=0;
        uint64_t l=0,r=0;
        while(l<keys.size()){
            Key last_key=seg_id<(lmpt_segments.size()-1)?lmpt_segments[seg_id+1].first_key:UINT64_MAX;
            while(r<keys.size() && keys[r]<last_key){
                ++r;
            }
            slogger.func_log(0, "Segment-%d. num_of keys: %ld, %ld - %ld = %ld", seg_id, r-l, keys[r], keys[l], keys[r]-keys[l]);            
            if(lmpt_segments[seg_id].slope==0){
                // accurate segment
                slogger.func_log(0, ", [0] table_size: %ld\n", r-l);
                // Value* sub_table=(Value*)malloc((r-l)*sizeof(Value));
                Value* sub_table=new Value[r-l];

                for(uint64_t i=l; i<r; ++i){
                    Key k=keys[i];
                    uint64_t pos=k-lmpt_segments[seg_id].first_key;
                    sub_table[pos]=values[i];
                    // if(pos==0){
                    //     printf("seg_id: %d, key: %ld", seg_id, k);
                    //     values[i].output();
                    //     sub_table[pos].output();
                    //     printf("\n");
                    // }
                }
                lmpt_segments[seg_id].next_addr=(uint64_t)sub_table;
                lmpt_segments[seg_id].size=r-l;
                lmpt_segments[seg_id].seg_type=0;
            }else{
                // approximate segment
                uint64_t table_size=pthash_map[seg_id].table_size();
                slogger.func_log(0, ", [%u] table_size: %ld\n",2-(pthash_map[seg_id].is_linear_mapping()), table_size); 
                // Value* sub_table=(Value*)malloc(table_size*sizeof(Value));
                Value* sub_table=new Value[table_size];
                std::unordered_set<Key> keys_set;
                for(uint64_t i=l; i<r; ++i){
                    Key k=keys[i];
                    uint64_t pos=pthash_map[seg_id](k);

                    //check if pos is identical
                    if(keys_set.find(pos)!=keys_set.end()){
                        printf("error: pos is identical\n");
                        return -1;
                    }
                    keys_set.insert(pos);

                    sub_table[pos]=values[i];
                }
                lmpt_segments[seg_id].next_addr=(uint64_t)sub_table;
                lmpt_segments[seg_id].size=table_size;
                lmpt_segments[seg_id].seg_type=2-(pthash_map[seg_id].is_linear_mapping());
            }
            // printf("seg-%d, first_key-%lu, seg_type-%u\n", seg_id, lmpt_segments[seg_id].first_key, lmpt_segments[seg_id].seg_type);
            ++seg_id;
            l=r;
        }
        int seg_type[3]={0,0,0};
        for(int i=0; i<lmpt_segments.size(); ++i){
            if(lmpt_segments[i].seg_type==0){
                seg_type[0]++;
            }else if(lmpt_segments[i].seg_type==1){
                seg_type[1]++;
            }else if(lmpt_segments[i].seg_type==2){
                seg_type[2]++;
            }
        }
        printf("    # seg_type: %ld, %ld, %ld\n", seg_type[0], seg_type[1], seg_type[2]);
        return 0;
    }

    int Cleaning(){
        printf("## Cleaning ##\n");
        uint32_t seg_num=lmpt_segments.size();
        printf("    # segment table size: %ld\n", seg_num);
        for(uint32_t i=0; i<seg_num; i++){
            if(i==seg_num){
                // printf("Break\n");
                break;
            }
            // printf("    seg-%u/%u: %p\n", i,seg_num, lmpt_segments[i].next_addr);
            if(lmpt_segments[i].next_addr!=0){
                // free((void*)(lmpt_segments[i].next_addr));
                Value* t=(Value*)(lmpt_segments[i].next_addr);
                delete[] t;
            }
        }
        return 0;
    }


    int Learning(){
        printf("## Inner Indexing ##\n");
        printf("    # segment table size: %ld, left_epsilon: %d, right_epsilon: %d\n", lmpt_segments.size(), cfg.left_epsilon, cfg.right_epsilon);
        std::vector<Key> segment_offsets;
        for(uint32_t i=0; i<lmpt_segments.size(); i++){
            segment_offsets.push_back(lmpt_segments[i].first_key);
        }
        // make sure segment_offsets are sorted
        assert(std::is_sorted(segment_offsets.begin(), segment_offsets.end()));

        // binary search for a minimum PGM height
        slogger.func_log(0, "    binary search for a minimum PGM height\n");
        pgm::PGMIndex<Key, 64,4,uint32_t> pgm(segment_offsets, cfg.right_epsilon, cfg.right_epsilon);

        int min_height=pgm.height();
        // int l=cfg.left_epsilon, r=cfg.right_epsilon;
        // while(l<=r){
        //     int mid=(l+r)/2;
        //     pgm::PGMIndex<Key, 64,4,uint32_t> pgm(segment_offsets, mid, mid);
        //     int height=pgm.height();
        //     slogger.func_log(0, "    mid: %d, height: %d\n", mid, height);
        //     if(height<min_height){
        //         min_height=height;
        //         l=mid+1;
        //     }else{
        //         r=mid-1;
        //     }
        // }

        // binary search for a minimum epsilon that generates the minimum height
        slogger.func_log(0, "    binary search for a minimum epsilon that generates the minimum height\n");
        int l=cfg.left_epsilon, r=cfg.right_epsilon;
        int ep=-1;
        while(l<=r){
            int mid=(l+r)/2;
            pgm::PGMIndex<Key, 64, 4, uint32_t> pgm(segment_offsets, mid, mid);
            int height=pgm.height();
            slogger.func_log(0, "    mid: %d, height: %d\n", mid, height);
            if(height==min_height){
                ep=mid;
                r=mid-1;
                pgm_index=pgm;
            }else{
                l=mid+1;
            }
        }
        printf("    min_height: %d, min_epsilon: %d\n", min_height, ep);
        if(ep==-1){
            printf("        WARNING: can't find a valid epsilon\n");
            return -1;
        }
        // pgm_index.output_levels();
        // pgm_index.output_segments();
        epsilon=ep;
        return 0;
    }

    int Verifing(std::vector<Key>& keys, std::vector<Value>& values){
        printf("## Verifing ##\n");
        for(uint64_t i=0; i<keys.size(); i++){
            slogger.func_log(2, "    %ld-th key: %lu\n", i, keys[i]);
            Key k=keys[i];
            Value v=values[i];
            // k=1354511;
            
            auto range=pgm_index.search(k);
            auto lo=range.lo;
            auto hi=range.hi;
            if(hi >= lmpt_segments.size()) hi=lmpt_segments.size()-1;
            int ans=-1;
            // binary search for a maximum segment whose first_key <= k
            while(lo<=hi){
                auto mid=(lo+hi)/2;
                slogger.func_log(2, " seg id: %ld, first_key: %ld, k: %ld\n",mid, lmpt_segments[mid].first_key, k);
                if(lmpt_segments[mid].first_key>k){
                    hi=mid-1;
                }else{
                    lo=mid+1;
                    ans=mid;
                }
            }
            if(ans==-1){
                printf("error: segment not found\n");
                return 1;
            }
            slogger.func_log(2, "    found seg_id: %ld\n", ans);
            int64_t pos=-1;
            if(lmpt_segments[ans].seg_type==0){
                pos=k-lmpt_segments[ans].first_key;
                slogger.func_log(2, "    [%d] pos: %u\n",lmpt_segments[ans].seg_type, pos);
            }else if(lmpt_segments[ans].seg_type==1 || lmpt_segments[ans].seg_type==2){
                pos=pthash_map[ans](k);
            }

            if(pos!=-1){
                Value* t=(Value*)(lmpt_segments[ans].next_addr);
                Key _lva=0;
                for(int j=0;j<8;j++){
                    // _lva=(_lva<<8)+t[pos].data[j];
                    // _lva=(_lva<<8)+v.data[j];
                    _lva<<=8;
                    _lva+=v.data[j];
                }

                if(_lva != k){
                    printf("%d-th key is wrong, should be %lx but got %lx\n", i, k, _lva);
                    return 1;
                }

                if(t[pos]!=v){
                    printf("%d-th key is wrong\n", i);
                    printf("error: wrong value\n");
                    printf("shoud be: ");
                    v.output();
                    printf(" but got:");
                    t[pos].output();
                    printf("\n");
                    return 1;
                }
            }else{
                printf("%d-th key is wrong\n", i);
                printf("wrong type\n");
                return 1;
            }

            // break;
        }
        printf("    all keys are correct\n");
        return 0;
    }

    int Compacting(void* ptr){
        printf("## Compacting ##\n");
        uint32_t inner_index_size=32+pgm_index.segments.size()*sizeof(clmpthash_pgm_segment)+(lmpt_segments.size()+1)/2*2*sizeof(uint64_t)+lmpt_segments.size()*16;
        printf("    # Inner Index Size: %f MB\n", inner_index_size/1024.0/1024.0);
        if(inner_index_size > 1024*1024){
            printf("Error while Compacting: inner_index_size > 1MB. %f MB\n", inner_index_size/1024.0/1024.0);
            return -1;
        }

        // first_key and last_key takes the first 16 bytes
        uint64_t* p64=(uint64_t*)ptr;
        p64[0]=lmpt_segments[0].first_key;
        ptr+=16;

        // num level of pgm_index take the next 2 bytes
        uint16_t* p16=(uint16_t*)ptr;
        p16[0]=pgm_index.height();
        if(p16[0]>3){
            printf("Error while Compacting: pgm height > 3\n");
            return -1;
        }
        ptr+=2;

        // variable_epsilon_value takes the next 2 bytes
        p16=(uint16_t*)ptr;
        p16[0]=pgm_index.variable_epsilon_value;
        ptr+=2;
        // array of level_offsets takes the next 12 bytes, each one takes 2 bytes
        p16=(uint16_t*)ptr;
        std::vector<size_t> lofs;
        pgm_index.get_level_offsets(lofs);
        
        int num_level_offsets=lofs.size();
        for(int i=0;i<num_level_offsets;i++){
            p16[i]=lofs[i]+2;
        }
        // p16[num_level_offsets]=p16[num_level_offsets-1]+(lmpt_segments.size()+1)/2;
        p16[num_level_offsets]=p16[num_level_offsets-1]+(lmpt_segments.size()+2)/2; // add the boundary to first_key
        p16[num_level_offsets+1]=p16[num_level_offsets]+lmpt_segments.size();

        printf("    level_offsets: ");
        for(int i=0; i<=num_level_offsets+1;i++){
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
            segs.push_back(ps);
        }
        memcpy(ptr, segs.data(), segs.size()*sizeof(clmpthash_pgm_segment));
        ptr+=segs.size()*sizeof(clmpthash_pgm_segment);

        // each first_key of clmpthash_htl_segment takes 8 bytes
        p64=(uint64_t*)ptr;
        for(int i=0; i<lmpt_segments.size();i++){
            p64[i]=lmpt_segments[i].first_key;
        }
        // add a max bound to the end of first_key
        p64[lmpt_segments.size()]=UINT64_MAX;
        // round up to 16 bytes
        ptr+=(lmpt_segments.size()+2)/2*2*sizeof(uint64_t);

        // each clmpthash_htl_segment takes 16 bytes
        p64=(uint64_t*)ptr;
        for(int i=0; i<lmpt_segments.size();i++){
            if(lmpt_segments[i].seg_type==0){
                // accurate segment
                // allocate table data
                uint64_t table_size=(lmpt_segments[i].size)*sizeof(Value)+sizeof(uint64_t);
                uint8_t* raw_table=new uint8_t[table_size];
                p64[i*2]=(uint64_t)raw_table;
                memcpy(raw_table, &table_size, sizeof(uint64_t));
                memcpy(raw_table+sizeof(uint64_t), (void*)(lmpt_segments[i].next_addr), table_size-sizeof(uint64_t));                
            }else{
                // approximate segment
                // allocate table data
                if(pthash_map.count(i)<=0){
                    printf("Error while Compacting: pthash_map.count(i)<=0\n");
                    return -1;
                }
                uint64_t table_size=(pthash_map[i].table_size())*sizeof(Value);
                uint64_t pilot_size=pthash_map[i].pilot_bytes();

                uint64_t total_size=table_size+pilot_size+sizeof(uint64_t);
                printf("seg-%ld, table_size: %f MB, pilot_size: %ld, total_size: %ld\n", i, (double)table_size/1024/1024, pilot_size, total_size);
                uint8_t* raw_table;
                try{
                    raw_table = new uint8_t[total_size];
                }catch(const std::bad_alloc& e){
                    std::cout<<"Error while Compacting: bad_alloc"<<e.what()<<std::endl;
                    return -1;
                }
                p64[i*2]=(uint64_t)raw_table;
                // copy table data
                memcpy(raw_table, &total_size, sizeof(uint64_t));
                memcpy(raw_table+sizeof(uint64_t), (void*)(lmpt_segments[i].next_addr), table_size);
                pthash_map[i].get_pilots(raw_table+sizeof(uint64_t)+table_size);


                uint64_t pthash_meta=0;
                // fisrt 2 bits for seg_type
                pthash_meta=lmpt_segments[i].seg_type;
                pthash_meta<<=62;
                // next 30 bits of meta is slope
                if(lmpt_segments[i].slope > 1<<30){
                    printf("Error while Compacting: lmpt_segments[i].slope > 1<<30\n");
                    return -1;
                }
                pthash_meta|=(uint64_t)(lmpt_segments[i].slope)<<32;
                // next 8 bits of meta is width of pilot
                uint32_t pilot_width=pthash_map[i].get_pilot_width();
                if(pilot_width > 255){
                    printf("Error while Compacting: pilot_width > 255\n");
                    return -1;
                }
                pthash_meta|=(uint64_t)(pilot_width<<24) & 0xffffffff;
                // next 24 bits of meta table_size
                uint32_t table_size_value=pthash_map[i].table_size();
                if(table_size_value > 16777215){
                    printf("Error while Compacting: table_size > 16777215\n");
                    return -1;
                }
                pthash_meta|=(uint64_t)table_size_value & 0xffffffff;

                // check if values in pthash_meta is right
                {
                    uint8_t _seg_type=pthash_meta>>62;
                    uint32_t _slope=(pthash_meta>>32)&0x3fffffff;
                    uint32_t _pilot_width=(pthash_meta&0xff000000)>>24;
                    uint32_t _table_size=pthash_meta&0xffffff;
                    if(_seg_type!=lmpt_segments[i].seg_type){
                        printf("Error while Compacting: _seg_type!=lmpt_segments[i].seg_type\n");
                        return -1;
                    }
                    if(_slope!=lmpt_segments[i].slope){
                        printf("Error while Compacting: _slope!=lmpt_segments[i].slope\n");
                        return -1;
                    }
                    if(_pilot_width!=pilot_width){
                        printf("Error while Compacting: _pilot_width!=pilot_width\n");
                        return -1;
                    }
                    if(_table_size!=table_size_value){
                        printf("Error while Compacting: _table_size!=table_size_value\n");
                        return -1;
                    }
                }

                p64[i*2+1]=pthash_meta;
            }
        }
        return 0;
    }
};


void parse_MSR_Cambridge(std::vector<uint64_t>&  uniq_lpn, std::vector<uint64_t>& lpns, std::string filename){
    // hash table for unique lpn
    std::unordered_set<uint64_t> ht;

    // open file with name "filename", and read by lines
    std::ifstream file(filename);
    std::string line;
    uint64_t timestamp, offset, size, t0;
    char trace_name[100];
    char op[100];
    int trace_id;
    uint64_t lpn;
    while (std::getline(file, line)) {
        sscanf(line.c_str(), "%lu,%100[^,],%d,%100[^,],%lu,%lu,%lu\n", &timestamp,trace_name,&trace_id,op,&offset,&size,&t0);
        for(int i=0;i<size/page_size;i++){
            lpn=offset/page_size+i;
            lpns.push_back(lpn);
            ht.insert(lpn);
        }
    }

    // get unique lpn
    for (auto it = ht.begin(); it != ht.end(); it++) {
        uniq_lpn.push_back(*it);
    }
    // sort uniq_lpn
    std::sort(uniq_lpn.begin(), uniq_lpn.end());
    file.close();
    return;   
}

void parse_femu(std::vector<uint64_t>&  uniq_lpn, std::vector<uint64_t>& lpns, std::string filename){
    // hash table for unique lpn
    std::unordered_set<uint64_t> ht;

    // open file with name "filename", and read by lines
    std::ifstream file(filename, std::ios::binary); // 以二进制方式打开文件
    if (!file) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        return;
    }
    // 读取整个文件并存储在uint64_t数组中
    // 获取文件大小
    file.seekg(0, std::ios::end);
    std::streampos filesize = file.tellg();
    file.seekg(0, std::ios::beg);
    // printf("file size: %f MB\n", filesize/1024.0/1024.0);
    // 确定文件中包含多少个uint64_t
    size_t num_elements = filesize / sizeof(uint64_t);
    // 创建一个vector用于存储文件内容
    std::vector<uint64_t> data(num_elements);
    // 读取文件内容到vector中
    if (!file.read(reinterpret_cast<char*>(data.data()), filesize)) {
        std::cerr << "Failed to read file: " << filename << std::endl;
        return;
    }
    file.close();

    for(size_t i=0;i<data.size();i++){
        // printf("%lx\n", data[i]);
        lpns.push_back(data[i]);
        ht.insert(data[i]);
    }
    // get unique lpn
    for (auto it = ht.begin(); it != ht.end(); it++){
        uniq_lpn.push_back(*it);
    }
    // sort uniq_lpn
    std::sort(uniq_lpn.begin(), uniq_lpn.end());
    return;
}


void output_query_to_file(std::vector<uint64_t>& query, uint64_t n){
    // create output csv file
    std::ofstream outfile;
    outfile.open("query.txt");
    for(int i=0;i<n;i++){
        outfile<<query[i]<<std::endl;
    }
    outfile.close();
    return;
}