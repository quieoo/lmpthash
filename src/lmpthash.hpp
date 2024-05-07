#include <vector>
#include <queue>
#include <iostream>
#include <vector>
#include "include/pthash.hpp"

typedef pthash::single_phf<pthash::murmurhash2_64, pthash::compact, false> pthash_type;



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
        // sort keys
        std::sort(keys.begin(), keys.end());
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
        printf("load keys: %lu, get %lu monotonical segments\n", keys.size(), segments.size());
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
        printf("Current number of segments: %u, merge to (P=%u) segments\n", N, P);
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
        
        printf("num_denseg_segs: %ld, num_sparse_segs: %ld, num_keys_dense: %ld, num_keys_sparse: %ld, num_slots_sparse: %ld, max_dense_seg: %ld, max_sparse_seg: %ld\n", num_denseg_segs, num_sparse_segs , num_keys_dense, num_keys_sparse, num_slots_sparse, max_dense_seg, max_sparse_seg);

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
        }
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
};

struct LMPTSegment{
    uint64_t first_key;
    uint32_t slope;
    uint64_t next_addr;
};


template <typename Key, typename Value>
struct LMPTHashBuilder{
    std::vector<LMPTSegment> lmpt_segments;
    MonoSegmentMerger<Key> ms_merger;
    lmpthash_config cfg;


    LMPTHashBuilder(lmpthash_config _cfg){
        cfg=_cfg;
        ms_merger=MonoSegmentMerger<Key>(cfg.alpha, cfg.beta, cfg.gamma, cfg.P);
    }

    int Segmentation(std::vector<Key>& keys){
        printf("## Segmentation ##\n");
        printf("   # alpha: %f, beta: %f, gamma: %f, P: %d\n", cfg.alpha, cfg.beta, cfg.gamma, cfg.P);
        // merge segments
        ms_merger.LoadKeys(keys);
        ms_merger.GreedyMerge();
        ms_merger.ScoreSegs();        
        return 0;
    }

    int Bucketing(){
        printf("## Bucketing ##\n");
        printf("    # hashed_num_bucket_c: %f, table_size_alpha: %f, max_bucket_size: %d, pilot_search_threshold: %d, dynamic_alpha: %d, alpha_limits: %f\n", cfg.hashed_num_bucket_c, cfg.table_size_alpha, cfg.max_bucket_size, cfg.pilot_search_threshold, cfg.dynamic_alpha, cfg.alpha_limits);
        uint32_t seg_num=ms_merger.segments.size();
        lmpt_segments.resize(seg_num);

        pthash::build_configuration config;
        config.c = cfg.hashed_num_bucket_c;
        config.alpha = cfg.table_size_alpha;
        config.minimal_output = true;  // mphf
        config.verbose_output = false;

        config.LinearMapping=true;
        config.max_bucket_size = cfg.max_bucket_size;
        config.pilot_search_threshold = cfg.pilot_search_threshold;

        config.dynamic_alpha = cfg.dynamic_alpha;
        config.alpha_limits = cfg.alpha_limits;

        for(uint32_t i=0; i<seg_num; ++i){
            MonoSegment<Key>& seg=ms_merger.segments[i];
            if(seg.first_key.size()==1){
                // accurate segment
                lmpt_segments[i].first_key=seg.first_key[0];
                lmpt_segments[i].slope=0;
            }else{
                printf("Segment-%d\n", i);
                // approximate segment
                lmpt_segments[i].first_key=seg.first_key[0];
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
                    printf("Failed to create a linear mapping pthash, transfer to hash mapping \n");
                    config.LinearMapping=false;
                    ret=f.build_in_internal_memory(keys.begin(), keys.size(), config);
                }
                // break;
            }
        }

        return 0;
    }
};