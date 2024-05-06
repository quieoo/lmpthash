#include <vector>
#include <queue>


#include <iostream>
#include <vector>


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
        next.first_key.clear();
        next.last_key.clear();
    }
};

template <typename Key>
struct MonoSegmentMerger{
    float alpha;
    float beta;
    float gamma;
    uint32_t P;
    std::vector<MonoSegment<Key>> segments;

    MonoSegmentMerger(float _alpha, float _beta, float _gamma, uint32_t _P) :
        alpha(_alpha), beta(_beta), gamma(_gamma), P(_P){
            printf("alpha: %f, beta: %f, gamma: %f, P: %d\n", alpha, beta, gamma, P);
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
        printf("load keys: %lu, get %lu segments\n", keys.size(), segments.size());
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
        printf("Current number of segments: %u, P: %u\n", N, P);
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
        
        printf("num_denseg_segs: %ld, num_sparse_segs: %ld, num_keys_dense: %ld, num_keys_sparse: %ld, num_slots_sparse: %ld, max_dense_seg: %ld, max_sparse_seg: %ld\n", num_denseg_segs, num_sparse_segs , num_keys_sparse, num_slots_sparse, max_dense_seg, max_sparse_seg);

    }
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

    LMPTHashBuilder(){
        lmpt_segments.clear();
        ms_merger=MonoSegmentMerger<Key>(0.5, 0.5, 0.5, 65536*9/10);
    }

    int Segmentation(std::vector<Key>& keys){
        
        // merge segments
        ms_merger.LoadKeys(keys);
        ms_merger.GreedyMerge();
        ms_merger.ScoreSegs();
        
        // create lmpt segments
        uint32_t seg_num=ms_merger.segments.size();
        lmpt_segments.resize(seg_num);
        for(uint32_t i=0; i<seg_num; ++i){
            MonoSegment<Key>& seg=ms_merger.segments[i];
            if(seg.first_key.size()==1){
                // accurate segment
                lmpt_segments[i].first_key=seg.first_key[0];
                lmpt_segments[i].slope=0;
            }else{
                // approximate segment
                lmpt_segments[i].first_key=seg.first_key[0];
                lmpt_segments[i].slope=1;
            }
        }
        return 0;
    }

    int Bucketing(){

    }
};