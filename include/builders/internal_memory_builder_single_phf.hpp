#pragma once

#include "include/builders/util.hpp"
#include "include/builders/search.hpp"
#include "include/utils/bucketers.hpp"
#include "include/utils/logger.hpp"
#include "include/utils/hasher.hpp"
#include <unordered_map>

namespace pthash {

template <typename Hasher>
struct internal_memory_builder_single_phf {
    typedef Hasher hasher_type;
    simple_logger sloger;


    template <typename RandomAccessIterator>
    build_timings build_from_keys(RandomAccessIterator keys, uint64_t num_keys, build_configuration const& config) {
        if(config.LinearMapping){
            std::vector<uint64_t> keys_vector(keys, keys + num_keys);
            auto ret= build_with_linear_mapping(keys_vector, num_keys, config);
            return ret;
        }
        if (config.seed == constants::invalid_seed) {
            for (auto attempt = 0; attempt < 10; ++attempt) {
                m_seed = random_value();
                try {
                    return build_from_hashes(hash_generator<RandomAccessIterator>(keys, m_seed), num_keys, config);
                } catch (seed_runtime_error const& error) {
                    std::cout << "attempt " << attempt + 1 << " failed" << std::endl;
                }
            }
            throw seed_runtime_error();
        }
        m_seed = config.seed;
        build_configuration new_config = config;
        while(1){
            auto time=build_from_hashes(hash_generator<RandomAccessIterator>(keys, m_seed), num_keys, new_config);
            if(time.searching_seconds==-1){
                new_config.alpha*=new_config.alpha;
                continue;
            }
            return time;
        }
        // return build_from_hashes(hash_generator<RandomAccessIterator>(keys, m_seed), num_keys, config);
    }

    template <typename RandomAccessIterator>
    build_timings build_with_linear_mapping(RandomAccessIterator keys, uint64_t num_keys, build_configuration const& config){
        sloger.func_log(1, "number keys: %lu\n", num_keys);
        m_num_keys = num_keys;
        m_seed=config.seed;
        
        build_timings time;
        // set to 2*num_keys to avoid over boundary
        uint64_t* temp_pilots=new uint64_t[m_num_keys*2];
        //printf("searching for a maximum bucket size with a searching threshold: %ld\n", config.pilot_search_threshold);
        // binary search for a maximum bucket size that can succssfully find out pilots
        uint64_t lo = 1, hi = config.max_bucket_size;
        m_num_buckets=0;
        /* DEBUG USE*/
        // lo=1;
        // hi=1;
        while (lo <= hi) {
            uint64_t mid=(lo+hi)/2;
            uint64_t num_bucket= num_keys/mid;
            if(num_bucket==0)    num_bucket=1;
            sloger.func_log(1, "try B: %ld, bucket number %ld\n", mid, num_bucket);  

            // get bucket ids 
            std::pair<uint64_t, uint64_t> bucketed;
            buckets_t buckets;
            {
                std::vector<pairs_t> pairs_blocks;
                bucketed=linear_map(keys, num_keys, pairs_blocks, config, num_bucket, config.max_bucket_size);
                if(bucketed.first!=0){
                    // actual used bucket
                    num_bucket=merge(pairs_blocks, buckets, config.verbose_output)+1;
                }
            }
            if(bucketed.first==0){
                // a bucket is larger than max_bucket_size
                // printf("    a bucket is larger than max_bucket_size\n");
                hi=mid-1;
                continue;
            }

            uint8_t found_in_alpha_limits=0;
            double alpha=config.alpha;
            // dynamic search for pilots
            while(1){
                uint64_t table_size = static_cast<double>(num_keys) / alpha;
                if ((table_size & (table_size - 1)) == 0) table_size += 1;
                sloger.func_log(1, "    table size: %lu, alpha: %lf\n", table_size, alpha);
                auto buckets_iterator = buckets.begin();
                memset(temp_pilots, 0, num_bucket*sizeof(uint64_t));
                // std::vector<uint64_t> temp_pilots;
                // temp_pilots.resize(num_bucket);
                // std::fill(temp_pilots.begin(), temp_pilots.end(), 0);
                int ret;
                {
                    bit_vector_builder taken(table_size);
                    uint64_t num_non_empty_buckets = buckets.num_buckets();
                    // pilots_wrapper_t pilots_wrapper(temp_pilots);
                    ret=linear_search(m_num_keys, num_bucket, num_non_empty_buckets, m_seed, config, buckets_iterator, taken, temp_pilots);
                }
                if(ret==0){
                    found_in_alpha_limits=1;
                    m_table_size=table_size;
                }else if(config.dynamic_alpha){
                    alpha*=alpha;
                    if(alpha >= config.alpha_limits){
                        continue;
                    }
                }
                break;
            }
            if(found_in_alpha_limits){
                sloger.func_log(1, "    current num_bucket works: %d\n", num_bucket);
                m_num_buckets = num_bucket;
                m_bucketer.set_divisor(bucketed.first);
                m_bucketer.set_min_key(bucketed.second);
                m_bucketer.set_num_buckets(num_bucket);
                m_pilots.resize(num_bucket);
                for(int i=0;i<num_bucket;++i){
                    m_pilots[i]=temp_pilots[i];
                }
                lo=mid+1;
            }else{
                sloger.func_log(1, "    searching threshold failed\n");
                hi=mid-1;
            }
            sloger.func_log(1, "finished searching B: %ld\n", mid);
        }
        sloger.func_log(1, "final num_bucket: %ld\n", m_num_buckets);
        // free(temp_pilots);
        delete[] temp_pilots;
        
        sloger.func_log(1, "build time: %lf seconds\n", time.searching_seconds);
        if(m_num_buckets==0){
            sloger.func_log(1, "failed to find a mimimum bucket number\n");
            time.searching_seconds=-1;
        }else{
            time.searching_seconds=0;
        }
        return time;
    }

    template <typename RandomAccessIterator>
    build_timings build_from_hashes(RandomAccessIterator hashes, uint64_t num_keys,
                                    build_configuration const& config) {
        assert(num_keys > 1);
        util::check_hash_collision_probability<Hasher>(num_keys);

        if (config.alpha == 0 or config.alpha > 1.0) {
            throw std::invalid_argument("load factor must be > 0 and <= 1.0");
        }

        clock_type::time_point start;

        start = clock_type::now();

        build_timings time;

        uint64_t table_size = static_cast<double>(num_keys) / config.alpha;
        if ((table_size & (table_size - 1)) == 0) table_size += 1;
        uint64_t num_buckets = (config.num_buckets == constants::invalid_num_buckets)
                                   ? (std::ceil((config.c * num_keys) / std::log2(num_keys)))
                                   : config.num_buckets;

        m_num_keys = num_keys;
        m_table_size = table_size;
        m_num_buckets = num_buckets;
        m_bucketer.init(m_num_buckets);
    

        if (config.verbose_output) {
            std::cout << "c = " << config.c << std::endl;
            std::cout << "alpha = " << config.alpha << std::endl;
            std::cout << "num_keys = " << num_keys << std::endl;
            std::cout << "table_size = " << table_size << std::endl;
            std::cout << "num_buckets = " << num_buckets << std::endl;
        }

        buckets_t buckets;
        {
            auto start = clock_type::now();
            std::vector<pairs_t> pairs_blocks;
            map(hashes, num_keys, pairs_blocks, config);
            auto elapsed = seconds(clock_type::now() - start);
            if (config.verbose_output) {
                std::cout << " == map+sort took: " << elapsed << " seconds" << std::endl;
            }

            start = clock_type::now();
            merge(pairs_blocks, buckets, config.verbose_output);
            elapsed = seconds(clock_type::now() - start);
            if (config.verbose_output) {
                std::cout << " == merge+check took: " << elapsed << " seconds" << std::endl;
            }
        }

        auto buckets_iterator = buckets.begin();
        time.mapping_ordering_seconds = seconds(clock_type::now() - start);
        if (config.verbose_output) {
            std::cout << " == mapping+ordering took " << time.mapping_ordering_seconds
                      << " seconds " << std::endl;
            uint64_t max_bucket_size = (*buckets_iterator).size();
            std::cout << " == max bucket size = " << max_bucket_size << std::endl;

            // avg. bucket size
            double lambda = std::log2(num_keys) / config.c;
            // avg. bucket size in first p2=b*m buckets containing p1=a*n keys
            double lambda_1 = constants::a / constants::b * lambda;
            // avg. bucket size in the other m-p2=(1-b)*m buckets containing n-p1=(1-a)*n keys
            double lambda_2 = (1 - constants::a) / (1 - constants::b) * lambda;
            std::cout << " == lambda = " << lambda << std::endl;
            std::cout << " == lambda_1 = " << lambda_1 << std::endl;
            std::cout << " == lambda_2 = " << lambda_2 << std::endl;
            buckets.print_bucket_size_distribution(max_bucket_size, num_buckets, lambda_1,
                                                   lambda_2);
        }

        start = clock_type::now();
        {
            m_pilots.resize(num_buckets);
            std::fill(m_pilots.begin(), m_pilots.end(), 0);
            bit_vector_builder taken(m_table_size);
            uint64_t num_non_empty_buckets = buckets.num_buckets();
            pilots_wrapper_t pilots_wrapper(m_pilots);
            int ret=search(m_num_keys, m_num_buckets, num_non_empty_buckets, m_seed, config, buckets_iterator, taken, pilots_wrapper);
            if (ret != 0) {
                time.searching_seconds=-1;
                return time;
            }
            if (config.minimal_output) {
                m_free_slots.clear();
                m_free_slots.reserve(taken.size() - num_keys);
                fill_free_slots(taken, num_keys, m_free_slots);
            }
        }
        time.searching_seconds = seconds(clock_type::now() - start);
        if (config.verbose_output) {
            std::cout << " == search took " << time.searching_seconds << " seconds" << std::endl;
        }

        return time;
    }

    uint64_t seed() const {
        return m_seed;
    }

    uint64_t num_keys() const {
        return m_num_keys;
    }

    uint64_t table_size() const {
        return m_table_size;
    }

    skew_bucketer bucketer() const {
        return m_bucketer;
    }

    std::vector<uint64_t> const& pilots() const {
        return m_pilots;
    }

    std::vector<uint64_t> const& free_slots() const {
        return m_free_slots;
    }

    void swap(internal_memory_builder_single_phf& other) {
        std::swap(m_seed, other.m_seed);
        std::swap(m_num_keys, other.m_num_keys);
        std::swap(m_num_buckets, other.m_num_buckets);
        std::swap(m_table_size, other.m_table_size);
        std::swap(m_bucketer, other.m_bucketer);
        m_pilots.swap(other.m_pilots);
        m_free_slots.swap(other.m_free_slots);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_seed);
        visitor.visit(m_num_keys);
        visitor.visit(m_num_buckets);
        visitor.visit(m_table_size);
        visitor.visit(m_bucketer);
        visitor.visit(m_pilots);
        visitor.visit(m_free_slots);
    }

    static uint64_t estimate_num_bytes_for_construction(uint64_t num_keys,
                                                        build_configuration const& config) {
        uint64_t table_size = static_cast<double>(num_keys) / config.alpha;
        if ((table_size & (table_size - 1)) == 0) table_size += 1;
        uint64_t num_buckets = (config.num_buckets == constants::invalid_num_buckets)
                                   ? (std::ceil((config.c * num_keys) / std::log2(num_keys)))
                                   : config.num_buckets;

        uint64_t num_bytes_for_map = num_keys * sizeof(bucket_payload_pair)          // pairs
                                     + (num_keys + num_buckets) * sizeof(uint64_t);  // buckets

        uint64_t num_bytes_for_search =
            num_buckets * sizeof(uint64_t)    // pilots
            + num_buckets * sizeof(uint64_t)  // buckets
            +
            (config.minimal_output ? (table_size - num_keys) * sizeof(uint64_t) : 0)  // free_slots
            + num_keys * sizeof(uint64_t)                                             // hashes
            + table_size / 8;  // bitmap taken

        return std::max<uint64_t>(num_bytes_for_map, num_bytes_for_search);
    }

private:
    uint64_t m_seed;
    uint64_t m_num_keys;
    uint64_t m_num_buckets;
    uint64_t m_table_size;
    skew_bucketer m_bucketer;
    std::vector<uint64_t> m_pilots;
    std::vector<uint64_t> m_free_slots;

    template <typename RandomAccessIterator>
    struct hash_generator {
        hash_generator(RandomAccessIterator keys, uint64_t seed) : m_iterator(keys), m_seed(seed) {}

        inline auto operator*() {
            return hasher_type::hash(*m_iterator, m_seed);
        }

        inline void operator++() {
            ++m_iterator;
        }

        inline hash_generator operator+(uint64_t offset) const {
            return hash_generator(m_iterator + offset, m_seed);
        }

    private:
        RandomAccessIterator m_iterator;
        uint64_t m_seed;
    };

    typedef std::vector<bucket_payload_pair> pairs_t;

    struct buckets_iterator_t {
        buckets_iterator_t(std::vector<std::vector<uint64_t>> const& buffers)
            : m_buffers_it(buffers.end() - 1), m_bucket_size(buffers.size()) {
            m_bucket.init(m_buffers_it->data(), m_bucket_size);
            skip_empty_buckets();
        }

        inline void operator++() {
            uint64_t const* begin = m_bucket.begin() + m_bucket_size;
            uint64_t const* end = m_buffers_it->data() + m_buffers_it->size();
            m_bucket.init(begin, m_bucket_size);
            if ((m_bucket.begin() - 1) == end and m_bucket_size != 0) {
                --m_bucket_size;
                --m_buffers_it;
                skip_empty_buckets();
            }
        }

        inline bucket_t operator*() const {
            return m_bucket;
        }

    private:
        std::vector<std::vector<uint64_t>>::const_iterator m_buffers_it;
        bucket_size_type m_bucket_size;
        bucket_t m_bucket;

        void skip_empty_buckets() {
            while (m_bucket_size != 0 and m_buffers_it->empty()) {
                --m_bucket_size;
                --m_buffers_it;
            }
            if (m_bucket_size != 0) m_bucket.init(m_buffers_it->data(), m_bucket_size);
        }
    };

    struct buckets_t {
        buckets_t() : m_buffers(MAX_BUCKET_SIZE), m_num_buckets(0) {}

        template <typename HashIterator>
        void add(bucket_id_type bucket_id, bucket_size_type bucket_size, HashIterator hashes) {
            assert(bucket_size > 0);
            uint64_t i = bucket_size - 1;
            m_buffers[i].push_back(bucket_id);
            for (uint64_t k = 0; k != bucket_size; ++k, ++hashes){
                m_buffers[i].push_back(*hashes);
            }

            ++m_num_buckets;
        }

        void reset(){
            m_num_buckets = 0;
            for (auto& buffer : m_buffers) {
                buffer.clear();
            }
        }

        uint64_t num_buckets() const {
            return m_num_buckets;
        };

        buckets_iterator_t begin() const {
            return buckets_iterator_t(m_buffers);
        }

        void print_bucket_size_distribution(uint64_t max_bucket_size, uint64_t num_buckets, double lambda_1, double lambda_2) {
            for (int64_t i = max_bucket_size - 1; i >= 0; --i) {
                uint64_t t = i + 1;
                uint64_t num_buckets_of_size_t = m_buffers[i].size() / (t + 1);
                uint64_t estimated_num_buckets_of_size_t =
                    (constants::b * poisson_pmf(t, lambda_1) +
                     (1 - constants::b) * poisson_pmf(t, lambda_2)) *
                    num_buckets;
                std::cout << " == num_buckets of size " << t << " = " << num_buckets_of_size_t
                          << " (estimated with Poisson = " << estimated_num_buckets_of_size_t << ")"
                          << std::endl;
            }
        }

    private:
        std::vector<std::vector<uint64_t>> m_buffers;
        uint64_t m_num_buckets;
    };

    struct pilots_wrapper_t {
        pilots_wrapper_t(std::vector<uint64_t>& pilots) : m_pilots(pilots) {}

        inline void emplace_back(bucket_id_type bucket_id, uint64_t pilot) {
            m_pilots[bucket_id] = pilot;
        }

    private:
        std::vector<uint64_t>& m_pilots;
    };

    template <typename RandomAccessIterator>
    void map_sequential(RandomAccessIterator hashes, uint64_t num_keys,
                        std::vector<pairs_t>& pairs_blocks, build_configuration const&) const {
        pairs_t pairs(num_keys);
        RandomAccessIterator begin = hashes;
        for (uint64_t i = 0; i != num_keys; ++i, ++begin) {
            auto hash = *begin;
            auto bucket_id = m_bucketer.bucket(hash.first());
            pairs[i] = {static_cast<bucket_id_type>(bucket_id), hash.second()};
        }
        std::sort(pairs.begin(), pairs.end());
        pairs_blocks.resize(1);
        pairs_blocks.front().swap(pairs);
    }

    template <typename RandomAccessIterator>
    void map_parallel(RandomAccessIterator hashes, uint64_t num_keys,
                      std::vector<pairs_t>& pairs_blocks, build_configuration const& config) const {
        pairs_blocks.resize(config.num_threads);
        uint64_t num_keys_per_thread = (num_keys + config.num_threads - 1) / config.num_threads;

        auto exe = [&](uint64_t tid) {
            auto& local_pairs = pairs_blocks[tid];
            RandomAccessIterator begin = hashes + tid * num_keys_per_thread;
            uint64_t local_num_keys = (tid != config.num_threads - 1)
                                          ? num_keys_per_thread
                                          : (num_keys - tid * num_keys_per_thread);
            local_pairs.resize(local_num_keys);

            for (uint64_t local_i = 0; local_i != local_num_keys; ++local_i, ++begin) {
                auto hash = *begin;
                auto bucket_id = m_bucketer.bucket(hash.first());
                local_pairs[local_i] = {static_cast<bucket_id_type>(bucket_id), hash.second()};
            }
            std::sort(local_pairs.begin(), local_pairs.end());
        };

        std::vector<std::thread> threads(config.num_threads);
        for (uint64_t i = 0; i != config.num_threads; ++i) threads[i] = std::thread(exe, i);
        for (auto& t : threads) {
            if (t.joinable()) t.join();
        }
    }

    template <typename RandomAccessIterator>
    void map(RandomAccessIterator hashes, uint64_t num_keys, std::vector<pairs_t>& pairs_blocks,
             build_configuration const& config) const {
        if (config.num_threads > 1 and num_keys >= config.num_threads) {
            map_parallel(hashes, num_keys, pairs_blocks, config);
        } else {
            map_sequential(hashes, num_keys, pairs_blocks, config);
        }
    }


    template <typename RandomAccessIterator>
    std::pair<uint64_t, uint64_t> linear_map(RandomAccessIterator keys, uint64_t num_keys, std::vector<pairs_t>& pairs_blocks, build_configuration const& config, uint32_t bucket_num, uint32_t max_bucket_size) const {
        
        uint64_t divisor;
        pairs_t pairs(num_keys);
        pairs_blocks.push_back(pairs);

        uint64_t min_key=keys[0];
        uint64_t max_key=keys[0];
        for(uint64_t i=1;i!=num_keys;++i){
            if(min_key>keys[i])
                min_key=keys[i];
            if(max_key<keys[i])
                max_key=keys[i];
        }

        divisor=(max_key-min_key+bucket_num-1)/bucket_num;
        // printf("    max_key=%lu,min_key=%lu,num_bucket:%u, divisor=%lu ",max_key,min_key,bucket_num,divisor);
        if(divisor==0)
            divisor=1;

        std::unordered_map<uint32_t, uint32_t> bucket_sizes;

        // generate pairs
        for(uint64_t i=0;i!=num_keys;++i){
            uint32_t bucket_id=(keys[i]-min_key)/divisor;
            // pairs[i] = {static_cast<bucket_id_type>(bucket_id), keys[i]};
            auto hash=hasher_type::hash(keys[i], config.seed);
            pairs_blocks[0][i]={static_cast<bucket_id_type>(bucket_id), hash.first()};
            // search bucket_sizes for bucket_id
            auto it=bucket_sizes.find(bucket_id);
            if(it==bucket_sizes.end())
                bucket_sizes[bucket_id]=1;
            else{
                it->second=it->second+1;
                // if((it->second)>max_bucket_size){
                //     // //printf("    bucket_id=%u, it->second=%u, max_bucket_size=%u\n",bucket_id,it->second,max_bucket_size);
                //     return {0,0};
                // }
            }
            // sloger.func_log(1, "[%u]%lu-%lu\n",bucket_id,keys[i], hash.first());
        }

        std::sort(pairs_blocks[0].begin(), pairs_blocks[0].end());
        // pairs_blocks.push_back(pairs);
        // pairs_blocks.resize(1);
        // pairs_blocks.front().swap(pairs);
        
        return {divisor,min_key};
    }
};

}  // namespace pthash