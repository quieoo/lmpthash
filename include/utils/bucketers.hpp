#pragma once

#include "include/utils/util.hpp"

namespace pthash {

struct skew_bucketer {
    skew_bucketer() {}

    void init(uint64_t num_buckets) {
        // m_num_dense_buckets = constants::b * num_buckets;
        // m_num_sparse_buckets = num_buckets - m_num_dense_buckets;
        // m_M_num_dense_buckets = fastmod::computeM_u64(m_num_dense_buckets);
        // m_M_num_sparse_buckets = fastmod::computeM_u64(m_num_sparse_buckets);
        m_num_buckets=num_buckets;
    }

    inline uint64_t bucket(uint64_t hash) const {
        if(m_num_buckets !=0){
            return hash % m_num_buckets;
        }else{
            printf("error: bucket_num=0\n");
            return 0;
        }

        // static const uint64_t T = constants::a * UINT64_MAX;
        // return (hash < T) ? fastmod::fastmod_u64(hash, m_M_num_dense_buckets, m_num_dense_buckets)
        //                   : m_num_dense_buckets + fastmod::fastmod_u64(hash, m_M_num_sparse_buckets,
        //                                                                m_num_sparse_buckets);
    }

    

    uint64_t num_buckets() const {
        return m_num_buckets;
    }

    size_t num_bits() const {
        return 8 * (sizeof(m_num_dense_buckets) + sizeof(m_num_sparse_buckets) +
                    sizeof(m_M_num_dense_buckets) + sizeof(m_M_num_sparse_buckets));
    }

    void swap(skew_bucketer& other) {
        std::swap(m_num_dense_buckets, other.m_num_dense_buckets);
        std::swap(m_num_sparse_buckets, other.m_num_sparse_buckets);
        std::swap(m_M_num_dense_buckets, other.m_M_num_dense_buckets);
        std::swap(m_M_num_sparse_buckets, other.m_M_num_sparse_buckets);
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_num_dense_buckets);
        visitor.visit(m_num_sparse_buckets);
        visitor.visit(m_M_num_dense_buckets);
        visitor.visit(m_M_num_sparse_buckets);
    }

    void set_divisor(uint64_t divisor) {
        this->divisor = divisor;
    }
    void set_min_key(uint64_t min_key) {
        this->min_key = min_key;
    }
    void set_num_buckets(uint64_t num_buckets) {
        this->m_num_buckets = num_buckets;
    }
    inline uint64_t linear_bucket(uint64_t key)const{
        // printf("key: %lu, min_key: %lu, divisor: %lu\n", key, min_key, divisor);
        return (key-min_key)/divisor;
    }

    uint64_t get_divisor(){
        return divisor;
    }

private:
    uint64_t m_num_dense_buckets, m_num_sparse_buckets;
    __uint128_t m_M_num_dense_buckets, m_M_num_sparse_buckets;
    uint64_t divisor;
    uint64_t min_key;
    uint64_t m_num_buckets=0;
};



struct uniform_bucketer {
    uniform_bucketer() {}

    void init(uint64_t num_buckets) {
        m_num_buckets = num_buckets;
        m_M_num_buckets = fastmod::computeM_u64(m_num_buckets);
    }

    inline uint64_t bucket(uint64_t hash) const {
        return fastmod::fastmod_u64(hash, m_M_num_buckets, m_num_buckets);
    }

    uint64_t num_buckets() const {
        return m_num_buckets;
    }

    size_t num_bits() const {
        return 8 * (sizeof(m_num_buckets) + sizeof(m_M_num_buckets));
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_num_buckets);
        visitor.visit(m_M_num_buckets);
    }

private:
    uint64_t m_num_buckets;
    __uint128_t m_M_num_buckets;
};

}  // namespace pthash