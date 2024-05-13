#pragma once

#include "include/utils/bucketers.hpp"
#include "include/builders/util.hpp"
#include "include/builders/internal_memory_builder_single_phf.hpp"
#include "include/builders/external_memory_builder_single_phf.hpp"

namespace pthash {

template <typename Hasher, typename Encoder, bool Minimal>
struct single_phf {
    typedef Encoder encoder_type;
    static constexpr bool minimal = Minimal;

    template <typename Iterator>
    build_timings build_in_internal_memory(Iterator keys, uint64_t n,
                                           build_configuration const& config) {
        internal_memory_builder_single_phf<Hasher> builder;
        auto timings = builder.build_from_keys(keys, n, config);
        if(timings.searching_seconds==-1)   return timings;
        timings.encoding_seconds = build(builder, config);
        return timings;
    }

    template <typename Iterator>
    build_timings build_in_external_memory(Iterator keys, uint64_t n,
                                           build_configuration const& config) {
        external_memory_builder_single_phf<Hasher> builder;
        auto timings = builder.build_from_keys(keys, n, config);
        timings.encoding_seconds = build(builder, config);
        return timings;
    }

    template <typename Builder>
    double build(Builder const& builder, build_configuration const& config) {
        // printf("encoding...\n");
        auto start = clock_type::now();
        linear_mapping=config.LinearMapping;
        m_seed = builder.seed();
        m_num_keys = builder.num_keys();
        m_table_size = builder.table_size();
        m_M = fastmod::computeM_u64(m_table_size);
        // printf("    table size: %lu\n", m_table_size);
        m_bucketer = builder.bucketer();
        // printf("    num buckets: %lu\n", m_bucketer.num_buckets());
        // printf("keys_per_bucket: %f\n", double(m_num_keys) / double(m_bucketer.num_buckets()));
        m_pilots.encode(builder.pilots().data(), m_bucketer.num_buckets());
        // printf("    pilots size: %lu B\n", (m_pilots.num_bits())/8);
        if (Minimal and m_num_keys < m_table_size) {
            m_free_slots.encode(builder.free_slots().data(), m_table_size - m_num_keys);
        }
        auto stop = clock_type::now();
        return seconds(stop - start);
    }

    template <typename T>
    uint64_t operator()(T const& key) const {
        if(linear_mapping){
            uint64_t bucket=m_bucketer.linear_bucket(key);
            // printf("key: %lu, bucket: %lu ", key, bucket);
            auto hash=Hasher::hash(key, m_seed);
            // printf("hash: %lu-%lu ", hash.first(), hash.second());
            uint64_t pilot = m_pilots.access(bucket);
            // printf("pilot: %lu ", pilot);
            uint64_t hashed_pilot = default_hash64(pilot, m_seed);
            uint64_t p = fastmod::fastmod_u64(hash.second() ^ hashed_pilot, m_M, m_table_size);
            if constexpr (Minimal) {
                if (PTHASH_LIKELY(p < num_keys())) return p;
                return m_free_slots.access(p - num_keys());
            }
            // printf("p: %lu\n", p);
            return p;
        }

        auto hash = Hasher::hash(key, m_seed);
        return position(hash);
    }

    uint64_t position(typename Hasher::hash_type hash) const {
        uint64_t bucket = m_bucketer.bucket(hash.first());
        uint64_t pilot = m_pilots.access(bucket);
        uint64_t hashed_pilot = default_hash64(pilot, m_seed);
        uint64_t p = fastmod::fastmod_u64(hash.second() ^ hashed_pilot, m_M, m_table_size);
        if constexpr (Minimal) {
            if (PTHASH_LIKELY(p < num_keys())) return p;
            return m_free_slots.access(p - num_keys());
        }
        return p;
    }

    size_t num_bits_for_pilots() const {
        return 8 * (sizeof(m_seed) + sizeof(m_num_keys) + sizeof(m_table_size) + sizeof(m_M)) +
               m_bucketer.num_bits() + m_pilots.num_bits();
    }

    size_t num_bits_for_mapper() const {
        return m_free_slots.num_bits();
    }

    size_t num_bits() const {
        return num_bits_for_pilots() + num_bits_for_mapper();
    }

    inline uint64_t num_keys() const {
        return m_num_keys;
    }

    inline uint64_t table_size() const {
        return m_table_size;
    }

    template <typename Visitor>
    void visit(Visitor& visitor) {
        visitor.visit(m_seed);
        visitor.visit(m_num_keys);
        visitor.visit(m_table_size);
        visitor.visit(m_M);
        visitor.visit(m_bucketer);
        visitor.visit(m_pilots);
        visitor.visit(m_free_slots);
    }

    uint64_t get_slope(){
        if(linear_mapping){
            return m_bucketer.get_divisor();
        }else{
            return m_bucketer.num_buckets();
        }
    }
    uint8_t is_linear_mapping(){
        if(linear_mapping){
            return 1;
        }else{
            return 0;
        }
    }

private:
    uint64_t m_seed;
    uint64_t m_num_keys;
    uint64_t m_table_size;
    __uint128_t m_M;
    skew_bucketer m_bucketer;
    Encoder m_pilots;
    ef_sequence<false> m_free_slots;

    bool linear_mapping;
};

}  // namespace pthash