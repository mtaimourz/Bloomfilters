
#ifndef BLOOM_FILTER_HPP
#define BLOOM_FILTER_HPP

#include <cstdint>
#include <vector>
#include <string>
#include <functional>
#include <cmath>
#include <fstream>

class BloomFilter {
public:
    // alias for hash function; input: string and int; returns vector | shorthand for a long type
    using HashFunction = std::function<std::vector<uint64_t>(const std::string&, unsigned int)>; 
    
    BloomFilter(uint64_t estimated_elements, double false_positive_rate, HashFunction hash_func = nullptr);
    
    explicit BloomFilter(const std::string& filepath);
    
    void add(const std::string& key);
    bool contains(const std::string& key) const;
    
    bool save(const std::string& filepath) const;
    bool load(const std::string& filepath);
    
    uint64_t elements_added() const { return m_elements_added; }
    double current_false_positive_rate() const;
    void print_stats() const;
    
    static BloomFilter union_filters(const BloomFilter& bf1, const BloomFilter& bf2);
    static BloomFilter intersect_filters(const BloomFilter& bf1, const BloomFilter& bf2);
    double jaccard_index(const BloomFilter& other) const;
    
private:
    std::vector<uint8_t> m_bloom;
    uint64_t m_estimated_elements;
    double m_false_positive_rate;
    unsigned int m_number_hashes;
    uint64_t m_number_bits;
    uint64_t m_elements_added;
    
    HashFunction m_hash_function;
    
    void calculate_optimal_parameters();
    void set_bit(uint64_t index);
    bool get_bit(uint64_t index) const;
    // fallback hash func if user does not provide one | uses FNV-1a
    std::vector<uint64_t> default_hash(const std::string& key) const; 
    
    // Default hash functions (FNV-1a based)
    static uint64_t fnv1a_hash(const std::string& str);
};

#endif // BLOOM_FILTER_HPP