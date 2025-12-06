#include "bloom_filter.hpp"
#include <cstring>
#include <stdexcept>
#include <iostream>
#include <vector>

using namespace std;

/*
  * -----------------------------
  * 
  * 
  * 
  * Constructors
  * 
  * 
  * 
  * -----------------------------
*/

// by default, hash_func is nullptr
BloomFilter::BloomFilter(uint64_t estimated_elements, double false_positive_rate, HashFunction hash_func)
    : m_estimated_elements(estimated_elements),
      m_false_positive_rate(false_positive_rate),
      m_elements_added(0),
      m_hash_function(hash_func)
{
    if (estimated_elements == 0) {
        throw std::invalid_argument("Estimated elements must be > 0");
    }
    
    if (false_positive_rate <= 0.0 || false_positive_rate >= 1.0) {
        throw std::invalid_argument("False positive rate must be between 0 and 1");
    }
    
    calculate_optimal_parameters();
    
    size_t num_bytes = (m_number_bits + 7) / 8; // m_number_bits are set in calculate_optimal_parameters()
    m_bloom.resize(num_bytes, 0);
}


// Constructor when we give the file
BloomFilter::BloomFilter(const std::string& filepath)
    : m_elements_added(0) {
    if (!load(filepath)) {
        throw std::runtime_error("Failed to load bloom filter from file");
    }
}


/*
  * -----------------------------
  * 
  * 
  * 
  * Helper Method -> calculate_optimal_parameters()
  *     sets: 1) m_number_bits   : determines Space usage
  *           2) m_number_hashes : determines False positive rate
  * 
  *     why:  1) Too few bits → many collisions → very high false positive rate
  *           2) Too many bits → wastes memory    
  * 
  * 
  * 
  * -----------------------------
*/

void BloomFilter::calculate_optimal_parameters() {
    // Formula: m = -n * ln(p) / (ln(2)^2)
    // Where: m = number of bits, n = estimated elements, p = false positive rate
    double ln2_squared = std::log(2.0) * std::log(2.0);
    m_number_bits = static_cast<uint64_t>(
        std::ceil(-1.0 * m_estimated_elements * std::log(m_false_positive_rate) / ln2_squared)
    );
    
    // Formula: k = (m/n) * ln(2)
    // Where: k = number of hash functions
    m_number_hashes = static_cast<unsigned int>(
        std::ceil((m_number_bits / static_cast<double>(m_estimated_elements)) * std::log(2.0))
    );
    
    // make sure we have at least 1 hash function
    if (m_number_hashes < 1) {
        m_number_hashes = 1;
    }
}

/*
  * -----------------------------
  * 
  * 
  * 
  *  Bit Manipulation
  *     
  *     set_bit -> Sets a specific bit to 1 using bitwise OR
  *     get_bit -> Checks whether a specific bit is 1 using bitwise AND
  * 
  * 
  *  eg; set bit at INDEX 30
  *                   ├─ byte_index = 30 ÷ 8 = 3
  *                   └─ bit_offset = 30 % 8 = 6
  * 
  *                     Byte 0: [bit 0]  [bit 1]   [bit 2]  [bit 3]   [bit 4]  [bit 5]   [bit 6]   [bit 7]
  *                     Byte 1: [bit 8]  [bit 9]   [bit 10] [bit 11]  [bit 12] [bit 13]  [bit 14]  [bit 15]
  *                     Byte 2: [bit 16] [bit 17]  [bit 18] [bit 19]  [bit 20] [bit 21]  [bit 22]  [bit 23]
  *                     Byte 3: [bit 24] [bit 25]  [bit 26] [bit 27]  [bit 28] [bit 29]  [bit 30]  [bit 31]
  *                                                                                          ┬
  *                                                                                      bit 30 is HERE!
  *                                                                                      Position 6 in byte 3
  *                    Step 1: Create mask
  *                            1 << 6 = 01000000
  *                    
  *                    Step 2: Apply mask to byte 3
  *                            m_bloom[3] = 00000000  (before)
  *                            mask       = 01000000
  *                                      OR ────────
  *                            m_bloom[3] = 01000000  (after) 
  * 
  *                     BEFORE: 00000000 00000000 00000000 00000000
  *                     AFTER:  01000000 00000000 00000000 00000000
  *
  *  eg; get bit at INDEX 30
  *     
  *                      Step 1: Create same mask
  *                              1 << 6 = 01000000
  *                      
  *                      Step 2: Check byte 3
  *                              m_bloom[3] = 01000000  (current value after setting)
  *                              mask       = 01000000
  *                                       AND ────────
  *                              result     = 01000000
  * 
  * 
  * -----------------------------
*/

void BloomFilter::set_bit(uint64_t index) {
    if (index >= m_number_bits) {
        throw std::out_of_range("Bit index out of range");
    }
    
    uint64_t byte_index = index / 8;
    uint8_t bit_offset = index % 8;
    
    m_bloom[byte_index] |= (1 << bit_offset);
}

bool BloomFilter::get_bit(uint64_t index) const {
    if (index >= m_number_bits) {
        throw std::out_of_range("Bit index out of range");
    }
    
    uint64_t byte_index = index / 8;
    uint8_t bit_offset = index % 8;
    
    return (m_bloom[byte_index] & (1 << bit_offset)) != 0;
}


/*
  * -----------------------------
  * 
  * 
  * 
  *  Hash Function Implementation
  *     
  *     generates different index in Bloom filter’s bit array (m_bloom, gets set in the constructor )
  * 
  * 
  * 
  * -----------------------------
*/

uint64_t BloomFilter::fnv1a_hash(const std::string& str) {
    // FNV-1a 64-bit constants
    const uint64_t FNV_OFFSET_BASIS = 14695981039346656037ULL;
    const uint64_t FNV_PRIME = 1099511628211ULL;
    
    uint64_t hash = FNV_OFFSET_BASIS;
    
    for (char c : str) {
        hash ^= static_cast<uint8_t>(c);
        hash *= FNV_PRIME;
    }
    
    return hash;
}

std::vector<uint64_t> BloomFilter::default_hash(const std::string& key) const {
    // eg: 
    //    key:  "https://example.com/page3"
    //    salt: "salt_https://example.com/page3"

    std::vector<uint64_t> hashes;
    hashes.reserve(m_number_hashes);
    
    uint64_t h1 = fnv1a_hash(key);
    
    
    std::string salted_key = "salt_" + key;
    uint64_t h2 = fnv1a_hash(salted_key);
    
    for (unsigned int i = 0; i < m_number_hashes; ++i) {
        uint64_t combined_hash = (h1 + i * h2) % m_number_bits; // m_number_bits: makes sure we are withing bounds
        hashes.push_back(combined_hash);
    }
    
    return hashes;
}


/*
  * -----------------------------
  * 
  * 
  * 
  *  Add() and Contains()
  * 
  * 
  * 
  * -----------------------------
*/


void BloomFilter::add(const std::string& key) {
    std::vector<uint64_t> hashes;
    
    // check if the hash_func exists or not
    if (m_hash_function) {
        hashes = m_hash_function(key, m_number_hashes);
    } else {
        hashes = default_hash(key);
    }
    
    for (uint64_t hash : hashes) {
        set_bit(hash);
    }
    
    m_elements_added++;
}

bool BloomFilter::contains(const std::string& key) const {
    
    std::vector<uint64_t> hashes;
    
    if (m_hash_function) {
        hashes = m_hash_function(key, m_number_hashes);
    } else {
        hashes = default_hash(key);
    }
    
    for (uint64_t hash : hashes) {
        if (!get_bit(hash)) {
            return false;  // Definitely not in the set
        }
    }
    
    return true;  // Probably in the set
}


/*
  * -----------------------------
  * 
  *  Serialization
  *             save() -> converts our in-memory data (Bloom filter object to sequence of bytes and write them on file
  * 
  *  Deserialization
  *             load() -> reads the bytes from our file and reconstructs the Bloom filter object in memory.
  *
  *  Stored bytes mapping
  *              Bytes 0-7   : estimated_elements (uint64_t)
  *              Bytes 8-15  : false_positive_rate (double)
  *              Bytes 16-19 : number_hashes (unsigned int)
  *              Bytes 20-27 : number_bits (uint64_t)
  *              Bytes 28-35 : elements_added (uint64_t)
  *              Bytes 36-43 : bloom_size (size_t)
  *              Bytes 44+   : bit array data 
  * 
  * -----------------------------
*/


bool BloomFilter::save(const std::string& filepath) const {
    std::ofstream file(filepath, std::ios::binary);
    if (!file) {
        return false;
    }
    
    file.write(reinterpret_cast<const char*>(&m_estimated_elements), sizeof(m_estimated_elements));
    file.write(reinterpret_cast<const char*>(&m_false_positive_rate), sizeof(m_false_positive_rate));
    file.write(reinterpret_cast<const char*>(&m_number_hashes), sizeof(m_number_hashes));
    file.write(reinterpret_cast<const char*>(&m_number_bits), sizeof(m_number_bits));
    file.write(reinterpret_cast<const char*>(&m_elements_added), sizeof(m_elements_added));
    
    size_t bloom_size = m_bloom.size();
    file.write(reinterpret_cast<const char*>(&bloom_size), sizeof(bloom_size));
    file.write(reinterpret_cast<const char*>(m_bloom.data()), bloom_size);
    
    return file.good();
}

bool BloomFilter::load(const std::string& filepath) {
    std::ifstream file(filepath, std::ios::binary);
    if (!file) {
        return false;
    }
    
    file.read(reinterpret_cast<char*>(&m_estimated_elements), sizeof(m_estimated_elements));
    file.read(reinterpret_cast<char*>(&m_false_positive_rate), sizeof(m_false_positive_rate));
    file.read(reinterpret_cast<char*>(&m_number_hashes), sizeof(m_number_hashes));
    file.read(reinterpret_cast<char*>(&m_number_bits), sizeof(m_number_bits));
    file.read(reinterpret_cast<char*>(&m_elements_added), sizeof(m_elements_added));
    
    size_t bloom_size;
    file.read(reinterpret_cast<char*>(&bloom_size), sizeof(bloom_size));
    m_bloom.resize(bloom_size);
    file.read(reinterpret_cast<char*>(m_bloom.data()), bloom_size);
    
    return file.good();
}

/*
  * -----------------------------
  * 
  * 
  * 
  *  Stats
  * 
  * 
  * 
  * -----------------------------
*/



double BloomFilter::current_false_positive_rate() const {
    if (m_elements_added == 0) {
        return 0.0;
    }
    
    // Formula: p = (1 - e^(-k*n/m))^k
    double exponent = -1.0 * m_number_hashes * m_elements_added / 
                      static_cast<double>(m_number_bits);
    return std::pow(1.0 - std::exp(exponent), m_number_hashes);
}

void BloomFilter::print_stats() const {
    uint64_t bits_set = 0;
    for (uint64_t i = 0; i < m_number_bits; ++i) {
        if (get_bit(i)) {
            bits_set++;
        }
    }
    std::cout << "Estimated elements:           " << m_estimated_elements << std::endl;
    std::cout << "Target false positive rate:   " << m_false_positive_rate << std::endl;
    std::cout << "Number of bits:               " << m_number_bits << std::endl;
    std::cout << "Number of bytes:              " << m_bloom.size() << std::endl;
    std::cout << "Number of hash functions:     " << m_number_hashes << std::endl;
    std::cout << "Elements added:               " << m_elements_added << std::endl;
    std::cout << "Current false positive rate:  " << (current_false_positive_rate() * 100) << "%" << std::endl;
    std::cout << "Bits set: " << bits_set << " / " << m_number_bits << " (" << (bits_set * 100.0 / m_number_bits) << "%)" << std::endl;
}

/*
  * -----------------------------
  * 
  * 
  * 
  *  Set Operrations
  * 
  *               1) union_filters  -> Creates a Bloom filter that represents the union of two sets
  *                             
  *                             eg: 
  *                                 bf1 bits: 1010 1100
  *                                 bf2 bits: 1100 1010
  *                                 union:    1110 1110
  * 
  *               2) intersect_filters -> Creates a Bloom filter that represents the intersection of two sets
  * 
  *                            eg;
  *                                 bf1 bits:      1010 1100
  *                                 bf2 bits:      1100 1010
  *                                 intersection:  1000 1000
  * 
  *               3) jaccard_index  -> Measures similarity between two Bloom filters (0.0 to 1.0).
  *                             
  *                             eg: 
  *                                 bf1 bits: 1010
  *                                 bf2 bits: 1100
  *                                 union: 1110 (3 bits)
  *                                 intersection: 1000 (1 bit)
  *                                 Jaccard = 1 / 3 ≈ 0.333  
  * 
  * -----------------------------
*/


BloomFilter BloomFilter::union_filters(const BloomFilter& bf1, const BloomFilter& bf2) {
    // must be compatible
    if (bf1.m_number_bits != bf2.m_number_bits || 
        bf1.m_number_hashes != bf2.m_number_hashes) {
        throw std::invalid_argument("Bloom filters must have same parameters");
    }
    
    BloomFilter result(bf1.m_estimated_elements, bf1.m_false_positive_rate);
    
    // OR all bits
    for (size_t i = 0; i < bf1.m_bloom.size(); ++i) {
        result.m_bloom[i] = bf1.m_bloom[i] | bf2.m_bloom[i];
    }
    
    result.m_elements_added = bf1.m_elements_added + bf2.m_elements_added;
    
    return result;
}

BloomFilter BloomFilter::intersect_filters(const BloomFilter& bf1, const BloomFilter& bf2) {
    if (bf1.m_number_bits != bf2.m_number_bits || 
        bf1.m_number_hashes != bf2.m_number_hashes) {
        throw std::invalid_argument("Bloom filters must have same parameters");
    }
    
    BloomFilter result(bf1.m_estimated_elements, bf1.m_false_positive_rate);
    
    // AND all bits
    for (size_t i = 0; i < bf1.m_bloom.size(); ++i) {
        result.m_bloom[i] = bf1.m_bloom[i] & bf2.m_bloom[i];
    }
    
    result.m_elements_added = std::min(bf1.m_elements_added, bf2.m_elements_added);
    
    return result;
}

double BloomFilter::jaccard_index(const BloomFilter& other) const {
    if (m_number_bits != other.m_number_bits) {
        throw std::invalid_argument("Bloom filters must have same size");
    }
    
    uint64_t bits_union = 0;
    uint64_t bits_intersection = 0;
    
    for (size_t i = 0; i < m_bloom.size(); ++i) {
        uint8_t union_byte = m_bloom[i] | other.m_bloom[i];
        uint8_t intersection_byte = m_bloom[i] & other.m_bloom[i];
        
        // Count set bits
        for (int j = 0; j < 8; ++j) {
            if (union_byte & (1 << j)) bits_union++;
            if (intersection_byte & (1 << j)) bits_intersection++;
        }
    }
    
    if (bits_union == 0) {
        return 0.0;
    }
    
    return static_cast<double>(bits_intersection) / bits_union;
}



/*
  * -----------------------------
  * 
  * 
  * 
  *  Main
  * 
  * 
  * 
  * -----------------------------
*/

int main() {
    try {
        BloomFilter url_filter(10000, 0.01); // 10,000 URLs, 1% false positive rate
        
        std::vector<std::string> urls = {
            "https://example.com/page1",
            "https://example.com/page2",
            "https://example.com/page3"
        };
        
        for (const auto& url : urls) {
            url_filter.add(url);
        }
        
        std::cout << "\nChecking membership:" << std::endl;
        std::cout << "Contains page1? " << (url_filter.contains("https://example.com/page1") ? "Yes" : "No") << std::endl;
        std::cout << "Contains page99? " << (url_filter.contains("https://example.com/page99") ? "Yes" : "No") << std::endl;
        
        std::cout << "\nSaving to disk..." << std::endl;
        if (url_filter.save("url_filter.bloom")) std::cout << "Saved successfully!" << std::endl;
        
        std::cout << "\nLoading from disk..." << std::endl;
        BloomFilter loaded_filter("url_filter.bloom");
        std::cout << "Loaded successfully!" << std::endl;
        std::cout << "Loaded filter contains page1? " 
                  << (loaded_filter.contains("https://example.com/page1") ? "Yes" : "No") << std::endl;
        
        // Print stats
        std::cout << std::endl;
        url_filter.print_stats();
        
        // Test false positive rate
        std::cout << "\nTesting false positive rate..." << std::endl;
        int false_positives = 0;
        int test_count = 10000;
        
        for (int i = 10000; i < 10000 + test_count; ++i) {
            std::string test_url = "https://example.com/page" + std::to_string(i);
            if (url_filter.contains(test_url)) {
                false_positives++;
            }
        }
        
        double actual_fpr = static_cast<double>(false_positives) / test_count;
        std::cout << "Actual false positive rate: " << (actual_fpr * 100) << "%" << std::endl;
        std::cout << "Theoretical false positive rate: " 
                  << (url_filter.current_false_positive_rate() * 100) << "%" << std::endl;
        
        // Set operations
        std::cout << "\n Set Operations" << std::endl;
        BloomFilter filter1(1000, 0.01);
        BloomFilter filter2(1000, 0.01);
        
        filter1.add("item1");
        filter1.add("item2");
        filter1.add("common");
        
        filter2.add("item3");
        filter2.add("item4");
        filter2.add("common");
        
        BloomFilter union_filter = BloomFilter::union_filters(filter1, filter2);
        std::cout << "Union contains 'item1': " << (union_filter.contains("item1") ? "Yes" : "No") << std::endl;
        std::cout << "Union contains 'item3': " << (union_filter.contains("item3") ? "Yes" : "No") << std::endl;
        
        BloomFilter intersect_filter = BloomFilter::intersect_filters(filter1, filter2);
        std::cout << "Intersection contains 'item1': " << (intersect_filter.contains("item1") ? "Yes" : "No") << std::endl;
        std::cout << "Intersection contains 'item3': " << (intersect_filter.contains("item3") ? "Yes" : "No") << std::endl;
        std::cout << "Intersection contains 'common': " << (intersect_filter.contains("common") ? "Yes" : "No") << std::endl;
        
        std::cout << "Jaccard index: " << filter1.jaccard_index(filter2) << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}