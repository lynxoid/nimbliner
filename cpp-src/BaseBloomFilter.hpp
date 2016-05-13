#ifndef LIB_BASE_BLOOM_FILTER
#define LIB_BASE_BLOOM_FILTER

#include <string>

#include <bf.h>

using namespace std;

// kmer type
typedef uint64_t kmer_t;

class BaseBloomFilter {
protected:
	bf::basic_bloom_filter bf_;	// generic bloom filter for holding the entries
	const int k; // kmer length

public:

	BaseBloomFilter(const int k, const size_t num_elems=1024*1024*32)
	: bf_(bf::make_hasher(1), num_elems), k(k)
	{};

	void add(const kmer_t & kmer) {
		bf_.add(kmer);
	}

	virtual bool contains(const kmer_t & kmer) {
		return bf_.lookup(kmer);
	}

	virtual void populate(const unordered_set<kmer_t> & kmer_set){
		for(auto km : kmer_set){
			add(km);
		}
	}

};

#endif
