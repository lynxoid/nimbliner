//
// IndexReader
//

#ifndef REFERENCE_INDEX_READER
#define REFERENCE_INDEX_READER

#include <vector>
#include <string>
#include <unordered_map>
#include <fstream>
#include <sstream>
#include <chrono>

// #include "BaseBloomFilter.hpp"
// #include "bit_tree_binary.hpp"
#include "JellyfishUtil.h"
#include "definitions.hpp"

using namespace std;

class ReferenceIndex {

public:

	ReferenceIndex() {};

	virtual int getK() =0;

	// read kmers (or pdBG describing the reference) and anchor locations
	virtual void readIndex(const string & kmers_path, const string & stars_path) =0;

	/* returns true if this kmer was present in the reference sequence, false otherwise */
	virtual bool has_kmer(const bin_kmer_t kmer) const =0;

	/* returns true is this kmer is found among anchors, false otherwise*/
	virtual bool is_anchor(const bin_kmer_t kmer) const =0;

	// TODO: what does this & do? do we use it?
	virtual vector<genomic_coordinate_t> & get_anchor_locations(const bin_kmer_t & kmer) const =0;

};

#endif // REFERENCE_INDEX_READER
