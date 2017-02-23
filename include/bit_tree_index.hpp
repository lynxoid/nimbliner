#ifndef BIT_TREE_INDEX_LIB
#define BIT_TREE_INDEX_LIB

#include <string>
#include <vector>
#include <memory>
#include <boost/progress.hpp>

#include "reference_index.hpp"
#include "anchor_index.hpp"
#include "bit_tree_binary.hpp"

class BitTreeIndex : public ReferenceIndex {

	int K;

	shared_ptr<BitTreeBin> _bit_tree;

	// shared_ptr<unordered_map<kmer_t, vector<genomic_coordinate_t>>> _stars;
    nimble::AnchorIndex _anchorIndex;

	shared_ptr<BitTreeBin> readBitTree(const string & kmers_path, const uint K) {
		shared_ptr<BitTreeBin> bit_tree = shared_ptr<BitTreeBin>(new BitTreeBin());
		int k;
		bit_tree->read(kmers_path, k);
		assert(k == K);
		assert(K > 0 && K <= 32);
		return bit_tree;
	}

public:

	uint getK() {return K;}

	// read kmers (or pdBG describing the reference) and anchor locations
	void readIndex(const string & kmers_path, const string & stars_path) {
		boost::timer t;
		_bit_tree = readBitTree(kmers_path, K);
		cerr << "Reading took: " << t.elapsed() << "s" << endl;
		t.restart();
		_bit_tree->build_index(10);
		cerr << "Building btTree index took: " << t.elapsed() << "s" << endl;
		_anchorIndex.readStarLocations(stars_path, K);
	}

	/* returns true if this kmer was present in the reference sequence, false otherwise */
	bool has_kmer(const bin_kmer_t kmer) const {
		return _bit_tree->contains(kmer);
	}

	/* returns true is this kmer is found among anchors, false otherwise*/
	bool is_anchor(const bin_kmer_t kmer) const {
		// return _stars->find(kmer) != _stars->end();
        return _anchorIndex.is_anchor(kmer);
	}

	// TODO: what does this & do? do we use it?
	vector<genomic_coordinate_t> & get_anchor_locations(const bin_kmer_t & kmer) const {
		// return (*_stars)[kmer];
        return _anchorIndex.get_anchor_locations(kmer);
	}

	shared_ptr<vector<kmer_type>> get_kmers() {
		 auto kmers = _bit_tree->decode();
		 return kmers;
	}
};

#endif // BIT_TREE_INDEX_LIB
