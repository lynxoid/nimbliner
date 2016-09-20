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

// #include "BaseBloomFilter.hpp"
// #include "bit_tree_binary.hpp"
#include "JellyfishUtil.h"
#include "definitions.hpp"

using namespace std;

class ReferenceIndex {

public:

	ReferenceIndex() {};

	// read kmers (or pdBG describing the reference) and anchor locations
	virtual void readIndex(const string & kmers_path, const string & stars_path, const uint K) =0;

	/* returns true if this kmer was present in the reference sequence, false otherwise */
	virtual bool has_kmer(const bin_kmer_t & kmer) const =0;

	/* returns true is this kmer is found among anchors, false otherwise*/
	virtual bool has_anchor(const bin_kmer_t & kmer) const =0;

	// TODO: what does this & do? do we use it?
	virtual vector<genomic_coordinate_t> & get_anchor_locations(const bin_kmer_t & kmer) const =0;

	/*
	 *
	 * read kmers that serve as anchors
	 *
	 */
	shared_ptr<unordered_map<kmer_t, vector<genomic_coordinate_t>>> readStarLocations(const string & path, const int K) {
		cerr << "reading stars from " << path << endl;
		auto start = std::chrono::system_clock::now();

		shared_ptr<unordered_map<kmer_t, vector<genomic_coordinate_t>>> kmer_locations = 
			shared_ptr<unordered_map<kmer_t, vector<genomic_coordinate_t>>>(new unordered_map<kmer_t, vector<genomic_coordinate_t>>());
		ifstream in(path);
		string line;
		while (getline(in, line)) {
			// parse the line now
			std::stringstream ss(line);
			string kmer, pos;
			getline(ss, kmer, ' ');
			uint64_t bin_kmer = stol(kmer);
			kmer_locations->emplace(bin_kmer, vector<genomic_coordinate_t>{});
			while(getline(ss, pos, ' ')) {
	    		// int location = stoul(pos);
	    		int location = stoi(pos);
	    		(*kmer_locations)[bin_kmer].push_back(location);
			}
			if ((*kmer_locations)[bin_kmer].size() > 10000) {
				cerr << nimble::mer_binary_to_string(bin_kmer, K) << " " << 
					(*kmer_locations)[bin_kmer].size() << endl;
				kmer_locations->erase( kmer_locations->find(bin_kmer) );
			}
		}
		in.close();

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		cerr << "reading anchors took: " << elapsed_seconds.count() << "s" << endl;

		return kmer_locations;
	}
};

#endif // REFERENCE_INDEX_READER