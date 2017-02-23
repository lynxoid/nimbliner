/*
*
* Anchor index: selecting kmers to serve as anchors, read and write an anchor
* collection; anchor locations and is_anchor queries
*
* Author: Darya Filippova
* Dec 2016
*/

#ifndef ANCHOR_INDEX
#define ANCHOR_INDEX

#include <unordered_map>
#include <fstream>
#include <chrono>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <list>

#include "JellyfishUtil.h"
#include "definitions.hpp"

using namespace std;

namespace nimble {

class AnchorIndex {

    shared_ptr<unordered_map<kmer_t, vector<genomic_coordinate_t>>> _anchors;

public:

    AnchorIndex() {
    }

    static constexpr const char* EXT = ".star";

    /*
     * returns True is the kmer is in the anchor set; false otherwise
     */
    bool is_anchor(const bin_kmer_t kmer) const {
		return _anchors->find(kmer) != _anchors->end();
    }

    /*
     * returns a list of genomic_coordinate_t locations corresponding to all
     * occurrences of this anchor. If kmer is not an anchor, returns an emplty
     * vector
     */
    vector<genomic_coordinate_t> & get_anchor_locations(const bin_kmer_t & kmer) const {
        return (*_anchors)[kmer];
    }

    uint size() {
        return _anchors->size();
    }

    /*
	 *
	 * read kmers that serve as anchors
	 *
	 */
    // TODO: rename to readAnchors
	void readStarLocations(const string & prefix, const int K) {
		cerr << "[AnchorIndex] reading anchors from " << prefix << EXT << endl;
		auto start = std::chrono::system_clock::now();

		_anchors = shared_ptr<unordered_map<kmer_t, vector<genomic_coordinate_t>>>(
                new unordered_map<kmer_t, vector<genomic_coordinate_t>>()
            );
		ifstream in(prefix + EXT);

        if (!in) {
            cerr << "[ERROR] [AnchorIndex] Could not open anchor index file " << prefix << EXT << endl;
            exit(1);
        }

		string line;
		while (getline(in, line)) {
			// parse the line now
			std::stringstream ss(line);
			string kmer, pos;
			getline(ss, kmer, ' ');
			uint64_t bin_kmer = stol(kmer);
			_anchors->emplace(bin_kmer, vector<genomic_coordinate_t>{});
			while(getline(ss, pos, ' ')) {
	    		// int location = stoul(pos);
	    		int location = stoi(pos);
	    		(*_anchors)[bin_kmer].push_back(location);
			}
            // TODO: this should be done during anchor selection phase, not here
            // remove anchors that occur more than 10K times
			if ((*_anchors)[bin_kmer].size() > 10000) {
				cerr << nimble::mer_binary_to_string(bin_kmer, K) << " " <<
					(*_anchors)[bin_kmer].size() << endl;
				_anchors->erase( _anchors->find(bin_kmer) );
			}
		}
		in.close();

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		cerr << "[AnchorIndex] reading anchors took: " << elapsed_seconds.count() << "s" << endl;
		// return kmer_locations;
	}

    /*
	 * TODO: sort and delta encode offsets; write out to a gzip
	 */
	static void write_anchors(shared_ptr<unordered_map<kmer_t, list<seed_position_t>>> anchors,
        const string & output_prefix) {
		cerr << "[AnchorIndex] writing anchors" << endl;

		ofstream star_locations_out(output_prefix + EXT);
        // TODO: check if can write
		auto start = std::chrono::system_clock::now();

        reference_id_t prev_ref_id = -1;
		for (auto & anchor_pair : *anchors) {
			star_locations_out << anchor_pair.first << " ";
			for (const seed_position_t & loc : anchor_pair.second) {
                if (prev_ref_id != loc.first) {
                    star_locations_out << "*" << loc.first << " ";
                    prev_ref_id = loc.first;
                }
                star_locations_out << loc.second << " ";
            }
			star_locations_out << endl;
		}
		star_locations_out.close();

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds_str = end - start;
	    cerr << "[AnchorIndex] Saving anchors took: " << elapsed_seconds_str.count() << "s" << endl;
	}
};

}

#endif // ANCHOR_INDEX
