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

    shared_ptr<unordered_map<kmer_t, vector<seed_position_t>>> _anchors;

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
    vector<seed_position_t> & get_anchor_locations(const bin_kmer_t & kmer) const {
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
	void readAnchors(const string & prefix, const int K) {
		cerr << "[AnchorIndex] reading anchors from " << prefix << EXT << endl;
		auto start = std::chrono::system_clock::now();

		_anchors = shared_ptr<unordered_map<kmer_t, vector<seed_position_t>>>(
                new unordered_map<kmer_t, vector<seed_position_t> >() );

		ifstream in(prefix + EXT);
        can_read_or_quit(in, prefix+EXT, true);
        // if (!in) {
            // cerr << "[ERROR] [AnchorIndex] Could not open anchor index file " << prefix << EXT << endl;
            // exit(1);
        // }

		string line;
		while (getline(in, line)) {
			// parse the line now
			std::stringstream ss(line);
			string kmer, pos;
			getline(ss, kmer, ' ');
			kmer_t bin_kmer = stol(kmer);
			_anchors->emplace(bin_kmer, vector<seed_position_t>() );

            reference_id_t ref_id = -1;
			while(getline(ss, pos, ' ')) {
	    		// int location = stoul(pos);
                if (pos[0] == '*') {
                    ref_id = stoi(pos.substr(1));
                }
                else {
	    		    genomic_coordinate_t location = stoi(pos);
                    (*_anchors)[bin_kmer].emplace_back(ref_id, location);
                }
			}
		}
		in.close();

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		cerr << "[AnchorIndex] reading anchors took: " << elapsed_seconds.count() << "s" << endl;
	}

    /*
	 * TODO: sort and delta encode offsets; write out to a gzip
   * Even before gzipping --- at least writing this out in binary would be much faster.
	 */
	static void write_anchors(shared_ptr<unordered_map<kmer_t, list<seed_position_t>>> anchors,
        const string & output_prefix) {
		cerr << "[AnchorIndex] writing anchors" << endl;

		ofstream star_locations_out(output_prefix + EXT);
        // TODO: check if can write
		auto start = std::chrono::system_clock::now();

    size_t i{0};
		for (auto & anchor_pair : *anchors) {
      //if (i % 10000 == 0) { std::cerr << "writing anchor pair " << i << "\n"; }
            reference_id_t prev_ref_id = -1;
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
