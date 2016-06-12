//
// Aligner -- the main class, i guess...
//
// dasha.filippova@gmail.com
//
// 2016

#ifndef ALIGNER
#define ALIGNER

#include "reference_index.hpp"
#include "definitions.hpp"

class Aligner {

	ReferenceIndex _index;

	vector<int> findAllMatchingAnchorPositions(const vector<pair<kmer_t,int>> & matched_stars, 
		const ReferenceIndex & index) {
		assert(matched_stars.size() > 0);

		vector<int> mappings;
		auto first_star = matched_stars[0];
		auto second_star = matched_stars[1];
		int delta = second_star.second - first_star.second;

		auto A = index.get_anchor_locations(first_star.first);
		auto B = index.get_anchor_locations(second_star.first);
		// auto A = star_locations[first_star.first];
		// auto B = star_locations[second_star.first];
		int i = 0, j = 0;
		while (i < A.size() && j < B.size() ) {
			if (A[i] < B[j]) {
				if (B[j] - A[i] == delta) {
					// found one match
					mappings.push_back(A[i] - first_star.second);
					i++;
					j++;
				} 
				else if (B[j] - A[i] < delta) {
					j++;
				}
				else {
					i++;
				}
			}
			else {
				j++;
			}
		}
		return mappings;
	}

	////////////////////////////////////////////////////////
	//	matched_stars -- in order in which they appear in the read
	//	resolve which of the star co-locations are the same distance
	//	apart as the stars in the read
	//  0123456789
	//	--***-***-
	//
	////////////////////////////////////////////////////////
	vector<int> resolve_mapping_locations(const vector<pair<kmer_t,int>> & matched_stars, 
		// unordered_map<kmer_t, vector<int>> & star_locations,
		int & extend) {
		vector<int> mappings;

		// what if no stars mapped?
		if (matched_stars.size() == 0) {
			// need to extend the read until we hit some star
			extend++;
		}
		if (matched_stars.size() == 1) {
			// that's the only thing we got going
			// TODO: extend to find more stars
			auto reference_locations = _index.get_anchor_locations(matched_stars[0].first);
			for (auto loc : reference_locations)
				mappings.push_back(loc - matched_stars[0].second);
		}
		else if (matched_stars.size() >= 2) {
			// if two stars	or more
			// TODO: take into account all the stars

			mappings = findAllMatchingAnchorPositions(matched_stars, _index);
		}

		return mappings;
	}

	void get_next_kmer(kmer_t & bin_kmer, const char next_base, const int K) {
		kmer_t mask = ( ((kmer_t)1) << ( 2*(K-1) ) ) - 1;
        bin_kmer = (bin_kmer & mask) << 2;
        // append a new char on the right
        bin_kmer = bin_kmer | dna_codes[ (size_t) seq->seq.s[i] ];
	}

public:

	Aligner(const ReferenceIndex & index) {
		_index = index;
	}

	////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////
	void alignReads(const string & path, int K) {
		FastaReader fr(path.c_str());
	    kseq_t * seq;
	    int passed_cutoff = 0;
	    int read_count = 0;
	    int need_to_extend_read = 0;

	    std::chrono::duration<double> elapsed_seconds_str;
	    while ( (seq = fr.nextSequence() ) ) {
	    	read_count++;
	    	if (read_count % 100000 == 0) cerr << read_count / 100000 << "00K ";
	        if (seq->seq.l < K) {
	        	continue;
	        }
	        vector<pair<kmer_t,int>> matched_stars;
	        int cnt_matched = 0;
	        int L = seq->seq.l - K + 1;
	        // only allow reads shorter than 2^16
	        assert(L < 65536);
	        auto bin_kmer = mer_string_to_binary(seq->seq.s, K);
	        //if (bf.contains(bin_kmer) ) {
	        //        cnt_matched++;
	                // check if hit a star kmer
	                if ( _index.has_anchor(bin_kmer) )
	                        matched_stars.emplace_back(bin_kmer, 0);
	        //}
	        // now go through the rest of the read
	        unsigned short i = 1;
	        while (i < L) {
                // update prev kmer // mask the leftmost character
        		get_next_kmer(bin_kmer, seq->seq.s[i], K);
                
                if (bf.contains(bin_kmer) ) {
                    cnt_matched++;
                    // check if hit a star kmer
                    if ( _index.has_anchor(bin_kmer) )
                            matched_stars.emplace_back(bin_kmer, i - K + 1);
                    // if (matched_stars.size() >= 2) break;
                }
                else {
                	// mismatch? indel?
                	// was prev kmer there?
                	// is next kmer there too?
                }
                i++;
	        }
	        // require at least 50% of all kmers to match
	        // if (cnt_matched / L < 0.45) continue;
	        // passed_cutoff++;
			cerr << seq->name.s << "\tmatched stars: " << matched_stars.size() << endl;
			auto start = std::chrono::system_clock::now();

			// resolve star kmers to get an exact mapping location
			auto mapping_locations = resolve_mapping_locations(matched_stars, need_to_extend_read);

			// output all potential locations for this read
			cout << seq->name.s << "\t";
			for (const auto & loc : mapping_locations) {
				cout << loc << " ";
			}
			cout << endl;

			auto end = std::chrono::system_clock::now();
			elapsed_seconds_str += end - start;
	    }
	    cerr << "Resolving ops: " << elapsed_seconds_str.count() << "s" << endl;
	    // cerr << "String ops: " << elapsed_seconds_str.count() << "s" << endl;
	    // cerr << "Passed 45% cutoff: " << passed_cutoff << endl;
	    cerr << "Times needed to extend the read: " << need_to_extend_read << endl;
	}	
};

#endif