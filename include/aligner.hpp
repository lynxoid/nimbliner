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

	////////////////////////////////////////////////////////
	//	matched_stars -- in order in which they appear in the read
	//	resolve which of the star co-locations are the same distance
	//	apart as the stars in the read
	//  0123456789
	//	--***-***-
	//
	////////////////////////////////////////////////////////
	vector<int> resolve_mapping_locations(vector<pair<kmer_t,int>> & matched_stars, 
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
			// if two stars	
			auto first_star = matched_stars[0];
			auto second_star = matched_stars[1];
			int x = second_star.second - first_star.second;

			auto A = _index.get_anchor_locations(first_star.first);
			auto B = _index.get_anchor_locations(second_star.first);
			// auto A = star_locations[first_star.first];
			// auto B = star_locations[second_star.first];
			int i = 0, j = 0;
			while (i < A.size() && j < B.size() ) {
				if (A[i] < B[j]) {
					if (B[j] - A[i] == x) {
						// found one match
						i++;
						j++;
						mappings.push_back(A[i] - first_star.second);
					} 
					else if (B[j] - A[i] < x) {
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
		}
		// else {
			// what if more than 2 stars?
			// OMG, what do we do now?!!
		// }

		return mappings;
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
	    	// cerr << seq->seq.l << " " << K << endl;
	        if (seq->seq.l < K) {
	        	continue;
	        }
	        // stream through BF
	        vector<pair<kmer_t,int>> matched_stars;
	        float cnt_matched = 0.0f;
	        int L = seq->seq.l - K + 1;
	        // only allow reads shorter than 2^16
	        assert(L < 65536);
	        // int true_pos = stoi(seq->name.s);
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
	                // update prev kmer
	                // mask the leftmost character
	                kmer_t mask = ( ((kmer_t)1) << ( 2*(K-1) ) ) - 1;
	                bin_kmer = (bin_kmer & mask) << 2;
	                // append a new char on the right
	                bin_kmer = bin_kmer | dna_codes[ (size_t) seq->seq.s[i] ];
	                // print for debug
	                // unsigned short j = 0;
	                // while (j < i) {cerr << " "; j++;}
	                // cerr << mer_binary_to_string(bin_kmer, K) << endl;
	                //if (bf.contains(bin_kmer) ) {
	                //        cnt_matched++;
	                        // check if hit a star kmer
	                        if ( _index.has_anchor(bin_kmer) )
	                                matched_stars.emplace_back(bin_kmer, i - K + 1);
	                        if (matched_stars.size() >= 2) break;
	                //}
	                i++;
	        }
	        // require at least 50% of all kmers to match
	        // if (cnt_matched / L < 0.45) continue;
	        // passed_cutoff++;
			// resolve star kmers to get an exact mapping location
			auto start = std::chrono::system_clock::now();
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