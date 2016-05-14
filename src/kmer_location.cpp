#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <chrono>

#include "FastaReader.h"
#include "SeqBFUtil.hpp"

#include "IndexBuilder.hpp"
#include "IndexReader.hpp"

using namespace std;

////////////////////////////////////////////////////////
// BaseBloomFilter readKmers_bf(const string & path, int K) {
// 	BaseBloomFilter bf(K, 113590577 * 10);
// 	ifstream in(path);
// 	string kmer;
// 	int i = 0;
// 	while (getline(in, kmer)) {
// 		// cerr << kmer << endl;
// 		// cerr << i++ << " ";
// 		uint64_t bin_kmer = stol(kmer);
// 		bf.add(bin_kmer);
// 		i++;
// 	}
// 	in.close();
// 	cerr << "read " << i << " kmers" << endl;
// 	return bf;
// }


////////////////////////////////////////////////////////
//	matched_stars -- in order in which they appear in the read
//	resolve which of the star co-locations are the same distance
//	apart as the stars in the read
//  0123456789
//	--***-***-
//
////////////////////////////////////////////////////////
vector<int> resolve_mapping_locations(vector<pair<kmer_t,int>> & matched_stars, 
	unordered_map<kmer_t, vector<int>> & star_locations,
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
		auto locations = star_locations[matched_stars[0].first];
		for (auto loc : locations)
			mappings.push_back(loc - matched_stars[0].second);
	}
	else if (matched_stars.size() >= 2) {
		// if two stars	
		auto first_star = matched_stars[0];
		auto second_star = matched_stars[1];
		int x = second_star.second - first_star.second;

		auto A = star_locations[first_star.first];
		auto B = star_locations[second_star.first];
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

////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////
void queryReads(const string & path, int K, BaseBloomFilter & bf,
	unordered_map<kmer_t, vector<int>> & star_locations) {
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
                if (star_locations.find(bin_kmer) != star_locations.end() )
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
                        if (star_locations.find(bin_kmer) != star_locations.end() )
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
		auto mapping_locations = resolve_mapping_locations(matched_stars, 
			star_locations, need_to_extend_read);

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

////////////////////////////////////////////////////////
//
// Improvements: 
//	- skip flanking kmers (expect low quals on the ends of the reads)
//	- afford for errors --- jump in de Bruijn graph? test for all 4*k variants of the string
//  - pick stars in a principled way
//	- if no stars matched -- extend the (read minus flanking seq) on either side until hit a star
//	- handle reverse-complimented reads
//
////////////////////////////////////////////////////////
int main(int argc, char * argv []) {
	string mode = argv[1];
	int K = stoi(argv[2]);
	string path = argv[3];

	if (mode == "index") {
		IndexBuilder index;
		index.getAllKmersAndStars(path, K);
	}
	else if (mode == "query") {
		string kmers_path = argv[4];
		string stars_path = argv[5];
		
		IndexReader index_reader;
		auto bf = index_reader.readKmers_bf_faster(kmers_path, K);
		auto star_locations = index_reader.readStarLocations(stars_path, K);

		cerr << "========================" << endl;
		cerr << "CAN NOW TEST THE MAPPING" << endl;
		// end = std::chrono::system_clock::now();
		// elapsed_seconds = end-start;
	    // cerr << "setup took: " << elapsed_seconds.count() << "s" << endl;

	    ////////////////////////////////////////////////////////////////////////
		// read one read at a time
		////////////////////////////////////////////////////////////////////////

		auto start = std::chrono::system_clock::now();
		queryReads(path, K, bf, star_locations);
	    auto end = std::chrono::system_clock::now();
	    std::chrono::duration<double> elapsed_seconds = end - start;
	    cerr << "querying: " << elapsed_seconds.count() << "s" << endl;
		
	}
}
