#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include "FastaReader.h"
#include "SeqBFUtil.hpp"

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
BaseBloomFilter readKmers_bf_faster(const string & path, int K) {
	ifstream in(path);
	if (!in) {
		cerr << "[ERROR] Could not open the file: " << path << endl;
		exit(1);
	}
	string line;
	getline(in, line);
	uint64_t kmer_count = stol(line);
	// cerr << "Expected kmer count: " << kmer_count << endl;
	BaseBloomFilter bf(K, kmer_count * 10);

	static const auto BUFFER_SIZE = (size_t)pow(2,22); // do not overwhelm the stack :)
	// cerr << "Buffer size, bytes: " << BUFFER_SIZE << endl;
    int fd = open(path.c_str(), O_RDONLY);
    if (fd == -1) {
    	cerr << "[ERROR] Could not open file: " << path << endl;
    	exit(1);
    }
 //    /* Advise the kernel of our access pattern.  */
 //    // posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL
    char buf[BUFFER_SIZE + 1];
    size_t i = 0, trailing = 0;
    size_t bytes_read = read(fd, buf, BUFFER_SIZE);
    size_t total_bytes_read = bytes_read;
    while (bytes_read != 0) {
        if (bytes_read == (size_t)-1) {
        	cerr << "Could not read" << endl;
        	exit(1);
        }
        char *p = buf;
        trailing = 0;
        while ( (p - buf) < bytes_read) {
        	char * p_next = (char*) memchr(p, '\n', (buf + bytes_read) - p);
        	if (p_next == NULL) {

        		trailing = buf + bytes_read - p;
        		for (int j = 0; j < trailing; j++)
        			buf[j] = *(p + j);
        		break;
        	}
        	else {
	        	if (i != 0) {
	        		uint64_t bin_kmer = strtol(p, &p_next, 10);
	        		// if (bytes_read < BUFFER_SIZE/2)
	        			// cerr << (p - buf) << " " << bin_kmer << endl;
	        		bf.add(bin_kmer);
	        	}
	        	p = p_next + 1;

	        	if (i % 5000000 == 0) cerr << i/1000000 << "m ";
	        	if (i > 1 + kmer_count) {
	        		cerr << "More kmers than expected: " << kmer_count << " vs. " << i << " bytes read: " << total_bytes_read << endl;
	        		cerr << (p-buf) << " " << bytes_read << endl;
	        		exit(1);
	        	}
	        	i++;
	        }
        }
        bytes_read = read(fd, buf + trailing, BUFFER_SIZE - trailing);
        // cerr << bytes_read << "b ";
        total_bytes_read += bytes_read;
    }
	close(fd);
	cerr << endl << "read " << (i-1) << " kmers" << endl;
	// cerr << "Trailing: " << trailing << " " << bytes_read << endl;
	if (i < kmer_count) {
		cerr << "Expected: " << kmer_count << endl;
		// exit(1);
	}
	return bf;
}

////////////////////////////////////////////////////////
unordered_map<kmer_t, vector<int>> readStarLocations(const string & path, int K) {
	unordered_map<kmer_t, vector<int>> kmer_locations;
	ifstream in(path);
	string line;
	while (getline(in, line)) {
		// parse the line now
		std::stringstream ss(line);
		string kmer, pos;
		getline(ss, kmer, ' ');
		uint64_t bin_kmer = stol(kmer);
		kmer_locations.emplace(bin_kmer, vector<int>{});
		while(getline(ss, pos, ' ')) {
    		int location = stoi(pos);
    		kmer_locations[bin_kmer].push_back(location);
		}
		if (kmer_locations[bin_kmer].size() > 10000) {
			cerr << mer_binary_to_string(bin_kmer, K) << " " << kmer_locations[bin_kmer].size() << endl;
			kmer_locations.erase( kmer_locations.find(bin_kmer) );
		}
	}
	in.close();
	return kmer_locations;
}

////////////////////////////////////////////////////////
// build index stage
////////////////////////////////////////////////////////
void getAllKmersAndStars(const string & ref_path, int K) {
	auto chromosomes = parseFasta(ref_path);
	// count kmers, record their locations
	unordered_map<kmer_t, vector<int>> kmer_locations;

	unordered_set<kmer_t> star_kmers;
	assert(chromosomes.size() > 0);
	auto chr = chromosomes[0];

	int c = 0;
	for (int i = 0; i < chr.size() - K + 1; i++) {
		kmer_t kmer = mer_string_to_binary(&chr[i], K);
		if ( kmer_locations.find(kmer) == kmer_locations.end() ) {
			kmer_locations.emplace(kmer, vector<int>{i});
		}
		else {
			 // can delta encode here and fit into less space technically
			kmer_locations[kmer].push_back(i);
		}
		// pick stars at every 50 bases
		if (i % 50 == 0) star_kmers.emplace(kmer);
		if (i % 1000000 == 0) cerr << i/1000000 << "Mbp ";
	}
	cerr << "(" << kmer_locations.size() << " kmers)" << endl;
	cerr << "(" << star_kmers.size() << " star kmers)" << endl;


	ofstream star_locations_out("star_locations.txt");
	ofstream all_kmers("all_kmers.txt");
	all_kmers << kmer_locations.size() << endl;
	cerr << "saving star locations" << endl;
	for (auto star : star_kmers) {
		star_locations_out << star << " ";
		for (auto loc : kmer_locations[star]) star_locations_out << loc << " ";
		star_locations_out << endl;
		all_kmers << star << endl;
		kmer_locations.erase(star);
	}
	cerr << "saving all kmers" << endl;
	for (auto p : kmer_locations) {
		all_kmers << p.first << endl;
	}
	all_kmers.close();
	star_locations_out.close();
	kmer_locations.clear();
}

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

	// what if only one star?
	if (matched_stars.size() == 0) {
		// need to extend the read until we hit some star
		extend++;
	}
	if (matched_stars.size() == 1) {
		// that's the only thing we got going
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
        if (seq->seq.l < K) continue;
        // stream through BF
        vector<pair<kmer_t,int>> matched_stars;
        float cnt_matched = 0.0f;
        int L = seq->seq.l - K + 1;
        assert(L < 65536);
        int true_pos = stoi(seq->name.s);
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
                                matched_stars.emplace_back(bin_kmer, i);
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
		auto end = std::chrono::system_clock::now();
		elapsed_seconds_str += end - start;
    }
    cerr << "Resolving ops: " << elapsed_seconds_str.count() << "s" << endl;
    // cerr << "String ops: " << elapsed_seconds_str.count() << "s" << endl;
    cerr << "Passed 45% cutoff: " << passed_cutoff << endl;
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
		getAllKmersAndStars(path, K);
	}
	else if (mode == "query") {
		string kmers_path = argv[4];
		string stars_path = argv[5];
		auto start = std::chrono::system_clock::now();
		cerr << "reading kmers from " << kmers_path << endl;
		auto bf = readKmers_bf_faster(kmers_path, K);

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		cerr << "reading kmers took: " << elapsed_seconds.count() << "s" << endl;

		cerr << "reading stars" << endl;
		auto star_locations = readStarLocations(stars_path, K);
		cerr << "CAN NOW TEST THE MAPPING" << endl;
		end = std::chrono::system_clock::now();
		elapsed_seconds = end-start;
	    cerr << "setup: " << elapsed_seconds.count() << "s" << endl;

	    ////////////////////////////////////////////////////////////////////////
		// read one read at a time
		////////////////////////////////////////////////////////////////////////

		start = std::chrono::system_clock::now();
		queryReads(path, K, bf, star_locations);
	    end = std::chrono::system_clock::now();
	    elapsed_seconds = end-start;
	    cerr << "querying: " << elapsed_seconds.count() << "s" << endl;
		
	}
}
