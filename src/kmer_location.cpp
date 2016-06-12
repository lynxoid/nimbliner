#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <chrono>

#include "FastaReader.h"
#include "SeqBFUtil.hpp"

#include "reference_index.hpp"
#include "aligner.hpp"
#include "definitions.hpp"

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
		ReferenceIndex index;
		index.buildIndex(path, K);
	}
	else if (mode == "query") {
		string kmers_path = argv[4];
		string stars_path = argv[5];
		
		ReferenceIndex index;
		index.readIndex(kmers_path, stars_path, K);
		// auto bf = index_reader.readKmers_bf_faster(kmers_path, K);
		// auto star_locations = index_reader.readStarLocations(stars_path, K);

		cerr << "========================" << endl;
		cerr << "CAN NOW TEST THE MAPPING" << endl;
		// end = std::chrono::system_clock::now();
		// elapsed_seconds = end-start;
	    // cerr << "setup took: " << elapsed_seconds.count() << "s" << endl;

	    ////////////////////////////////////////////////////////////////////////
		// read one read at a time
		////////////////////////////////////////////////////////////////////////

		auto start = std::chrono::system_clock::now();
		Aligner aligner(index);
		// TODO: separate sequence reads and aligner -- make aligner pull things off the queue
		aligner.alignReads(path, K);
	    auto end = std::chrono::system_clock::now();
	    std::chrono::duration<double> elapsed_seconds = end - start;
	    cerr << "querying: " << elapsed_seconds.count() << "s" << endl;
		
	}
}
