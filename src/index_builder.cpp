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
	int K = stoi(argv[1]);
	string path = argv[2];

	ReferenceIndex index;
	index.buildIndex(path, K);
}