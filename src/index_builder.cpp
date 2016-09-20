#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <chrono>

#include "FastaReader.h"
#include "SeqBFUtil.hpp"

// #include "reference_index.hpp"
#include "reference_index_builder.hpp"
#include "aligner.hpp"
#include "definitions.hpp"

using namespace std;

////////////////////////////////////////////////////////
//
// Improvements: 
//	- skip flanking kmers (expect low quals on the ends of the reads)
//
////////////////////////////////////////////////////////
int main(int argc, char * argv []) {
	int K = stoi(argv[1]);
	string path = argv[2];

	ReferenceIndexBuilder index;
	index.buildIndex(path, K);
}