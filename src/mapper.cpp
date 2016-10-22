#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <chrono>
#include <memory>

#include <tclap/CmdLine.h>

// #include "FastaReader.h"
#include "SeqBFUtil.hpp"

#include "reference_index.hpp"
#include "bloom_reference_index.hpp"
#include "bit_tree_index.hpp"
#include "aligner.hpp"
// #include "parallel_aligner.hpp"
#include "definitions.hpp"

using namespace std;

struct input_parameters {
	string input_fasta;
	string input_anchors;
	string input_index;
	string output_path;
	// TODO: embed this into the index so that user does not have to enter it
	int K;		// kmer size to use (assumes this K was used in index)
};

input_parameters parse_arguments(const int argc, char * argv []) {
	TCLAP::CmdLine cmd("Align reads faster and in less memory", ' ', "0.1");
	TCLAP::ValueArg<std::string> input("i","input",
		"Input reads (fasta or fastq)", true, "?", "string");
	cmd.add( input );

	TCLAP::ValueArg<std::string> output("o","output",
		"Filename to write alignments to", false, "?", "string");
	cmd.add( output );
	cmd.parse( argc, argv );

	// TCLAP::ValueArg<int> output("k","kmer-length",
		// "Kmer legnth to use for ", false, "?", "string");
	// cmd.add( output );
	// cmd.parse( argc, argv );

	// Get the value parsed by each arg.
	input_parameters ip;
	ip.input_fasta = input.getValue();
	// ip.input_anchors = anchors.getValue();
	// ip.input_index = index.getValue();
	ip.K = 20;
	ip.output_path = output.getValue();

	return ip;
}

////////////////////////////////////////////////////////
//
// Improvements:
//	- afford for errors --- jump in de Bruijn graph? test for all 4*k variants of the string
//  - pick stars in a principled way
//	- if no stars matched -- extend the (read minus flanking seq) on either side until hit a star
//	- handle reverse-complimented reads
//
////////////////////////////////////////////////////////
int main(int argc, char * argv []) {
	// input_parameters ip = parse_arguments(argc, argv);

	// string mode = argv[1];
	int K = stoi(argv[1]);
	string path = argv[2];

	// if (mode == "index") {
		// ReferenceIndex index;
		// index.buildIndex(path, K);
	// }
	// else if (mode == "query")
	{
		string kmers_path = argv[3];
		string stars_path = argv[4];

		shared_ptr<ReferenceIndex> index = shared_ptr<BloomReferenceIndex>(new BloomReferenceIndex() );
		// shared_ptr<ReferenceIndex> index = shared_ptr<BitTreeIndex>(new BitTreeIndex() );
		index->readIndex(kmers_path, stars_path, K);

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
		// ParallelAligner aligner(index);
		// TODO: separate sequence reads and aligner -- make aligner pull things off the queue
		aligner.alignReads(path, K, false /* debug */ );
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end - start;
		cerr << "querying: " << elapsed_seconds.count() << "s" << endl;
	}
}
