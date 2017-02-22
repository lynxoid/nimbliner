#include <iostream>
#include <fstream>
#include <cassert>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <chrono>

#include <tclap/CmdLine.h>
#include "FastaReader.h"
// #include "SeqBFUtil.hpp"

// #include "reference_index.hpp"
#include "reference_index_builder.hpp"
#include "aligner.hpp"
#include "definitions.hpp"

using namespace std;

struct input_parameters {
	int K;
	string fasta;
    string output;
};

input_parameters parse_args(int argc, char * argv[]) {
	TCLAP::CmdLine cmd("Build an index given the reference sequence", ' ', "0.1");

    TCLAP::ValueArg<std::string> input("i","input",
                "Reference sequence (fasta or fastq)",
                true, " ", "string");
    cmd.add( input );
    TCLAP::ValueArg<std::string> output("o","output",
                "Use as output prefix",
                false, " ", "string");
    cmd.add( output );

	TCLAP::ValueArg<int> klen("k","kmer-length","Seed length to sample",
                true, 20, "int");
    cmd.add( klen );

	cmd.parse(argc, argv);
    input_parameters ip;
    if (output.isSet() )
        ip.output = output.getValue();
    else
        ip.output = input.getValue();
	ip.K = klen.getValue();
	ip.fasta = input.getValue();
	return ip;
}

////////////////////////////////////////////////////////
//
// Improvements:
//	- skip flanking kmers (expect low quals on the ends of the reads)
//
////////////////////////////////////////////////////////
int main(int argc, char * argv []) {
	auto ip = parse_args(argc, argv);
	ReferenceIndexBuilder index;
	index.buildIndex(ip.fasta, ip.output, ip.K);
}
