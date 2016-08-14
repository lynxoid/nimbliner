#include <iostream>
#include <fstream>
#include <random>
#include <cassert>

#include <tclap/CmdLine.h>

#include "FastaReader.h"

using namespace std;

struct input_parameters {
	string path;
	string output_path;
	int R;		// read length to sample
	int N;		// number of reads to sample
	float m = 0.0f;	// mismatch rate
	float d = 0.0f;	// indel rate
};

void print_parameters(const input_parameters & p) {
	cerr << "Input: " << p.path << endl;
	cerr << "Output: " << p.output_path << endl;
	cerr << "Read length: " << p.R << endl;
	cerr << "Read count: " << p.N << endl;
	cerr << "Mismatch rate: " << p.m << endl;
	cerr << "Indel rate: " << p.d << endl;
}

////////////////////////////////////////////////////////////////
// parse a fasta file and return vector of reads
////////////////////////////////////////////////////////////////
vector<string> parseFasta(string const & path) {
    vector<string> reads;
    FastaReader fr(path.c_str());
    kseq_t * seq;
    size_t cnt = 0;
    while ( (seq = fr.nextSequence() ) ) {
      // cerr << seq->seq.s << endl;
      reads.push_back(seq->seq.s);
      cnt++;
    }
    cerr << "(Parsed " << cnt << " input sequence) " << endl;
    return reads;
}

const char alphabet[] = {'A', 'C', 'G', 'T'};

int indexOf(const char c) {
	switch(c) {
		case 'A': return 0;
		case 'C': return 1;
		case 'G': return 2;
		case 'T': return 3;
	}
	return -1;
}

char get_next_base(const char base) {
	return alphabet[ (indexOf(base) + 1) % 4 ];
}

void select_next_base(string & read, const int j) {
	char modified_base = get_next_base(read[j]);
	read[j] = modified_base;
}

vector<int> add_mismatches(string & read, const float rate) {
	std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution d(rate);
	// flip a coin for every base w/ probability mismatch_rate
	vector<int> flipped;
	for (int j = 0; j < read.size(); j++) {
		if (d(gen)) {// returns true w/ probability rate
			// modify this base
			select_next_base(read, j);
			flipped.push_back(j);
		}
	}
	return flipped;
}

void insert_base(string & read, const int j) {
	read.erase(j, 1);
}

vector<int> add_indels(string & read, const float rate) {
	std::random_device rd;
    std::mt19937 gen(rd());
    std::bernoulli_distribution d(rate);
    vector<int> flipped;
    bool insertion = false;
    for (int j = 0; j < read.size(); j++) {
    	// toss a coin for each base -- shall we remove it or not?
		if (d(gen)) {
			flipped.push_back(j);
			if (insertion) {
				// insert a new base to the read at pos j
				read.insert(j, 1, get_next_base(read[j]) );
			}
			else {
				// delete a base from the read at pos j
				read.erase(j, 1);
			}
			insertion = !insertion;
		}
	}
    return flipped;
}

void sample_reads(const input_parameters & ip) {
	// parse a FASTA file containing the sequence
	auto chromosomes = parseFasta(ip.path);
	assert(chromosomes.size() > 0);
	string chr = chromosomes[0];
	assert(chr.size() > 0);

	ofstream r_out(ip.output_path);
	// if (ip.output_path.size() == 0)
		// r_out = cout;
	// else
		// r_out = open(ip.output_path);
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distro(0, chr.size() - 1 - ip.R );
	for (int i = 0; i < ip.N; i++) {
		// sample the read -- 
		// get a random number in [0, |chr| - R] range -- begining of the read
		int start = distro(gen);
		string read = chr.substr(start, ip.R);
		// write read name
		r_out << ">" << i << "_" << start;

		if (ip.m > 0) {
			// modify the read and save locations of errors (write as a read name)
			auto modified_bases = add_mismatches(read, ip.m);
			// write out modified mismatch bases
			if (modified_bases.size() > 0) {
				r_out << "_m=";
				for (auto const & b : modified_bases)
					r_out << b << "_";
			}
		}
		if (ip.d > 0) {
			// add / remove bases (indels)
			auto modified_bases = add_indels(read, ip.d);
			if (modified_bases.size() > 0) {
				r_out << "_i=";
				for (auto const & b : modified_bases)
					r_out << b << "_";	
			}
		}
		r_out << endl << read << endl;
	}
	// if (ip.output_path.size() > 0)
		r_out.close();
	cerr << "Sampled reads written to " << ip.output_path << endl;
}

input_parameters parse_arguments(const int argc, char * argv []) {
	TCLAP::CmdLine cmd("Sample reads w/ various error models", ' ', "0.1");
	TCLAP::ValueArg<std::string> input("i","input",
		"Input sequence to sample (e.g. chromosome sequence)", 
		true, "homer", "string");
	cmd.add( input );

	TCLAP::ValueArg<std::string> output("o","output",
		"Path to a file where we will write the reads", 
		false, "sampled_reads.fa", "string");
	cmd.add( output );

	TCLAP::ValueArg<int> rlen("l","length","Read length to sample", 
		true, 100, "int");
	cmd.add( rlen );

	TCLAP::ValueArg<int> number("n","ncount","Number of reads to sample", 
		true, 1000, "int");
	cmd.add( number );

	TCLAP::ValueArg<float> mm_rate("m","mismatch",
		"Mismatch rate per 100 nucleotides",
		false, 1.0f, "float");
	cmd.add( mm_rate );

	TCLAP::ValueArg<float> indel_rate("d", "indels",
		"Indel rate per 100 nucleotides", 
		false, 0.1f, "float");
	cmd.add( indel_rate );

	cmd.parse( argc, argv );

	// Get the value parsed by each arg. 
	input_parameters ip;
	ip.path = input.getValue();
	
	ip.R = rlen.getValue();
	ip.N = number.getValue();
	ip.m = mm_rate.getValue() / 100.0f;
	ip.d = indel_rate.getValue() / 100.0f;
	ip.output_path = output.getValue();

	print_parameters(ip);

	return ip;
}

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
int main(int argc, char * argv []) {

	// parse arguments
	input_parameters ip = parse_arguments(argc, argv);
	sample_reads(ip);
}
