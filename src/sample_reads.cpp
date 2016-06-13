#include <iostream>
#include <fstream>
#include <random>
#include <cassert>
#include "FastaReader.h"

using namespace std;

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
    cerr << "(" << cnt << " reads) ";
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

void select_next_base(string & read, const int j) {
	char base = read[j];
	char modified_base = alphabet[ (indexOf(base) + 1) % 4 ];
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

////////////////////////////////////////////////////////////////
//
////////////////////////////////////////////////////////////////
int main(int argc, char * argv []) {
	// read length
	int R = stoi(argv[1]);
	// number of reads to sample
	int N = stoi(argv[2]);
	// reference sequence to sample from
	auto chromosomes = parseFasta(argv[3]);

	float mismatch_rate = 0;
	if (argc > 4) {
		// mismatch rate [0;100], float
		mismatch_rate = stof(argv[4]) / 100.0f;
	}

	assert(chromosomes.size() > 0);
	// assert(mismatch_rate < 1.0f);
	string chr = chromosomes[0];

	// ofstream r_out("/data/chr10_reads.fa");
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> distro(0, chr.size() - 1 - R );
	for (int i = 0; i < N; i++) {
		// get a random number in [0, |chr| - R] range -- begining of the read
		int start = distro(gen);
		string read = chr.substr(start, R);
		// cerr << "get read" << endl;
		
		// write to stdout
		cout << ">" << i << "_" << start;

		if (mismatch_rate > 0) {
			// modify the read and save locations of errors (write as a read name)
			auto modified_bases = add_mismatches(read, mismatch_rate);
			// write out modified bases
			for (auto const & b : modified_bases)
				cout << "_" << b;
		}
		cout << endl << read << endl;
	}
	// r_out.close();
}
