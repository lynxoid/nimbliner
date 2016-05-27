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
	assert(chromosomes.size() > 0);
	string chr = chromosomes[0];

	// ofstream r_out("/data/chr10_reads.fa");
	std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<> dis(0, chr.size() - 1 - R );
	for (int i = 0; i < N; i++) {
		// get a random number in [0, |chr| - R] range
		int start = dis(gen);
		// r_out << ">" << i << endl << chr.substr(start, R) << endl;
		cout << ">" << i << "_" << start << endl << chr.substr(start, R) << endl;
	}
	// r_out.close();
}
