#include <iostream>
#include <fstream>

#include <boost/progress.hpp>

#include "bit_tree_binary.hpp"

using namespace std;

//////////////////////////////////////////////////////////////////////////////////////////
// 
//////////////////////////////////////////////////////////////////////////////////////////
shared_ptr<vector<kmer_type>> parseStrings(string const filename, int & k) {
	ifstream in_file(filename);
	if (!in_file) {
		cerr << "File does not exist:" << filename << endl;
		exit(0);
	}
	string kmer;
	shared_ptr<vector<kmer_type>> kmers(new vector<kmer_type>());
	getline(in_file, kmer);
	long expected_count = stol(kmer);
	cerr << "Expected count: " << expected_count << endl;
	k = 20;
	while ( getline(in_file, kmer) ) {
		kmer_type bin_kmer = stol(kmer);
		// k = kmer.size();
		// kmers->push_back( mer_string_to_binary(kmer.c_str(), kmer.size()) );
		kmers->push_back(bin_kmer);
	}
	return kmers;
}

/****************************************************************
**
**
**
****************************************************************/
int main(int argc, char * argv []) {
	const string mode = argv[1];
	string input_file = argv[2];
	bool unsorted = (stoi(argv[3]) == 1);
	int k = 0;
	shared_ptr<vector<kmer_type>> kmers(new vector<kmer_type>());

	BitTreeBin bit_tree;
	// DAWGS: https://code.google.com/p/dawgdic/
	if (mode.compare("encode") == 0) {
		cerr << "encode" << endl << "Reading kmers ";
		{
		boost::timer t;
		kmers = parseStrings(input_file, k);
		cerr << "(" << t.elapsed() << " s)" << endl;
		}
		cerr << "Found " << kmers->size() << " kmers" << endl;
		assert(k > 0 && k <= 32);

		assert(kmers->size() > 0);
		cerr << "Read " << kmers->size() << " strings" << endl;
		cerr << "k = " << k << endl;
		if (unsorted) { // sort the little guys
			cerr << "Sorting kmers" << endl;
			std::sort(kmers->begin(), kmers->end());
		}


		boost::timer t;
		cerr << "encoding" << endl;
		bit_tree.encode(kmers, k);
		cerr << "(" << t.elapsed() << " s)" << endl;
		t.restart();
		cerr << "Writing to a binary file" << endl;
		bit_tree.write(input_file + ".btbin");
		cerr << "(" << t.elapsed() << " s)" << endl;
	}
	else if (mode.compare("decode") == 0 ) {
		cerr << "decode mode " << input_file << endl;
		boost::timer t;
		bit_tree.read(input_file, k);
		cerr << "Reading took: " << t.elapsed() << endl;
		assert(k > 0 && k <= 32);
		t.restart();

		// TEST 1000000mln kmers for set membership


		kmers = bit_tree.decode();
		cerr << "(" << t.elapsed() << " s)" << endl;
		cerr << "Recovered " << kmers->size() << " kmers" << endl;
		t.restart();
		cerr << "Writing kmers to stdin... ";
		for (auto kmer : *kmers) {
			// cout << mer_binary_to_string(kmer, k) << endl;
			cout << kmer << endl;
		}
		cerr << "(" << t.elapsed() << " s)" << endl;
	}
	else {
		cerr << "Unknown mode: " << mode << ". Stopping." << endl;
	}
}
