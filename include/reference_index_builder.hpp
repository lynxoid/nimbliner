// reference index builder
#ifndef INDEX_BUILDER
#define INDEX_BUILDER

#include <string>
#include <vector>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <memory>
#include <algorithm>

#include <boost/timer.hpp>

#include "bit_tree_binary.hpp"
#include "definitions.hpp"

class ReferenceIndexBuilder {

	// TODO: build reference index and pick anchors based on the
	// Kingsford (de Bruijn cover), Pachter (kallisto)
	void build_index2() {

	}

	/*
	 */
	void write_anchors(unordered_set<kmer_t> & star_kmers, 
		const unordered_map<kmer_t, vector<genomic_coordinate_t>> & kmer_locations) {
		cerr << "saving anchors" << endl;

		ofstream star_locations_out("anchors.txt");
		auto start = std::chrono::system_clock::now();

		for (auto star : star_kmers) {
			star_locations_out << star << " ";
			auto it = kmer_locations.find(star);
			assert(it != kmer_locations.end());
			for (auto loc : it->second) star_locations_out << loc << " ";
			star_locations_out << endl;
			// all_kmers << star << endl;
			// kmer_locations.erase(star);
		}
		star_locations_out.close();

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds_str = end - start;
	    cerr << "Saving anchors took: " << elapsed_seconds_str.count() << "s" << endl;
	}

	/*
	 */
	void write_index(unordered_map<kmer_t, vector<genomic_coordinate_t>> & kmer_locations) {
		cerr << "saving all kmers" << endl;
		auto start = std::chrono::system_clock::now();

		ofstream all_kmers("all_kmers.txt");
		// write hte # of kmers to expect
		all_kmers << kmer_locations.size() << endl;

		int i = 0;
		while (kmer_locations.size() > 0) {
			auto it = kmer_locations.begin();
			all_kmers << it->first << endl;
			kmer_locations.erase(it);
			i++;
		}
		all_kmers.close();
		cerr << "wrote " << i << " kmers" << endl;
		assert(kmer_locations.size() == 0);

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds_str = end - start;
	    cerr << "Saving kmers took: " << elapsed_seconds_str.count() << "s" << endl;
		kmer_locations.clear();
	}

	void write_bit_tree_index(unordered_map<kmer_t, vector<genomic_coordinate_t>> & kmer_locations,
			const uint k) {
		// dump to a vector
		cerr << "dumping kmers into a vector" << endl;
		boost::timer t;
		shared_ptr<vector<kmer_t>> kmers = shared_ptr<vector<kmer_t>>(new vector<kmer_t>());
		while (kmer_locations.size() != 0) {
			auto it = kmer_locations.begin();
			kmers->push_back(it->first);
			kmer_locations.erase(it);
		}
		assert(kmer_locations.size() == 0);
		t.restart();
		cerr << "Sorting kmers" << endl;
		std::sort(kmers->begin(), kmers->end() );
		cerr << "Sorting took " << t.elapsed() << " s" << endl;
		t.restart();
		cerr << "encoding" << endl;
		BitTreeBin bit_tree;
		bit_tree.encode(kmers, k);
		cerr << "encoding took " << t.elapsed() << " s" << endl;
		t.restart();
		cerr << "Writing to a binary file" << endl;
		// bit_tree.write(input_file + ".btbin");
		bit_tree.write("index.btbin");
		cerr << "(" << t.elapsed() << " s)" << endl;
	}

public:
	ReferenceIndexBuilder() {}

	// given a path to the reference and a kmer length K, parse the reference
	// fasta, obtain all kmers and select anchors using a heuristic
	// write kmers to an *.index file and write anchors to *.star file.
	void buildIndex(const string & ref_path, const uint K) {
		auto chromosomes = parseFasta(ref_path);
		// count kmers, record their locations
		unordered_map<kmer_t, vector<genomic_coordinate_t>> kmer_locations;

		unordered_set<kmer_t> star_kmers;
		assert(chromosomes.size() > 0);
		auto chr = chromosomes[0];

		int c = 0;
		cerr << "gathering all kmers... ";
		for (genomic_coordinate_t i = 0; i < chr.size() - K + 1; i++) {
			kmer_t kmer = nimble::mer_string_to_binary(&chr[i], K);
			if ( kmer_locations.find(kmer) == kmer_locations.end() ) {
				kmer_locations.emplace(kmer, vector<genomic_coordinate_t>{i});
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
		cerr << "(" << star_kmers.size() << " anchors)" << endl;

		write_anchors(star_kmers, kmer_locations);
		// switch between different implementations of the index
		// write_index(kmer_locations);
		// bit tree representation will take less space
		write_bit_tree_index(kmer_locations, K);
		return;
	}
};

#endif // INDEX_BUILDER
