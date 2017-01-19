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
#include <list>

#include <boost/timer.hpp>

#include "FastaReader.h"
#include "anchor_index.hpp"
#include "bloom_reference_index.hpp"
#include "bit_tree_binary.hpp"
#include "definitions.hpp"

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

class ReferenceIndexBuilder {

	// TODO: build reference index and pick anchors based on the
	// Kingsford (de Bruijn cover), Pachter (kallisto)
	void build_index2() {

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
		cerr << "Dump to vector: " << t.elapsed() << "s" << endl;
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
		cerr << "Writing BitTree took " << t.elapsed() << " s" << endl;
	}

public:
	ReferenceIndexBuilder() {}

	// given a path to the reference and a kmer length K, parse the reference
	// fasta, obtain all kmers and select anchors using a heuristic
	// write kmers to an *.index file and write anchors to *.star file.
	// void buildIndex(const string & ref_path, const uint K) {
	// 	auto chromosomes = parseFasta(ref_path);
	// 	// count kmers, record their locations
	// 	unordered_map<kmer_t, vector<genomic_coordinate_t>> kmer_locations;
	//
	// 	unordered_set<kmer_t> star_kmers;
	// 	assert(chromosomes.size() > 0);
	// 	auto chr = chromosomes[0];
	//
	// 	int c = 0;
	// 	cerr << "gathering all kmers... ";
	// 	for (genomic_coordinate_t i = 0; i < chr.size() - K + 1; i++) {
	// 		kmer_t kmer = nimble::mer_string_to_binary(&chr[i], K);
	// 		if ( kmer_locations.find(kmer) == kmer_locations.end() ) {
	// 			kmer_locations.emplace(kmer, vector<genomic_coordinate_t>{i});
	// 		}
	// 		else {
	// 			 // can delta encode here and fit into less space technically
	// 			kmer_locations[kmer].push_back(i);
	// 		}
	// 		// pick stars at every 50 bases
	// 		if (i % 50 == 0) star_kmers.emplace(kmer);
	// 		if (i % 1000000 == 0) cerr << i/1000000 << "Mbp ";
	// 	}
	// 	cerr << "(" << kmer_locations.size() << " kmers)" << endl;
	// 	cerr << "(" << star_kmers.size() << " anchors)" << endl;
	//
	// 	write_anchors(star_kmers, kmer_locations);
	// 	// switch between different implementations of the index
	// 	write_index(kmer_locations);
	// 	// bit tree representation will take less space
	// 	// write_bit_tree_index(kmer_locations, K);
	// 	return;
	// }

    /*
     * build index in 2 passes
     */
	void buildIndex(const string & ref_path, const uint K) {
        auto chromosomes = parseFasta(ref_path);

		// naive counter
		// TODO: use counting BF
		unordered_map<kmer_t, uint8_t> kmer_counts;
		assert(chromosomes.size() > 0);
		// TODO: index all chromosomes given as input
		auto chr = chromosomes[0];

		int c = 0;
		cerr << "Gathering kmers: pass 1" << endl;
		for (genomic_coordinate_t i = 0; i < chr.size() - K + 1; i++) {
			kmer_t kmer = nimble::mer_string_to_binary(&chr[i], K);
			if ( kmer_counts.find(kmer) == kmer_counts.end() ) {
					kmer_counts[kmer] = 1;
			}
			else {
				// can delta encode here and fit into less space technically
				if (kmer_counts[kmer] < 255)
					kmer_counts[kmer]++;
				// else -- do not increment to avoid overflow
			}
			if (i % 1000000 == 0) cerr << i/1000000 << "Mbp ";
		}
		cerr << "(" << kmer_counts.size() << " kmers)" << endl;
		cerr << "Gathering kmers: pass 2" << endl;

		const int x = 5;

		// TODO: can keep linked lists since expect these lists to be short
		unordered_map<kmer_t, list<genomic_coordinate_t>	> anchors;
		for (genomic_coordinate_t i = 0; i < chr.size() - K + 1; i++) {
			kmer_t kmer = nimble::mer_string_to_binary(&chr[i], K);
			if (kmer_counts[kmer] < x) {
				// add to anchors
				if (anchors.find(kmer) == anchors.end() ) {
					anchors.emplace(kmer, list<genomic_coordinate_t>{i});
				}
				else
					anchors[kmer].push_back(i);
				i += 50;
			}
			if (i % 1000000 == 0) cerr << i/1000000 << "Mbp ";
		}

		cerr << "(" << anchors.size() << " anchors)" << endl;

		nimble::AnchorIndex::write_anchors(anchors);
		// switch between different implementations of the index
        // nimble::ReferenceIndex::write_index(kmer_counts);
		nimble::BloomReferenceIndex::write_index(kmer_counts);
		// bit tree representation will take less space
		// write_bit_tree_index(kmer_locations, K);
		return;
	}

};

#endif // INDEX_BUILDER
