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

namespace nimble {

inline bool can_read_or_quit(const ifstream & in, const string & name, const bool quit = true) {
    if (!in) {
        cerr << "[ERROR] Can not read from file " << name << endl;
        if (quit)
            exit(1);
        return false;
    }
    return true;
}

inline void can_write_or_quit(const ifstream & in) {

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

    shared_ptr<unordered_map<kmer_t, uint8_t>> build_pdBG(const string & fofn_path, const uint K) {
        // open a fofn file, read one line at a time, parse fasta seqeunces at
        // every line
        ifstream fofn_in(fofn_path);
        // quit if cant read from this file
        can_read_or_quit(fofn_in, fofn_path, true);

    	// TODO: use counting BF/perfect hash here
        shared_ptr<unordered_map<kmer_t, uint8_t>> kmer_counts(new unordered_map<kmer_t, uint8_t>() );

        string fasta_path;
        while (getline(fofn_in, fasta_path)) {
            FastaReader fr(fasta_path.c_str());
            if (!fr.is_open() ) {
                cerr << "[ERROR] Can not read from " << fasta_path <<
                    ". Skipping..." << endl;
                continue;
            }
            kseq_t * seq;
            while ( (seq = fr.nextSequence() ) ) {
                string chr = seq->seq.s;
                // get all kmers from this sequence
        		for (genomic_coordinate_t i = 0; i < chr.size() - K + 1; i++) {
                    // TODO: use a rolling binary transform here to get the next kmer
        			kmer_t kmer = nimble::mer_string_to_binary(&chr[i], K);
        			if ( kmer_counts->find(kmer) == kmer_counts->end() ) {
                        (*kmer_counts)[kmer] = 1;
        			}
        			else {
        				// can delta encode here and fit into less space technically
        				if ( (*kmer_counts)[kmer] < 255)
        					(*kmer_counts)[kmer]++;
        				// else -- do not increment to avoid overflow
        			}
        			if (i % 1000000 == 0) cerr << i/1000000 << "Mbp ";
        		}
                cerr << "(" << kmer_counts->size() << " kmers) ";
            }
        }
        cerr << endl;
        fofn_in.close();
        return kmer_counts;
    }

    shared_ptr<unordered_map<kmer_t, list<seed_position_t> > > select_anchors(
        const string & fofn_path, const uint K, uint freq_cutoff,
        shared_ptr<unordered_map<kmer_t, uint8_t>> kmer_counts) {
        ifstream fofn_in(fofn_path);
        can_read_or_quit(fofn_in, fofn_path, true);
        shared_ptr<unordered_map<kmer_t, list<seed_position_t> > > anchors(
            new unordered_map<kmer_t, list<seed_position_t> >());

        // TODO: expose this variable as a parameter
        int offset = 50;

        reference_id_t ref_id = 0;
        string fasta_path;
        while (getline(fofn_in, fasta_path)) {
            FastaReader fr(fasta_path.c_str());
            if (!fr.is_open() ) {
                cerr << "[ERROR] Can not read from " << fasta_path <<
                    ". Skipping..." << endl;
                continue;
            }
            kseq_t * seq;
            while ( (seq = fr.nextSequence() ) ) {
                string chr = seq->seq.s;
                for (genomic_coordinate_t i = 0; i < chr.size() - K + 1; i++) {
                    // TODO: use a rolling approach to get next kmer
        			kmer_t kmer = nimble::mer_string_to_binary(&chr[i], K);
        			if ( (*kmer_counts)[kmer] < freq_cutoff) {
        				// add to anchors
        				if (anchors->find(kmer) == anchors->end() ) {
        					anchors->emplace(kmer, list<seed_position_t>{{ref_id,i}});
        				}
        				else
        					(*anchors)[kmer].emplace_back(ref_id, i);
        				i += offset - 1;
        			}
        			if (i % 1000000 == 0) cerr << i/1000000 << "Mbp ";
        		}
                // increment reference sequence index
                ref_id++;
                cerr << "(" << anchors->size() << " anchors from " << ref_id << " reference sequences)" << endl;
            }
        }

        return anchors;
    }

    /*
     * build index in 2 passes
     */
	void buildIndex(const string & fofn_path, const string & output_prefix,
        const uint K) {

        cerr << "Pass 1: building background pdBG" << endl;
        shared_ptr<unordered_map<kmer_t, uint8_t>> kmer_counts =
            build_pdBG(fofn_path, K);


		cerr << "Pass 2: selecting anchors" << endl;

        // TODO: expose as a paramter
		const int x = 5;

		// TODO: can keep linked lists since expect these lists to be short
        // is LL less overhead than vector?
		shared_ptr<unordered_map<kmer_t, list<seed_position_t> > > anchors =
            select_anchors(fofn_path, K, x, kmer_counts);


        // switch between different implementations of the index
        // nimble::ReferenceIndex::write_index(kmer_counts);
		nimble::BloomReferenceIndex::write_index(kmer_counts, output_prefix);
        // TODO: prune kmer_counts as we write them only keeping those that are
        // below frequency cutoff

		nimble::AnchorIndex::write_anchors(anchors, output_prefix);
		// bit tree representation will take less space
		// write_bit_tree_index(kmer_locations, K);
		return;
	}

};

}

#endif // INDEX_BUILDER
