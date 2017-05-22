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

// #include <boost/timer.hpp>

#include "bf/bloom_filter/counting.h"
#include "bf/hash.h"

#include "FastaReader.h"
#include "anchor_index.hpp"
#include "bloom_reference_index.hpp"
// #include "bit_tree_binary.hpp"
#include "fofn_reader.hpp"
#include "definitions.hpp"

namespace nimble {

    void increment_map(unordered_map<kmer_t, uint32_t> & m,
        const bin_kmer_t kmer) {
            if (m.find(kmer) == m.end()) {
                m.emplace(kmer, 1);
            }
            else
                m[kmer]++;
    }

    void add_to_map(shared_ptr<unordered_map<kmer_t, list<seed_position_t> > > anchors,
        const bin_kmer_t kmer, const reference_id_t ref_id, const genomic_coordinate_t pos) {
        if (anchors->find(kmer) == anchors->end() ) {
            anchors->emplace(kmer, list<seed_position_t>{{ref_id,pos}});
        }
        else
            (*anchors)[kmer].emplace_back(ref_id, pos);
    }


class ReferenceIndexBuilder {

	// TODO: build reference index and pick anchors based on the
	// Kingsford (de Bruijn cover), Pachter (kallisto)

public:
	ReferenceIndexBuilder() {}

    // open a fofn file, read one line at a time, parse fasta seqeunces at
    // every line
    // shared_ptr<unordered_map<kmer_t, uint8_t>>
    shared_ptr<bf::counting_bloom_filter> build_pdBG(const string & fofn_path, const uint K) {
        FOFNReader fofn_reader(fofn_path);
        // shared_ptr<unordered_map<kmer_t, uint8_t>> kmer_counts(new unordered_map<kmer_t, uint8_t>() );

        // we only care about very rare kmers, so anything above 5 will be
        // considered frequent -> ceil(log_2 (5)) = 3
        uint counting_bits = 3;
        // may expect around 3 * 10^9 kmers
        uint64_t elements = 3 * 10^9;
        shared_ptr<bf::counting_bloom_filter> cbf = shared_ptr<bf::counting_bloom_filter>(new
                bf::counting_bloom_filter(bf::make_hasher(2), /* hasing fun-ns */
                                          elements,           /* expected # elem */
                                          counting_bits) );   /* bit per counter */

        kseq_t * seq;
        while ( (seq = fofn_reader.getNextSequence() ) ) {
            string chr = seq->seq.s;
            cerr << "Chromo " << chr.size() << "bp " << K << endl;
            // get all kmers from this sequence
    		for (genomic_coordinate_t i = 0; i < chr.size() - K + 1; i++) {
                // TODO: use a rolling binary transform here to get the next kmer
    			kmer_t kmer = nimble::mer_string_to_binary(&chr[i], K);
    			// if (kmer_counts->find(kmer) == kmer_counts->end() ) {
                //     (*kmer_counts)[kmer] = 1;
    			// }
    			// else {
    			// 	if ( (*kmer_counts)[kmer] < 255)
    			// 		(*kmer_counts)[kmer]++;
                if (cbf->lookup(kmer) < 2^counting_bits) {
                    cbf->add(kmer);
    			// else -- do not increment to avoid overflow
    			}
    			if (i % 1000000 == 0) cerr << i/1000000 << "Mbp ";
    		}
            // cerr << "(" << kmer_counts->size() << " kmers) " << endl;
        }
        cerr << endl;
        // return kmer_counts;
        return cbf;
    }

    shared_ptr<unordered_map<kmer_t, list<seed_position_t> > > select_anchors(
        const string & fofn_path, const uint K, uint freq_cutoff,
        // shared_ptr<unordered_map<kmer_t, uint8_t>> kmer_counts) {
        shared_ptr<bf::counting_bloom_filter> kmer_counts,
        const int offset = 50) {
        FOFNReader fofn_reader(fofn_path);
        shared_ptr<unordered_map<kmer_t, list<seed_position_t> > > anchors(
            new unordered_map<kmer_t, list<seed_position_t> >());

        unordered_map<kmer_t,uint32_t> actual_counts;

        reference_id_t ref_id = 0;
        kseq_t * seq;
        while ( (seq = fofn_reader.getNextSequence() ) ) {
            string chr = seq->seq.s;
            for (genomic_coordinate_t i = 0; i < chr.size() - K + 1; i++) {
                // TODO: use a rolling approach to get next kmer
    			kmer_t kmer = nimble::mer_string_to_binary(&chr[i], K);
                increment_map(actual_counts, kmer);
    			if ( kmer_counts->lookup(kmer) < freq_cutoff) {
                    add_to_map(anchors, kmer, ref_id, i);
    				i += offset;
    			}
                else {
                    cerr << "kmer count above " << freq_cutoff << " true count: " <<
                        actual_counts[kmer] << endl;
                }
    			if (i % 1000000 == 0) cerr << i/1000000 << "Mbp ";
    		}
            // increment reference sequence index
            ref_id++;
            cerr << "(" << anchors->size() << " anchors from " << ref_id << " reference sequences)" << endl;
        }
        return anchors;
    }

    /*
     * build index in 2 passes
     */
	void buildIndex(const string & fofn_path, const string & output_prefix,
        const uint K) {
        cerr << "Pass 1: building background pdBG" << endl;
        // shared_ptr<unordered_map<kmer_t, uint8_t>> kmer_counts =
        shared_ptr<bf::counting_bloom_filter> kmer_counts =
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
		nimble::BloomReferenceIndex::write_index(kmer_counts, K, output_prefix);
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
