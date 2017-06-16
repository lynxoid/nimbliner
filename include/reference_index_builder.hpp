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

#include "kmer_stream.hpp"
#include "FastaReader.h"
#include "anchor_index.hpp"
#include "bloom_reference_index.hpp"
#include "kseq.h"
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

    // TODO document
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
    shared_ptr<unordered_map<kmer_t,uint8_t>> build_pdBG(const string & fofn_path, const uint K) {
      FOFNReader fofn_reader(fofn_path);
      // we only care about very rare kmers, so anything above 5 will be
      // considered frequent -> ceil(log_2 (5)) = 3 -- only need 3 bits for this
      shared_ptr<unordered_map<kmer_t,uint8_t>> map(new unordered_map<kmer_t,uint8_t>());
      kseq_t * seq;
      while ( (seq = fofn_reader.getNextSequence() ) ) {
        cerr << "Chromo " << seq->name.s << ", " << seq->seq.l << "bp " << K << endl;
        KmerStream kmer_stream(seq->seq.s, seq->seq.l, K);
        for (genomic_coordinate_t i = 0; i < seq->seq.l - K + 1; i++) {
          kmer_t kmer = kmer_stream.getNextBinKmer();
          if (map->find(kmer) == map->end() ) {
            (*map)[kmer] = 1;
          }
          else {
            if ( (*map)[kmer] < 255 )
              (*map)[kmer]++;
          }
          if (i % 1000000 == 0) cerr << i/1000000 << "Mbp, " << map->size() << " kmers | ";
        }
      }
      return map;
    }

    // keep kmers at 0, offset, 2offset, 3offset positions if their frequency is
    // below freq_cutoff; otherwise skip that kmer_t
    // TODO: may end up with long runs w/o kmers -- confirm that they are the
    // low complexity regions
    shared_ptr<unordered_map<kmer_t, list<seed_position_t> > > select_anchors(
      const string & fofn_path, const uint K, const uint freq_cutoff,
      shared_ptr<unordered_map<kmer_t,uint8_t>> kmer_counts,
      // shared_ptr<bf::counting_bloom_filter> kmer_counts,
      const int offset = 30) {
      FOFNReader fofn_reader(fofn_path);
      shared_ptr<unordered_map<kmer_t, list<seed_position_t> > > anchors(
        new unordered_map<kmer_t, list<seed_position_t> >());

      reference_id_t ref_id = 0;
      kseq_t * seq;
      while ( (seq = fofn_reader.getNextSequence() ) ) {
        cerr << seq->name.s << " ";
        KmerStream kmer_stream(seq->seq.s, seq->seq.l, K);
        for (genomic_coordinate_t i = 0; i < seq->seq.l - K + 1; i++) {
          kmer_t kmer = kmer_stream.getNextBinKmer();
          if ( (*kmer_counts)[kmer] < freq_cutoff) {
            add_to_map(anchors, kmer, ref_id, i);
            i += offset;
          }
          // else: try next available kmer
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
    shared_ptr<unordered_map<kmer_t,uint8_t>> kmer_counts = build_pdBG(fofn_path, K);

    cerr << "Pass 2: selecting anchors" << endl;
    // TODO: expose as a paramter
    const int x = 5;
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
