// parallel aligner
//
// dasha.filippova@gmail.com
//
// 2016

#ifndef PARALLEL_ALIGNER
#define PARALLEL_ALIGNER

#include <iostream>
#include <chrono>
#include <vector>
#include <memory>

#include <tbb/concurrent_queue.h>

#include "reference_index.hpp"
#include "definitions.hpp"
#include "FastaReader.h"

class ParallelAligner {

  shared_ptr<ReferenceIndex> _index;

  void align_read(const kseq_t * seq, const int K, const bool debug = false) {
    // vector<pair<kmer_t, int>> matched_stars = find_anchors(seq, K, debug);
    // // DEBUG info
    // if (debug) {
    //   cerr << seq->name.s << "\tmatched stars: " << matched_stars.size() << " ";
    //   // print anchor positions
    //   for (auto& p : matched_stars) {
    //       cerr << p.second << ",";
    //   }
    //   cerr << " ";
    // }
    //
    // // resolve star kmers to get an exact mapping location
    // auto start = std::chrono::system_clock::now();
    // auto mapping_locations = resolve_mapping_locations(matched_stars, need_to_extend_read, K);
    // auto end = std::chrono::system_clock::now();
    // if (debug) cerr << endl;
    //
    // // output all potential locations for this read
    // // TODO: generate CIGAR strings and all
    // cout << seq->name.s << "\t";
    // for (const auto & loc : mapping_locations) {
    //   cout << loc << " ";
    // }
    // cout << endl;
  }

public:

  ParallelAligner(const shared_ptr<ReferenceIndex> index) {
    _index = index;
  }

  ////////////////////////////////////////////////////////
  //
  ////////////////////////////////////////////////////////
  void alignReads(const string & path, const int K, bool debug = false) {
    FastaReader fr(path.c_str());
    kseq_t * seq;
    int passed_cutoff = 0;
    int read_count = 0;
    int need_to_extend_read = 0;

    std::chrono::duration<double> elapsed_seconds_str;
    extern concurrent_queue<kseq_t*> reads_queue;  
    while ( (seq = fr.nextSequence() ) ) {
      reads_queue.push_back();

      read_count++;
      if (read_count % 100000 == 0)
        cerr << read_count / 100000 << "00K ";
      if (seq->seq.l < K) {
        continue;
      }

      // auto start = std::chrono::system_clock::now();
      align_read(seq, K, debug);

      // elapsed_seconds_str += end - start;
    }
    // cerr << "Resolving ops: " << elapsed_seconds_str.count() << "s" << endl;
    // cerr << "Times needed to extend the read: " << need_to_extend_read << endl;
  }
};

#endif // PARALLEL_ALIGNER
