//
// Aligner -- the main class, i guess...
//
// dasha.filippova@gmail.com
//
// 2016

#ifndef ALIGNER
#define ALIGNER

// #include "bloom_reference_index.hpp"
#include "reference_index.hpp"
#include "definitions.hpp"

void print_read(const kseq_t * seq) {
    cerr << "read " << seq->seq.s << endl;
}

class Aligner {

  shared_ptr<ReferenceIndex> _index;
  // shared_ptr<BloomReferenceIndex> _index;

    /*
     */
  vector<genomic_coordinate_t> findAllMatchingAnchorPositions(const vector<pair<kmer_t,int>> & matched_stars,
    const shared_ptr<ReferenceIndex> index) {
    assert(matched_stars.size() > 0);

        // cerr << "(";
        // for (auto & star : matched_stars)
        //     cerr << star.second << " ";
        // cerr << ") ";

    vector<genomic_coordinate_t> mappings;
    auto first_star = matched_stars[0];
    auto second_star = matched_stars[1];
        // TODO: why +1 here? are matched_stars locations incorrect?
    int delta = second_star.second - first_star.second;
        // cerr << "delta=" << delta << " ";
        // if (delta < 0)

        // locations are sorted in increasing order
    auto A = index->get_anchor_locations(first_star.first);
        // cerr << "As ";
        // for (auto loc : A) cerr << loc << " ";

    auto B = index->get_anchor_locations(second_star.first);
        // cerr << "Bs ";
        // for (auto loc : B) cerr << loc << " ";
    int i = 0, j = 0;
    while (i < A.size() && j < B.size() ) {
      if (A[i] < B[j]) {
        if (B[j] - A[i] == delta) {
          // found a match
                    if (A[i] < first_star.second) {
                        cerr << "ERR:" << A[i] << " " << first_star.second << " " << second_star.second << endl;
                    }
          mappings.push_back(A[i] - first_star.second);
          i++;
          j++;
        }
        else if (B[j] - A[i] < delta) {
          j++;
        }
        else {
          i++;
        }
      }
      else {
        j++;
      }
    }
        // cerr << " | " << mappings.size() << " ";
    return mappings;
  }

  ////////////////////////////////////////////////////////
  //  matched_stars -- in order in which they appear in the read
  //  resolve which of the star co-locations are the same distance
  //  apart as the stars in the read
  //  0123456789
  //  --***-***-
  //
  ////////////////////////////////////////////////////////
  vector<genomic_coordinate_t> resolve_mapping_locations(const vector<pair<kmer_t,int>> & matched_stars,
    // unordered_map<kmer_t, vector<int>> & star_locations,
    int & extend,
    const int K) {
    vector<genomic_coordinate_t> mappings;

    // what if no stars mapped?
    if (matched_stars.size() == 0) {
      // need to extend the read until we hit some star
      extend++;
      // cerr << endl;
    }
    if (matched_stars.size() == 1) {
      // that's the only thing we got going
      // TODO: extend to find more stars
      auto reference_locations = _index->get_anchor_locations(matched_stars[0].first);
      for (auto loc : reference_locations) {
                mappings.push_back(loc - matched_stars[0].second);
                // cerr << mappings.back() << endl;
            }
    }
    else if (matched_stars.size() >= 2) {
      // if two stars  or more
      // TODO: take into account all the stars
      mappings = findAllMatchingAnchorPositions(matched_stars, _index);
            // for (auto p : mappings) cerr << p << " ";
            // cerr << endl;
    }

    return mappings;
  }

  /* return a 64-bit binary value that has 1nes for the 2K rightmost bits */
  kmer_t get_all_ones(const int K) {
      return ( ((kmer_t)1) << ( 2 * K ) ) - (kmer_t)1;
  }

  /*
   *
   */
  kmer_t get_next_kmer(kmer_t bin_kmer, const char next_base, const short K) {
        // mask is all ones
    kmer_t mask = get_all_ones(K - 1);
        bin_kmer = (bin_kmer & mask) << 2;
        // append a new char on the right
        bin_kmer = bin_kmer | (kmer_t)nimble::dna_codes[ next_base ];
        return bin_kmer;
  }

  /*
   * // cerr << "considering " << mer_binary_to_string(bin_kmer, K) << " ";
      // should deal w/ ends separately -- may need to trim or skip low-quality ends
      // TODO: clipping
      // TODO: check for mismatches in the first kmer (may affect the next K-1 kmers if we don't check)
   */
  kmer_t handle_first_kmer(kmer_t bin_kmer, vector<bool> & matched_kmers,
      vector<pair<kmer_t, int>> & matched_stars) {
      if (_index->has_kmer(bin_kmer) ) {
          matched_kmers.push_back(1);
      }
      else {
          matched_kmers.push_back(0);
      }
      if ( _index->is_anchor(bin_kmer) ) {
          matched_stars.emplace_back(bin_kmer, 0);
      }
      return bin_kmer;
  }

    /* a kmer in the current read is not part of the reference as it is -- we will try permuting
    the last base to see if an altered kmer is present in the reference. if it is, we will return
    the base that worked through reference_base and correct the kmer and return it through bin_kmer */
    bool is_mismatch(char * seq, const int pos, kmer_t bin_kmer, const shared_ptr<ReferenceIndex> _index,
        char & reference_base, const int K, const vector<bool> & matches, const int L) {
        if ( !matches.back() ) return false;
        // for K=4 mask is 11111100
        kmer_t mask = get_all_ones(K - 1) << 2;
        for (kmer_t i = 0; i < 4; i++) {// 0 - A, 1 - C, 2 - G, 3 - T
            // mask the last character
            bin_kmer = (bin_kmer & mask) | i;
            // cerr << "try " << mer_binary_to_string(bin_kmer, K) << " ";
            if (_index->has_kmer(bin_kmer) ) {
                // cerr << "YES ";
                // is there a branch in de Buiijn graph that allows for the following kmer as well?
                if (i + 1 >= L) {
                    // end of read -- can not verify whether this was a mismatch
                    // TODO: ignore the last base? (soft clip)
                    return false;
                }
                else {
                    kmer_t next_kmer = bin_kmer;
                    next_kmer = get_next_kmer(next_kmer, seq[ pos + K + 1], K);
                    // cerr << "next kmer " <<  mer_binary_to_string(next_kmer, K) << " ";
                    if (_index->has_kmer(next_kmer)) {
                        // cerr << "also YES ";
                        // found a base that works
                        reference_base = nimble::bases_to_bits[i];
                        return true;
                    }
                    else {
                        // cerr << "next !in_ref ";
                    }
                }
            }
            else {
                // cerr << "- NO, ";
            }
        }
        return false;
    }

    /*
     *
     */
    vector<pair<kmer_t, int>> find_anchors(const kseq_t * seq, const short K, bool debug = true) {
        if (debug) {
            cerr << "============================================================" << endl;
            print_read(seq);
        }

        vector<pair<kmer_t, int>> matched_stars;
        bool has_correction = false;
        int L = seq->seq.l - K + 1;
        assert(L < 65536); // only allow reads shorter than 2^16. WHY? aug 13 2016
        vector<bool> matched_kmers;

        kmer_t bin_kmer = nimble::mer_string_to_binary(seq->seq.s, K), bin_kmer_next;

        bin_kmer = handle_first_kmer(bin_kmer, matched_kmers, matched_stars);

        int i = 0;
        // now go through the rest of the read
        while (i + 1 < L) {
            // update prev kmer, mask the leftmost character
            bin_kmer = get_next_kmer(bin_kmer, seq->seq.s[i + K], K);

            // if (debug)
            // cerr << "i=" << i << " 0-coord: " << i+1 << " " << mer_binary_to_string(bin_kmer, K) << " ";

            // check if kmer is present in the reference -- if not, try to correct it assuming a
            // mismatch first, indels second
            if (_index->has_kmer(bin_kmer) ) {
                // TODO: could be a FP -- reduce FP rate w/ cascading BF?
                matched_kmers.push_back(1);
            }
            else {
                // if (debug)
                // cerr << "i=" << i << " 0-coord: " << i+1 << " " << mer_binary_to_string(bin_kmer, K) << " ";

                // mismatch? indel?
                // TODO: step over the base if not
                char reference_base;
                // corrects bin_kmer, returns the reference_base that worked
                bool mismatch = is_mismatch(seq->seq.s, i, bin_kmer, _index, reference_base, K, matched_kmers, seq->seq.l);

                if (!mismatch) {
                    if (debug)
                    // cerr << "could not resolve as MM ";
                    matched_kmers.push_back(0);
                    // TODO: is an indel?
                    // bool is_indel = is_indel(bin_kmer, _index);
                }
                else {
                    has_correction = true;
                    // fix the base, change the kmer, and move on
                    // seq->seq.s[i + K] = reference_base;
                    matched_kmers.push_back(1);
                    // correct the current kmer
                    kmer_t mask = get_all_ones(K-1) << 2;
                    bin_kmer = (bin_kmer & mask ) | (kmer_t)nimble::dna_codes[reference_base];
                    // if (debug)
                    // cerr << "\tcorrected kmer " << mer_binary_to_string(bin_kmer, K) << " ";
                }
                // cerr << endl;
            }

            if ( _index->is_anchor(bin_kmer) &&  matched_kmers.back() == 1 &&
                matched_kmers[matched_kmers.size() - 2] == 1 ) {
                // using i+1 since that makes it a 1-based coord system relative to read start
                matched_stars.emplace_back(bin_kmer, i + 1);
            }

            i++;
        }

        // build_cigar_string(matched_kmers);

        // print the matching pattern
        if (debug) {
            cerr << " ";
            for (auto b : matched_kmers) cerr << b;
            cerr << endl;
        }

        // if (has_correction) exit(1);

        return matched_stars;
    }

    /*
     * Find potential seeds for this read, extend the alignments, and produce
     * cigar strings, etc
     */
    vector<genomic_coordinate_t> align_single_read(const kseq_t * seq, const int K, const bool DEBUG) {
        vector<pair<kmer_t, int>> matched_stars = find_anchors(seq, K, DEBUG);
        // DEBUG info
        if (DEBUG) {
              cerr << seq->name.s << "\tmatched stars: " << matched_stars.size() << " ";
              // print anchor positions
              for (auto& p : matched_stars) {
                  cerr << p.second << ",";
              }
              cerr << " ";
        }

        // resolve star kmers to get an exact mapping location
	int need_to_extend_read = 0;
        auto mapping_locations = resolve_mapping_locations(matched_stars, need_to_extend_read, K);
        // output all potential locations for this read
        // TODO: generate CIGAR strings and all
        
        cout << seq->name.s << "\t";
        for (const auto & loc : mapping_locations) {
          cout << loc << " ";
        }
        cout << endl;
        
        if (DEBUG) {
            cerr << "mapping locations: " << mapping_locations.size() << " -- ";
            for (const auto & loc : mapping_locations) {
                cout << loc << " ";
            }
            cerr << endl;
        }
        return mapping_locations;
    }

public:

    Aligner(const shared_ptr<ReferenceIndex> index) {
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
    while ( (seq = fr.nextSequence() ) ) {
      read_count++;
      if (read_count % 100000 == 0)
        cerr << read_count / 100000 << "00K ";
      if (seq->seq.l < K) {
        continue;
      }
      auto start = std::chrono::system_clock::now();
      align_single_read(seq, K, true);
      auto end = std::chrono::system_clock::now();
      elapsed_seconds_str += end - start;
    }
    cerr << "Computing alignments (no IO): " << elapsed_seconds_str.count() << "s" << endl;
    // cerr << "Times needed to extend the read: " << need_to_extend_read << endl;
  }
};

#endif
