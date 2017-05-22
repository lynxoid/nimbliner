//
// IndexReader
//

#ifndef BLOOM_INDEX_READER
#define BLOOM_INDEX_READER

#include <unordered_map>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <memory>
#include <cassert>

// for low-level file IO
#include <fcntl.h>
#include <unistd.h>
#include <string.h>

#include "BaseBloomFilter.hpp"
#include "reference_index.hpp"
#include "anchor_index.hpp"
#include "definitions.hpp"

using namespace std;

namespace nimble {

/* generate all possible variants for a kmer of a given length w/ hamming
distance mm from the original */
shared_ptr<vector<kmer_t>> generate_all_variants(const kmer_t & kmer, const unsigned char K, const unsigned char mm = 1) {
	shared_ptr<vector<kmer_t>> variants( new vector<kmer_t>() );
	for (int i = 0; i < K; i++) {
		// create a mask w/ 11 in all positions but i, 00 in i position
		kmer_t mask = pow(2, K * 2) - 1;
		mask = mask ^ ( ( (kmer_t) 3) << (2*i) );
		for (kmer_t base = 0; base < 4; base++) {
			// apply mask to the kmer -- creates two 0s in the 2i position
			kmer_t masked_kmer = kmer & mask;
			// now OR the base at that position w/ the kmer
			variants->push_back( masked_kmer | (base << (2*i) ) );
		}
	}
	return variants;
}

class BloomReferenceIndex : public ReferenceIndex {

	// probabilistic dBG -- holds the dBG of the reference
	shared_ptr<BaseBloomFilter> _bf;

	// anchor locations (stars)
    nimble::AnchorIndex _anchorIndex;

	// kmer length
	uint K = 0;

    uint64_t observed_kmers = 0;

	/*
	 *
	 * read kmers describing the reference from the file
	 *
	 */
	shared_ptr<BaseBloomFilter> readKmers_bf_faster(const string & prefix) {
		cerr << "[BloomFilterIndex] reading index from " << prefix << EXT << endl;
		auto start = std::chrono::system_clock::now();

		// vector<kmer_bin_t> all_kmers;

		ifstream in(prefix + EXT);
		if (!in) {
			cerr << "[ERROR] [BloomFilterIndex] Could not open the file: " << prefix << endl;
			exit(1);
		}
		string line;
        // read first line: K -- kmer length
		getline(in, line);
		uint64_t k = stol(line);
		this->K = k;
		cerr << "[BloomFilterIndex] Kmer length: " << K << endl;
        // read second line: N -- number of kmers (lines) to expect
		getline(in, line);
		uint64_t kmer_count = stol(line);
		cerr << "[BloomFilterIndex] Expected kmer count: " << kmer_count << endl;

        // read the rest of the file -- one kmer per line
		shared_ptr<BaseBloomFilter> bloom =
			shared_ptr<BaseBloomFilter>(new BaseBloomFilter(K, (size_t)(kmer_count * 10) ) );

		static const auto BUFFER_SIZE = (size_t)pow(2,22); // do not overwhelm the stack :)
        // TODO why are we reading the file from the start? need to skip first 2 lines
	    int fd = open((prefix + EXT).c_str(), O_RDONLY);
	    if (fd == -1) {
	    	cerr << "[ERROR] [BloomFilterIndex] Could not open file: " << prefix << EXT << endl;
	    	exit(1);
	    }
		else {
			cerr << "Opened index file successfully" << endl;
		}
	 	//    /* Advise the kernel of our access pattern.  */
	 	//    // posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL
	    char buf[BUFFER_SIZE + 1];
        // i -- count the lines
	    size_t i = 0, trailing = 0;
	    size_t bytes_read = read(fd, buf, BUFFER_SIZE);
	    size_t total_bytes_read = bytes_read;
	    while (bytes_read != 0) {
	        if (bytes_read == (size_t)-1) {
	        	cerr << "[BloomFilterIndex] Could not read" << endl;
	        	exit(1);
	        }
	        char *p = buf;
	        trailing = 0;
	        while ( (p - buf) < bytes_read) {
	        	// find ends of lines
	        	char * p_next = (char*) memchr(p, '\n', (buf + bytes_read) - p);
	        	if (p_next == NULL) {
	        		trailing = buf + bytes_read - p;
	        		for (int j = 0; j < trailing; j++)
	        			buf[j] = *(p + j);
	        		break;
	        	}
	        	else {
		        	if (i != 0 && i != 1) { // skip the first 2 lines
		        		uint64_t bin_kmer = strtol(p, &p_next, 10);
		        		// if (bytes_read < BUFFER_SIZE/2)
		        			// cerr << (p - buf) << " " << bin_kmer << endl;
		        		// Test: time the speen w/o adding stuff to BF
		        		bloom->add(bin_kmer);

		        		// all_kmers.push_back(bin_kmer);
		        	}
		        	p = p_next + 1;

		        	if (i % 5000000 == 0) cerr << i/1000000 << "m ";
		        	if (i > 1 + kmer_count) {
		        		cerr << "[BloomFilterIndex] Read more kmers than expected: "
                            << kmer_count << " vs. " << i <<
                            " bytes read: " << total_bytes_read << endl;
		        		cerr << (p-buf) << " " << bytes_read << endl;
		        		exit(1);
		        	}
		        	i++; // increment the line count
		        }
	        }
	        bytes_read = read(fd, buf + trailing, BUFFER_SIZE - trailing);
	        // cerr << bytes_read << "b ";
	        total_bytes_read += bytes_read;
	    }
		close(fd);
        observed_kmers = i - 2;
		cerr << endl << "read " << observed_kmers << " kmers" << endl;
		cerr << "[BloomFilterIndex] Trailing: " << trailing << " " << bytes_read << endl;
		if (i < kmer_count) {
			cerr << "[BloomFilterIndex] Read fewer kmers than expected: " << i <<
                "vs. " << kmer_count << endl;
		}
		// debug info
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		cerr << "[BloomFilterIndex] reading kmers took: " << elapsed_seconds.count() << "s" << endl;

		// cerr << "sorting" << endl;
		// std::sort(all_kmers.begin(), all_kmers.end());

		return bloom;
	}

public:

	BloomReferenceIndex() {};

    const constexpr static char* EXT = ".idx";

	uint getK() {return K;}

    uint64_t size() {
        return observed_kmers;
    }

	// read kmers (or pdBG describing the reference) and anchor locations
	void readIndex(const string & index_prefix) {

		this->K = 0;
		_bf = readKmers_bf_faster(index_prefix);
        // TODO: K is not set here -- ro anywhere
        _anchorIndex.readAnchors(index_prefix, this->K);
	}

	/* returns true if this kmer was present in the reference sequence, false otherwise */
	bool has_kmer(const bin_kmer_t kmer) const {
		return _bf->contains(kmer);
	}

	/* returns true is this kmer is found among anchors, false otherwise*/
	bool is_anchor(const bin_kmer_t kmer) const {
		// return _anchors->find(kmer) != _anchors->end();
        return _anchorIndex.is_anchor(kmer);
	}

	// TODO: what does this & do? do we use it?
	vector<seed_position_t> & get_anchor_locations(const bin_kmer_t & kmer) const {
		return _anchorIndex.get_anchor_locations(kmer);
	}

    /*
     * Writes all kmers in the index in the format that is agreed upon between
     * this function and readIndex()
     */
    static void write_index(shared_ptr<unordered_map<kmer_t,uint8_t>> kmer_counts,
        const uint K,
        const string & output_prefix) {
		auto start = std::chrono::system_clock::now();
		ofstream all_kmers(output_prefix + EXT);
        // TODO: check that can write
        cerr << "[BloomFilterIndex] writing index kmers to " << output_prefix + EXT << endl;
        // write out the K -- kmer length
        all_kmers << K << endl;
		// write the # of kmers to expect
		all_kmers << kmer_counts->size() << endl;

		int i = 0;
        // write and prune
		while (kmer_counts->size() > 0) {
			auto it = kmer_counts->begin();
			all_kmers << it->first << endl;
			kmer_counts->erase(it);
			i++;
		}
		all_kmers.close();
		cerr << "[BloomFilterIndex] wrote " << i << " kmers" << endl;
		assert(kmer_counts->size() == 0);

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds_str = end - start;
	    cerr << "[BloomFilterIndex] Saving kmers took: " << elapsed_seconds_str.count() << "s" << endl;
		kmer_counts->clear();
    }

    static void write_index(shared_ptr<bf::counting_bloom_filter> kmer_counts,
        const uint K,
        const string & output_prefix) {
        cerr << "TODO" << endl;
	}
};

}

#endif // BLOOM_INDEX_READER
