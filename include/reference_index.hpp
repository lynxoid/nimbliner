//
// IndexReader
//

#ifndef INDEX_READER
#define INDEX_READER

#include <unordered_map>
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <memory>

#include "BaseBloomFilter.hpp"
#include "definitions.hpp"

using namespace std;


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

class ReferenceIndex {

	// probabilistic dBG -- holds the dBG of the reference
	shared_ptr<BaseBloomFilter> _bf;

	// anchor locations (stars)
	shared_ptr<unordered_map<kmer_t, vector<genomic_coordinate_t>>> _stars;

	// read kmers describing the reference from the file
	shared_ptr<BaseBloomFilter> readKmers_bf_faster(const string & path, int K) {
		cerr << "reading kmers from " << path << endl;
		auto start = std::chrono::system_clock::now();

		ifstream in(path);
		if (!in) {
			cerr << "[ERROR] Could not open the file: " << path << endl;
			exit(1);
		}
		string line;
		getline(in, line);
		uint64_t kmer_count = stol(line);
		// cerr << "Expected kmer count: " << kmer_count << endl;
		// BaseBloomFilter bf(K, kmer_count * 10);
		shared_ptr<BaseBloomFilter> bloom = shared_ptr<BaseBloomFilter>(new BaseBloomFilter(K, kmer_count * 10) );

		static const auto BUFFER_SIZE = (size_t)pow(2,22); // do not overwhelm the stack :)
		// cerr << "Buffer size, bytes: " << BUFFER_SIZE << endl;
	    int fd = open(path.c_str(), O_RDONLY);
	    if (fd == -1) {
	    	cerr << "[ERROR] Could not open file: " << path << endl;
	    	exit(1);
	    }
	 //    /* Advise the kernel of our access pattern.  */
	 //    // posix_fadvise(fd, 0, 0, 1);  // FDADVICE_SEQUENTIAL
	    char buf[BUFFER_SIZE + 1];
	    size_t i = 0, trailing = 0;
	    size_t bytes_read = read(fd, buf, BUFFER_SIZE);
	    size_t total_bytes_read = bytes_read;
	    while (bytes_read != 0) {
	        if (bytes_read == (size_t)-1) {
	        	cerr << "Could not read" << endl;
	        	exit(1);
	        }
	        char *p = buf;
	        trailing = 0;
	        while ( (p - buf) < bytes_read) {
	        	char * p_next = (char*) memchr(p, '\n', (buf + bytes_read) - p);
	        	if (p_next == NULL) {

	        		trailing = buf + bytes_read - p;
	        		for (int j = 0; j < trailing; j++)
	        			buf[j] = *(p + j);
	        		break;
	        	}
	        	else {
		        	if (i != 0) {
		        		uint64_t bin_kmer = strtol(p, &p_next, 10);
		        		// if (bytes_read < BUFFER_SIZE/2)
		        			// cerr << (p - buf) << " " << bin_kmer << endl;
		        		bloom->add(bin_kmer);
		        	}
		        	p = p_next + 1;

		        	if (i % 5000000 == 0) cerr << i/1000000 << "m ";
		        	if (i > 1 + kmer_count) {
		        		cerr << "More kmers than expected: " << kmer_count << " vs. " << i << " bytes read: " << total_bytes_read << endl;
		        		cerr << (p-buf) << " " << bytes_read << endl;
		        		exit(1);
		        	}
		        	i++;
		        }
	        }
	        bytes_read = read(fd, buf + trailing, BUFFER_SIZE - trailing);
	        // cerr << bytes_read << "b ";
	        total_bytes_read += bytes_read;
	    }
		close(fd);
		cerr << endl << "read " << (i-1) << " kmers" << endl;
		// cerr << "Trailing: " << trailing << " " << bytes_read << endl;
		if (i < kmer_count) {
			cerr << "Expected: " << kmer_count << endl;
			// exit(1);
		}
		// debug info
		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		cerr << "reading kmers took: " << elapsed_seconds.count() << "s" << endl;
		return bloom;
	}

	// read kmers that serve as anchors
	shared_ptr<unordered_map<kmer_t, vector<genomic_coordinate_t>>> readStarLocations(const string & path, const int K) {
		cerr << "reading stars from " << path << endl;
		auto start = std::chrono::system_clock::now();

		shared_ptr<unordered_map<kmer_t, vector<genomic_coordinate_t>>> kmer_locations = 
			shared_ptr<unordered_map<kmer_t, vector<genomic_coordinate_t>>>(new unordered_map<kmer_t, vector<genomic_coordinate_t>>());
		ifstream in(path);
		string line;
		while (getline(in, line)) {
			// parse the line now
			std::stringstream ss(line);
			string kmer, pos;
			getline(ss, kmer, ' ');
			uint64_t bin_kmer = stol(kmer);
			kmer_locations->emplace(bin_kmer, vector<genomic_coordinate_t>{});
			while(getline(ss, pos, ' ')) {
	    		// int location = stoul(pos);
	    		int location = stoi(pos);
	    		(*kmer_locations)[bin_kmer].push_back(location);
			}
			if ((*kmer_locations)[bin_kmer].size() > 10000) {
				cerr << mer_binary_to_string(bin_kmer, K) << " " << (*kmer_locations)[bin_kmer].size() << endl;
				kmer_locations->erase( kmer_locations->find(bin_kmer) );
			}
		}
		in.close();

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		cerr << "reading anchors took: " << elapsed_seconds.count() << "s" << endl;

		return kmer_locations;
	}

public:

	ReferenceIndex() {};

	// read kmers (or pdBG describing the reference) and anchor locations
	void readIndex(const string & kmers_path, const string & stars_path, const uint K) {
		_bf = readKmers_bf_faster(kmers_path, K);
		_stars = readStarLocations(stars_path, K);
	}

	/* returns true if this kmer was present in the reference sequence, false otherwise */
	bool has_kmer(const kmer_t & kmer) const {
		return _bf->contains(kmer);
	}

	/* returns true is this kmer is found among anchors, false otherwise*/
	bool has_anchor(const kmer_t & kmer) const {
		return _stars->find(kmer) != _stars->end();
	}

	/* allow anchors w/ 1 mm  */
	// bool has_anchor(const kmer_t & kmer, const unsigned char K = 20) const {
	// 	// TODO: generate all variants of this kmer
	// 	// and try them all
	// 	// TODO: SIMD optimization
	// 	shared_ptr<vector<kmer_t>> kmer_variants = generate_all_variants(kmer, K, 1);
	// 	for (const kmer_t & kmer : *kmer_variants) {
	// 		if (_stars->find(kmer) != _stars->end()) return true;
	// 	}
	// 	return false;
	// }

	// TODO: what does this & do? do we use it?
	vector<genomic_coordinate_t> & get_anchor_locations(const kmer_t & kmer) const {
		return (*_stars)[kmer];
	}
};

#endif