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

class ReferenceIndex {

	// probabilistic dBG -- holds the dBG of the reference
	shared_ptr<BaseBloomFilter> _bf;

	// anchor locations (stars)
	// TODO: how to choose these intelligently
	shared_ptr<unordered_map<kmer_t, vector<int>>> _stars;

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
	shared_ptr<unordered_map<kmer_t, vector<int>>> readStarLocations(const string & path, int K) {
		cerr << "reading stars from " << path << endl;
		auto start = std::chrono::system_clock::now();

		shared_ptr<unordered_map<kmer_t, vector<int>>> kmer_locations = 
			shared_ptr<unordered_map<kmer_t, vector<int>>>(new unordered_map<kmer_t, vector<int>>());
		ifstream in(path);
		string line;
		while (getline(in, line)) {
			// parse the line now
			std::stringstream ss(line);
			string kmer, pos;
			getline(ss, kmer, ' ');
			uint64_t bin_kmer = stol(kmer);
			kmer_locations->emplace(bin_kmer, vector<int>{});
			while(getline(ss, pos, ' ')) {
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
		cerr << "reading stars took: " << elapsed_seconds.count() << "s" << endl;

		return kmer_locations;
	}

public:

	ReferenceIndex() {};

	// read kmers (or pdBG describing the reference) and anchor locations
	void readIndex(const string & kmers_path, const string & stars_path, const uint K) {
		_bf = readKmers_bf_faster(kmers_path, K);
		_stars = readStarLocations(stars_path, K);
	}

	// given a path to the reference and a kmer length K, parse the reference
	// fasta, obtain all kmers and select anchors using a heuristic
	// write kmers to an *.index file and write anchors to *.star file.
	void buildIndex(const string & ref_path, int K) {
		auto chromosomes = parseFasta(ref_path);
		// count kmers, record their locations
		unordered_map<kmer_t, vector<int>> kmer_locations;

		unordered_set<kmer_t> star_kmers;
		assert(chromosomes.size() > 0);
		auto chr = chromosomes[0];

		int c = 0;
		for (int i = 0; i < chr.size() - K + 1; i++) {
			kmer_t kmer = mer_string_to_binary(&chr[i], K);
			if ( kmer_locations.find(kmer) == kmer_locations.end() ) {
				kmer_locations.emplace(kmer, vector<int>{i});
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
		cerr << "(" << star_kmers.size() << " star kmers)" << endl;

		ofstream star_locations_out("star_locations.txt");
		ofstream all_kmers("all_kmers.txt");
		all_kmers << kmer_locations.size() << endl;
		cerr << "saving star locations" << endl;
		for (auto star : star_kmers) {
			star_locations_out << star << " ";
			for (auto loc : kmer_locations[star]) star_locations_out << loc << " ";
			star_locations_out << endl;
			all_kmers << star << endl;
			kmer_locations.erase(star);
		}
		cerr << "saving all kmers" << endl;
		for (auto p : kmer_locations) {
			all_kmers << p.first << endl;
		}
		all_kmers.close();
		star_locations_out.close();
		kmer_locations.clear();
		return;
	}

	bool has_anchor(const kmer_t & kmer) const {
		return _stars->find(kmer) != _stars->end();
	}

	// TODO: what does this & do?
	vector<int> & get_anchor_locations(const kmer_t & kmer) const {
		return (*_stars)[kmer];
	}
};

#endif