//
// IndexReader
//

#ifndef INDEX_READER
#define INDEX_READER

class IndexReader {
public:
	////////////////////////////////////////////////////////
	BaseBloomFilter readKmers_bf_faster(const string & path, int K) {
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
		BaseBloomFilter bf(K, kmer_count * 10);

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
		        		bf.add(bin_kmer);
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

		return bf;
	}

	////////////////////////////////////////////////////////
	unordered_map<kmer_t, vector<int>> readStarLocations(const string & path, int K) {
		cerr << "reading stars from " << path << endl;
		auto start = std::chrono::system_clock::now();

		unordered_map<kmer_t, vector<int>> kmer_locations;
		ifstream in(path);
		string line;
		while (getline(in, line)) {
			// parse the line now
			std::stringstream ss(line);
			string kmer, pos;
			getline(ss, kmer, ' ');
			uint64_t bin_kmer = stol(kmer);
			kmer_locations.emplace(bin_kmer, vector<int>{});
			while(getline(ss, pos, ' ')) {
	    		int location = stoi(pos);
	    		kmer_locations[bin_kmer].push_back(location);
			}
			if (kmer_locations[bin_kmer].size() > 10000) {
				cerr << mer_binary_to_string(bin_kmer, K) << " " << kmer_locations[bin_kmer].size() << endl;
				kmer_locations.erase( kmer_locations.find(bin_kmer) );
			}
		}
		in.close();

		auto end = std::chrono::system_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
		cerr << "reading stars took: " << elapsed_seconds.count() << "s" << endl;

		return kmer_locations;
	}
};

#endif