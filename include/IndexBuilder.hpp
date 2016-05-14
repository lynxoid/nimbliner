//
//
// Indexing the reference
//
//

#ifndef INDEX_BUILDER
#define INDEX_BUILDER

class IndexBuilder {

public:
	////////////////////////////////////////////////////////
	// build index stage
	////////////////////////////////////////////////////////
	void getAllKmersAndStars(const string & ref_path, int K) {
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
	}
};

#endif