#include <boost/progress.hpp>
#include "bit_tree_index.hpp"
#include "reference_index.hpp"

using namespace std;

int main(int argc, char * argv []) {
	shared_ptr<BitTreeIndex> index = shared_ptr<BitTreeIndex>(new BitTreeIndex() );
	string kmers_path = "data/output_mismatches/chr20/index/chr20.index.btbin";
	string stars_path = "data/output_mismatches/chr20/index/chr20.star";
	index->readIndex(kmers_path, stars_path, 20);

	shared_ptr<vector<kmer_type>> kmers = index->get_kmers();

	boost::timer t;
	int n = 100000;
	for (int i = 0; i < n; i++) {
		int j = kmers->size() / n * i;
		bin_kmer_t bin_kmer = (*kmers)[j];
		
		assert( index->has_kmer(bin_kmer) );	
	}
	cerr << "Testing for " << n << " kmers took " << t.elapsed() << "s" << endl;
}