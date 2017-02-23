// Binary BitTree library
#ifndef BIT_TREE_BINARY_LIB
#define BIT_TREE_BINARY_LIB

#include <vector>
#include <stack>
#include <unordered_map>
#include <fstream>

#include <boost/dynamic_bitset.hpp>
#include <boost/progress.hpp>

#include "JellyfishUtil.h"

using namespace std;

#define bin_kmer_t uint64_t

#define kmer_type bin_kmer_t

size_t getFileLength(ifstream & stream) {
  // get file length
  stream.seekg(0, stream.end);  // move to the end
  int length = stream.tellg();   // get the size in bytes?
  stream.seekg(0, stream.beg);  // move back to the beginning
  return length;
}

class BitTreeBin {
private:
	// contains binary data about the BitTree
	shared_ptr<boost::dynamic_bitset<>> bitstream;

	unordered_map<bin_kmer_t, size_t> dmer_index;

	uint D = 0;

	// kmer size
	int k = 0;

	//////////////////////////////////////////////////////////////////////////////////////////
	// given a range (start, end), find a subrange that has char c at position depth
	//////////////////////////////////////////////////////////////////////////////////////////
	pair<int,int> findRange(int start, int end, shared_ptr<vector<kmer_type>> kmers, int depth, int k, char c) {
		auto match = find(table, table + 4, c);
		int letter = match - table;
		assert(letter >= 0);
		// cerr << "letter " << c << "'s index " << letter << " ";

		int i = start;
		int shift = (k - depth - 1) * 2;
		// string prefix = mer_binary_to_string( (*kmers)[i] >> shift, depth+1);
		while ( ( ( (*kmers)[i] >> shift ) & 3) != letter && i < end) i++;
		// prefix = mer_binary_to_string( (*kmers)[i] >> shift, depth+1);
		int a = i;
		i = a;
		while ( ( ( (*kmers)[i] >> shift ) & 3) == letter && i < end) i++;

		int b = i;
		return pair<int,int>(a, b);
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	// given a range (start, end), find a subrange that has char c at position depth
	//////////////////////////////////////////////////////////////////////////////////////////
	void write64(ofstream & f_output, unsigned long long n) {
		f_output << (unsigned char) ((n >> 0) & (uint64_t)255);
		f_output << (unsigned char) ((n >> 8) & (uint64_t)255);
		f_output << (unsigned char) ((n >> 16) & (uint64_t)255);
		f_output << (unsigned char) ((n >> 24) & (uint64_t)255);
		f_output << (unsigned char) ((n >> 32) & (uint64_t)255);
		f_output << (unsigned char) ((n >> 40) & (uint64_t)255);
		f_output << (unsigned char) ((n >> 48) & (uint64_t)255);
	    f_output << (unsigned char) ((n >> 56) & (uint64_t)255);
	}

public:

	// reference table for translation from/to binary
	/*const*/ char table[4] = { 'A', 'C', 'G', 'T' };

	//////////////////////////////////////////////////////////////////////////////////////////
	// given a set of kmers (in binary), encode the trie they form in binary
	//////////////////////////////////////////////////////////////////////////////////////////
	void encode(shared_ptr<vector<kmer_type>> kmers, int const k) {
		stack<int> starts; starts.push(0);
		stack<int> ends; ends.push(kmers->size());
		stack<int> depths; depths.push(0);
		string alphabet = "TGCA";

		shared_ptr<boost::dynamic_bitset<>> bitstream( new boost::dynamic_bitset<>());
		int start, end, depth;

		while (starts.size() > 0) {
			start = starts.top(); starts.pop();
			end = ends.top(); ends.pop();
			depth = depths.top(); depths.pop();
			// cerr << "(" << start << "," << end << ") ";
			// if (bitstream->size() % 1000000 == 0) { cerr << bitstream->size() << " "; }
			if (start == end) {
				bitstream->push_back(0);
			}
			else {
				bitstream->push_back(1);
				if (depth < k) {
					for (auto c : alphabet) {
						pair<int,int> p = findRange(start, end, kmers, depth, k, c);
						// cerr << depth << ", " << c << ": (" << p.first << "," << p.second << ")" << endl;
						starts.push(p.first); ends.push(p.second); depths.push(depth+1);
					}
				}
			}
		}
		// string s;
		// boost::to_string(*bitstream, s);
		// cerr << s << endl;
		kmers->clear();

		this->bitstream = bitstream;
		this->k = k;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	// given binary stream of a trie, reconstruct the strings
	//////////////////////////////////////////////////////////////////////////////////////////
	shared_ptr<vector<kmer_type>> decode() {
		shared_ptr<vector<kmer_type>> kmers(new vector<kmer_type>());

		// skip the first 1 -- it is the root
		size_t index = bitstream->find_first() + 1;
		// string alphabet = "ACGT";

		int branch = 0;
		int depth = 0;
		bin_kmer_t prefix = 0;
		// put branches index on stack as we go down the tree, pop as we go back up the tree
		// keep track of unexplored branches
		stack<int> branches;

		while (index < bitstream->size() ) {
			// cerr << (*bitstream)[index] << " ";
			if ( (*bitstream)[index] ) { // current bit is 1
				prefix = (prefix << 2) | branch;
				if (depth == (k-1) ) { // reached a leaf
					kmers->push_back(prefix); // prefix is a decoded kmer
					// output some stats
					if (kmers->size() % 1000000 == 0) cerr << kmers->size() << " ";
					prefix = prefix >> 2;

					while (branch >= 3 && !branches.empty()) {
							depth--;
							branch = branches.top(); branches.pop();
							prefix = prefix >> 2;
					}
					if (branch < 3)
						branch++;
				}
				else {
					branches.push(branch); // remember the current branch
					branch = 0; // start from branch A at the next level
					depth++;
				}
			}
			else { // current bit is 0
				while (branch >= 3 && !branches.empty()) {
					depth--;
					branch = branches.top(); branches.pop();
					prefix = prefix >> 2;
				}
				if (branch < 3) //skip to the next branch
					branch++;
			}
			index++; // move to the next bit
		}
		// cerr << endl;
		return kmers;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	// write binary stream to a file
	//////////////////////////////////////////////////////////////////////////////////////////
	void write(string const & fname) {
		cerr << "k = " << k << endl;
		ofstream f_output(fname, ios::binary | ios::out);
		f_output << (char)k; // write out the value of k
		auto size = sizeof(uint64_t)*8;
		auto tail = bitstream->size() % size;
		// pad with 0s at the end so that the tail does not get messed up when writing
		bitstream->resize(bitstream->size() + size - tail);
		*bitstream <<= (size - tail);

		// convert into a vector of block_types
		vector<boost::dynamic_bitset<>::block_type> OutIter;
		to_block_range(*bitstream, back_inserter(OutIter) );
		// iterate and write one block_type at a time
		for (boost::dynamic_bitset<>::block_type block : OutIter) {
			write64(f_output, block);
		}
		f_output.close();
		cerr << "Output written to " << fname << endl;
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	// read binary stream into the dynamic_bitset
	//////////////////////////////////////////////////////////////////////////////////////////
	void read(string const & fname, int & k) {
		ifstream f_in(fname, ios::binary | ios::in);
		auto fsize = getFileLength(f_in);
		f_in.read(reinterpret_cast<char *>(&k), 1);
		cerr << "BitTree:: k = " << (int)k << endl;
		assert(k > 0 && k < 33);
		shared_ptr<boost::dynamic_bitset<>> bitstream(new boost::dynamic_bitset<>());
		// add buffer[1:] to the bitstream
		auto bytes_per_block = sizeof(boost::dynamic_bitset<>::block_type);
		cerr << "Reading bit tree at " << bytes_per_block << "bytes per block"
			<< endl;
		auto num_blocks = (fsize - 1) / bytes_per_block;
		// cerr << num_blocks << endl;
		boost::dynamic_bitset<>::block_type block;
		for (int i = 0; i < num_blocks; i++) {
			f_in.read(reinterpret_cast<char *>(&block), bytes_per_block);
			bitstream->append(block);
		}
		this->bitstream = bitstream;
		this->k = k;
	}

	/*
	 * build an index -- an offset to the bitstream for every D-mer we observe
	 * in the kmer set. Store the index in a hash map of (5mer, index) pairs
	 */
	void build_index(const uint D) {
		this->D = D;
		dmer_index.clear();

		// skip the first 1 -- it is the root
		size_t index = bitstream->find_first() + 1;
		// string alphabet = "ACGT";

		bin_kmer_t branch = 0;
		int depth = 0;
		bin_kmer_t prefix = 0;
		// put branches index on stack as we go down the tree, pop as we go back up the tree
		// keep track of unexplored branches
		stack<int> branches;

		while (index < bitstream->size() ) {
			// cerr << (*bitstream)[index] << " ";
			if ( (*bitstream)[index] ) { // current bit is 1
				prefix = (prefix << 2) | branch;
				if (depth == (D-1) ) { // reached level D
					dmer_index.emplace(prefix, index);
					prefix = prefix >> 2;

					while (branch >= 3 && !branches.empty()) {
							depth--;
							branch = branches.top(); branches.pop();
							prefix = prefix >> 2;
					}
					if (branch < 3)
						branch++;
				}
				else {
					branches.push(branch); // remember the current branch
					branch = 0; // start from branch A at the next level
					depth++;
				}
			}
			else { // current bit is 0
				while (branch >= 3 && !branches.empty()) {
					depth--;
					branch = branches.top(); branches.pop();
					prefix = prefix >> 2;
				}
				if (branch < 3) //skip to the next branch
					branch++;
			}
			index++; // move to the next bit
		}
		cerr << "created an index into the bit tree w/ " << dmer_index.size() <<
			" elements (max elements: " << pow(4, D) << ")" << endl;
	}

	/*
	 * returns True if the bit_tree (a set) contains the kmers, otherwise
	 * returns False. Computes the answer in O(L) where L is the length of the
	 * query (kmer size)
	 */
	bool contains(const bin_kmer_t query) {
		// compute offset into the stream to query for the first D letters
		bin_kmer_t prefix = query >> ( (k - D) * 2 );
		auto it = dmer_index.find(prefix);
		if (it == dmer_index.end() ) // prefix does not exist
			return false;
		size_t index = it->second;

		////////////////////////////////////////////////////////
		// now traverse the subtree starting at index to see if kmer is there
		////////////////////////////////////////////////////////

		bin_kmer_t branch = 0;
		// depth is 0-based, while D is 1-based. here we are at depth D+1 -- one
		// longer than prefix
		int depth = D;
		// put branches index on stack as we go down the tree, pop as we go back
		// up the tree keep track of unexplored branches
		stack<int> branches;

		while (index < bitstream->size() ) {
			// cerr << (*bitstream)[index] << " ";
			if ( (*bitstream)[index] ) { // current bit is 1
				prefix = (prefix << 2) | branch; // branch is always 0 here
				if (depth == (k - 1) ) { // reached level K
					if (prefix == query)
						return true;
					// else: reached a leaf, but did not leave the subtree that
					// may have the query
					// continue exploring the subtrees
					prefix = prefix >> 2; // up one level

					while (branch >= 3 && !branches.empty()) {
							depth--;
							// if (depth < D) return false;
							branch = branches.top(); branches.pop();
							prefix = prefix >> 2;
					}
					if (branch < 3)
						branch++;
					// did we go past the query already?
					bin_kmer_t query_i = (query >> ( (k - depth) * 2) ) & (bin_kmer_t)3;
					if (branch > query_i )
						return false;
				}
				else {
					branches.push(branch); // remember the current branch
					branch = 0; // start from branch A at the next level
					depth++;
				}
			}
			else { // current bit is 0
				while (branch >= 3 && !branches.empty()) {
					depth--;
					// if (depth < D) return false;
					branch = branches.top(); branches.pop();
					prefix = prefix >> 2;
				}
				if (branch < 3) //skip to the next branch
					branch++;
				// check if we already went past our query
				// if symbol at query[i] < branch -- symbol is not in the tree
				// compute symbol at query[i] in its binary form
				bin_kmer_t query_i = (query >> ( (k - depth) * 2) ) & (bin_kmer_t)3;
				if (branch > query_i )
					return false;
			}
			index++; // move to the next bit
		}
		return false;
	}

};

#endif // BIT_TREE_BINARY_LIB
