// Binary BitTree library
#ifndef BIT_TREE_BINARY_LIB
#define BIT_TREE_BINARY_LIB

#include <vector>
#include <stack>

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
	// contain binary data about the BitTree
	shared_ptr<boost::dynamic_bitset<>> bitstream;

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
		
		// skip the first 1
		size_t index = bitstream->find_first() + 1;
		// string alphabet = "ACGT";
		
		int branch = 0;
		int depth = 0;
		bin_kmer_t prefix = 0;
		// put branches index on stack as we go down the tree, pop as we go back up the tree
		stack<int> branches;

		while (index < bitstream->size() ) {
			// cerr << (*bitstream)[index] << " ";
			if ( (*bitstream)[index] ) {
				prefix = (prefix << 2) | branch;
				if (depth == (k-1) ) {
					kmers->push_back(prefix);
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
					branches.push(branch);
					branch = 0;
					depth++;
				}
			}
			else {
				while (branch >= 3 && !branches.empty()) {
					depth--;
					branch = branches.top(); branches.pop();
					prefix = prefix >> 2;
				}
				if (branch < 3)
					branch++;
			}
			index++;
		}
		// cerr << endl;
		return kmers;
	}

	/* returns True if the bit_tree (a set) contains the kmers, otherwise 
	returns False. Computes the answer in O(L) where L is the length of the query
	(kmer size) */
	bool contains(const bin_kmer_t bin_kmer) {
		// compute offset into the stream to query for the first letter, and so on
		// computation of ranges over a bit stream
		return false;
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
		cerr << "k = " << (int)k << endl;
		assert(k > 0 && k < 33);
		shared_ptr<boost::dynamic_bitset<>> bitstream(new boost::dynamic_bitset<>());
		// add buffer[1:] to the bitstream
		auto bytes_per_block = sizeof(boost::dynamic_bitset<>::block_type);
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

};

#endif // BIT_TREE_BINARY_LIB