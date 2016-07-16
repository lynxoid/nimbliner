#ifndef SAM_WRITER
#define SAM_WRITER

#include "buffered_output.hpp"

class SamWriter : public BufferedOutput {
	
public:
	SamWriter() {}

	SamWriter(const string path) {
		// todo: open an fd for file output
	}

	// void add_alignment(const vector<genomic_position> & positions, )
	void add_alignment(const genomic_position x, const vector<bool> & matched_kmers, const char * edited_sequence) {
		// TODO: generate CIGAR from this
		// TODO: push bytes onto the buffered output
	}
}

#endif