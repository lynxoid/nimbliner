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

		string rnext = "*", read_qualities = "*";
		int mapq = 0, pnext = 0, tlen = 0;

		// write to STDIN by defualt (can pipe into samtools)
		cout << read_id << "\t" << 
				flags << "\t" << 
				reference_name << "\t" << 
				reference_position << "\t" << 
				mapq /* mapq */ << "\t" <<
				cigar << "\t" <<
				rnext << "\t"
				pnext << "\t"
				tlen << "\t"
				read_sequence << "\t" <<
				read_qualities << "\t" << 
				endl;
	}
}

#endif