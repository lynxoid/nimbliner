#ifndef BUFFERED_OUTPUT
#define BUFFERED_OUTPUT

#include <iostream>
#include <fstream>

#include <vector>
#include <queue>

#include <fcntl.h>
#include <unistd.h>

// #include <compress.h>

// #include "IntervalTree.h"

const mode_t usr_rw = S_IRUSR | S_IWUSR;
const mode_t all_rw = usr_rw | S_IRGRP | S_IWGRP | S_IROTH | S_IWOTH;
mode_t outfd_mode = usr_rw;

#ifdef O_BINARY
const int o_binary = O_BINARY;
#else
const int o_binary = 0;
#endif

using namespace std;

class BufferedOutput {

	// vector<char> bytes_;
	deque<uint8_t> bytes_;

	int64_t total_bytes = 0;

	int out_fd; // output file descriptor

	int dictionary_size = 1<<23;

	int match_len_limit = 36; // equivalent to -6 option

	bool timeToDump() {
		return data.size() >= dictionary_size;
	}

	////////////////////////////////////////////////////////////////
	//
	////////////////////////////////////////////////////////////////
	/*
	void compressAndWriteOut() {
		int data_size = data.size();
		total_bytes += data_size;

		// TODO: is this usually one packet?
		// can we avoid memcpy and just give away this data vector to the packet?
		int max_size = 0;
		if (flush_all) 
			max_size = data_size;
		else 
			max_size = (data_size / dictionary_size) * dictionary_size;

		for (int ate = 0; ate < max_size; ate += dictionary_size) {
			// cerr << "sending the block to PLZIP" << endl;
			int remainder = data_size - ate;
			int size = std::min(dictionary_size, remainder);
			uint8_t * block_data = new( std::nothrow ) uint8_t[ size ];
			int i;
			for (i = 0; i < size; i++) block_data[i] = data[i];
			courier->receive_packet( block_data, size, out_fd ); // associate an output stream to the packet
			// erase the elements that we sent to a packet
			// pop_front is better for this
			i = 0;
			while (i < size) {
				i++;
				data.pop_front();
			}
		}
	}
	*/

public:
	
	BufferedOutput() {}

	BufferedOutput(int d = 1<<23, int match_len = 36): 
		dictionary_size(d),
		match_len_limit(match_len) {
		int flags = O_CREAT | O_WRONLY | o_binary;
		out_fd = open( (fn + suff).c_str(), flags, outfd_mode );
		cerr << "Opened stream for output data (fd=" << out_fd << ")" << endl;
	}

	///////////////////////////////////////////////////////////
	~BufferedOutput() {
		// cerr << "Closing fd=" << out_fd << ", processed " << total_bytes << " bytes; ";
		close(out_fd);
	}

	void add_to_buffer(vector<char> & v) {
		// add to buffer
		// if buffer full - flush to disk
	}

	////////////////////////////////////////////////////////////////
	// stream size
	int size() { return data.size(); }

	void flush() {
		// TODO: 
		// cerr << "flushing " << stream_suffix << endl;
		// compressAndWriteOut();
	}
}