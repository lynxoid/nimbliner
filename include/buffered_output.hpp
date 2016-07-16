#ifndef BUFFERED_OUTPUT
#define BUFFERED_OUTPUT

#include <iostream>
#include <fstream>

class BufferedOutput {

	vector<char> bytes_;

public:
	
	BufferedOutput() {}

	void add_to_buffer(vector<char> & v) {
		// add to buffer
		// if buffer full - flush to disk
	}
}

#endif // BufferedOutput