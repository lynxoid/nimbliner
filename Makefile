all: sample mapper
INCLUDES=-I /home/dpellow/libbf/src/
LIBS=/home/dpellow/libbf/build/src/bf/

sample:
	g++ -O3 -std=c++11 -o sample cpp-src/sample_reads.cpp -I /usr/local/include/ -L /usr/local/lib/

mapper:
	g++ -O3 -std=c++11 -o mapper cpp-src/kmer_location.cpp $(INCLUDES) -I /usr/local/include/ -L /usr/local/lib/ -L $(LIBS) -lbf
