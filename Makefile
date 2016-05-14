BIN=bin
SRC=cpp-src
INCLUDE=-I /usr/local/include/
LIB=-L /usr/local/lib/
FLAGS=-O3 -std=c++11

all: sample mapper

sample:
	g++ $(FLAGS) -o $(BIN)/sample $(SRC)/sample_reads.cpp $(INCLUDE) $(LIB)

mapper:
	g++ $(FLAGS) -o $(BIN)/mapper $(SRC)/kmer_location.cpp $(INCLUDE) $(LIB) -lbf

clean:
	rm -f sample mapper