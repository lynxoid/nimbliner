BIN=bin
SRC=src
INCLUDE=-I /usr/local/include/ -I include/
BOOST_INCLUDE=-I /Users/geet/boost_1_52_0/
LIB=-L /usr/local/lib/
BOOST_LIB=
FLAGS=-O3 -std=c++11

all: sample mapper indexer bit_tree

indexer:
	g++ $(FLAGS) -o $(BIN)/indexer $(SRC)/index_builder.cpp $(INCLUDE) $(BOOST_INCLUDE) $(LIB)

sample:
	g++ $(FLAGS) -o $(BIN)/sample $(SRC)/sample_reads.cpp $(INCLUDE) $(LIB)

mapper:
	g++ $(FLAGS) -o $(BIN)/mapper $(SRC)/mapper.cpp $(INCLUDE) $(BOOST_INCLUDE) $(LIB) -lbf

bit_tree:
	g++ $(FLAGS) -o $(BIN)/bit_tree $(SRC)/bit_tree_binary.cpp $(INCLUDE) $(BOOST_INCLUDE) $(LIB)

clean:
	rm -f sample mapper
