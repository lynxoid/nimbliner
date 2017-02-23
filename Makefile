BIN=bin
SRC=src
INCLUDE=-I /usr/local/include/ -I include/ # -I $(HOME)/tbb2017_20160916oss/include/
# BOOST_INCLUDE=-I /usr/include/
LIB=-L /usr/local/lib/ # -L $(HOME)/tbb2017_20160916oss/lib/
BOOST_LIB=/usr/lib/x86_64-linux-gnu/
FLAGS=-O3 -std=c++11
# DEBUGFLAGS=-DWITHGPERFTOOLS -g -fno-inline
# enable debug symbols and no inlining for better call stack info
# DEBUGFLAGS=-g -fno-inline
# DEBUGFLAGS=-g

all: create-dir nb-sample nb-mapper nb-indexer nb-bit_tree tests

create-dir:
	mkdir -p $(BIN)

nb-indexer: create-dir
	g++ $(FLAGS) -o $(BIN)/indexer $(SRC)/index_builder.cpp $(INCLUDE) $(BOOST_INCLUDE) $(LIB)

nb-sample: create-dir
	g++ $(FLAGS) -o $(BIN)/sample $(SRC)/sample_reads.cpp $(INCLUDE) $(LIB)

nb-mapper: create-dir
	# g++ $(FLAGS) -fopenmp -o $(BIN)/mapper $(SRC)/mapper.cpp $(INCLUDE) $(BOOST_INCLUDE) $(LIB) -lbf
	g++ $(FLAGS) $(DEBUGFLAGS) -o $(BIN)/mapper $(SRC)/mapper.cpp $(INCLUDE) $(BOOST_INCLUDE) $(LIB) -lbf

nb-bit_tree: create-dir
	g++ $(FLAGS) -o $(BIN)/bit_tree $(SRC)/bit_tree_binary.cpp $(INCLUDE) $(BOOST_INCLUDE) $(LIB)

clean:
	rm -f $(BIN)/*

install:
	cp $(BIN)/* /usr/local/bin/

tests: create-dir
	g++ $(FLAGS) -o $(BIN)/nb_tests $(INCLUDE) tests/*_test.cpp -lbf
