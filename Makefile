BIN=bin
SRC=src
INCLUDE=-I /usr/local/include/ -I include/ # -I $(HOME)/tbb2017_20160916oss/include/
# BOOST_INCLUDE=-I /home/ubuntu/boost_1.62.0/
LIB=-L /usr/local/lib/ # -L $(HOME)/tbb2017_20160916oss/lib/
BOOST_LIB=
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

run_test: nb-mapper
	./$(BIN)/mapper 20 data/output_mismatches/chr20/sampled/sampled_100_1000000.fa data/output_mismatches/chr20/index/chr20.index data/output_mismatches/chr20/index/chr20.star;

run_sample: nb-sample
	./$(BIN)/sample -m 3 -n 10000000 -l 100 -i /data/human/GRCh38/chr20.fna -o data/output_mismatches/chr20/sampled/sampled_100_10000000_m=3.0pct_d=0.0pct.fa

tests: create-dir
	g++ $(FLAGS) -o $(BIN)/nb_tests $(INCLUDE) tests/anchor_index_test.cpp -lbf
