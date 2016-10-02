BIN=bin
SRC=src
INCLUDE=-I /usr/local/include/ -I include/ # -I $(HOME)/tbb2017_20160916oss/include/
BOOST_INCLUDE=-I /Users/geet/boost_1_52_0/
LIB=-L /usr/local/lib/ # -L $(HOME)/tbb2017_20160916oss/lib/
BOOST_LIB=
FLAGS=-O3 -std=c++11

all: sample mapper indexer bit_tree

indexer:
	g++ $(FLAGS) -o $(BIN)/indexer $(SRC)/index_builder.cpp $(INCLUDE) $(BOOST_INCLUDE) $(LIB)

sample:
	g++ $(FLAGS) -o $(BIN)/sample $(SRC)/sample_reads.cpp $(INCLUDE) $(LIB)

mapper:
	# g++ $(FLAGS) -fopenmp -o $(BIN)/mapper $(SRC)/mapper.cpp $(INCLUDE) $(BOOST_INCLUDE) $(LIB) -lbf
	g++ $(FLAGS) -o $(BIN)/mapper $(SRC)/mapper.cpp $(INCLUDE) $(LIB) -lbf

bit_tree:
	g++ $(FLAGS) -o $(BIN)/bit_tree $(SRC)/bit_tree_binary.cpp $(INCLUDE) $(BOOST_INCLUDE) $(LIB)

clean:
	rm -f sample mapper

run_test: mapper
	./$(BIN)/mapper 20 data/output_mismatches/chr20/sampled/sampled_100_1000000.fa data/output_mismatches/chr20/index/chr20.index data/output_mismatches/chr20/index/chr20.star

docker-build:
	docker build -t nimbliner-dev:0.1 docker/

docker-run:
	docker run -v `pwd`:/nimbliner nimbliner-dev:0.1 make run_test
