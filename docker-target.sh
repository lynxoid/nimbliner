# docker shits

docker build -t lynxoid/nimbliner:0.1 -f docker/Dockerfile .

docker run -v $PWD/data:/nimbliner/data/ lynxoid/nimbliner:0.1 ./bin/nb_tests

# docker run  -v `pwd`:/nimbliner -v /Users/lynxoid/data/reference/:/data nimbliner-dev:0.1 make run_sample
