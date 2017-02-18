# docker shits

docker build -t nimbliner-dev:0.1 -f docker/Dockerfile.dev docker/

docker run -v `pwd`:/nimbliner nimbliner-dev:0.1 make run_test

docker run  -v `pwd`:/nimbliner -v /Users/lynxoid/data/reference/:/data nimbliner-dev:0.1 make run_sample
