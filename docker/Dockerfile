# this dockerfile creates a runnable docker image w/ nimbliner and all its
# dependencies installed

FROM ubuntu:14.04

MAINTAINER dasha.filippova@gmail.com

RUN apt-get update && apt-get install -y wget \
  build-essential \
  zlib1g-dev \
  cmake \
  libboost-all-dev

# install tclap
RUN wget https://sourceforge.net/projects/tclap/files/tclap-1.2.1.tar.gz && \
  tar -xvzf tclap-1.2.1.tar.gz && \
  cd tclap-1.2.1 && \
  ./configure && make && make install

# install libbf
RUN wget https://github.com/mavam/libbf/archive/v0.1-beta.tar.gz && \
  tar -xvzf v0.1-beta.tar.gz && \
  cd libbf-0.1-beta && \
  ./configure && make && make install

ENV LD_LIBRARY_PATH /usr/local/lib/

COPY . /nimbliner/

RUN cd /nimbliner/ && \
  make && make install

WORKDIR /nimbliner
