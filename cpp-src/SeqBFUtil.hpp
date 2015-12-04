#ifndef SEQBF_UTIL
#define SEQBF_UTIL

#include <unordered_set>

#include "BaseBloomFilter.hpp"
#include "FastaReader.h"
#include "JellyfishUtil.h"

bool contains(unordered_set<kmer_t> & set, kmer_t t) {
  return set.find(t) != set.end();
}

////////////////////////////////////////////////
// search for kmer extensions to the left of the input
// kmer in bf up to extensions of length len
////////////////////////////////////////////////
bool searchLeft(const kmer_t & kmer, unsigned len, const int k, bf::basic_bloom_filter * bf){
    if(len == 0)
        return true;
    bool contains = false;
    kmer_t suffix = kmer >> 2;
    for(size_t i = 0; i < 4; i++){
      kmer_t extension = suffix | (i<<((k-1)*2)); //TODO: make sure this works - can shift this much?
      if (bf->lookup(extension))
          contains |= searchLeft(extension, len-1, k, bf);
      if(contains)
        return contains;
    }
    return contains;
}


////////////////////////////////////////////////
// search for kmers that are extensions of the input
// of length dist
////////////////////////////////////////////////
bool skipAndSearchLeft(const kmer_t& kmer, unsigned dist, const int k, bf::basic_bloom_filter * bf){
  kmer_t suffix = kmer >> (2*dist); //TODO: put in check that this distance is not too far
  for(unsigned i = 0; i < pow(4,dist); i++){
    kmer_t extension = suffix | (i<<((k-dist)*2));
    if(bf->lookup(extension))
      return true;
  }
  return false;
}

//same for the right
bool searchRight(const kmer_t & kmer, unsigned len, const int k, bf::basic_bloom_filter * bf){
    if(len == 0)
        return true;
    bool contains = false;
    kmer_t prefix = kmer << 2;
    prefix &= (1<<(k*2))-1; //mask out the leftmost 2 bits
    for(size_t i = 0; i < 4; i++){
      kmer_t extension = prefix | i;
      if (bf->lookup(extension))
          contains |= searchRight(extension, len-1, k, bf);
      if(contains)
        return contains;
    }
    return contains;
}


////////////////////////////////////////////////
bool skipAndSearchRight(const kmer_t& kmer, unsigned dist, const int k, bf::basic_bloom_filter * bf){
  kmer_t prefix = kmer<<(2*dist);
  prefix &= (1<<(k*2))-1;
  for(unsigned i = 0; i < pow(4,dist); i++){
    kmer_t extension = prefix | i;
    if(bf->lookup(extension))
      return true;
  }
  return false;
}

////////////////////////////////////////////////
// write out kmers to binary file
// returns false if fopen or fwrite errors
////////////////////////////////////////////////
bool writeKmers(const unordered_set<kmer_t> & kmer_set, const string & fname){
  FILE* outfile;
  size_t buf_size = 1024*1024;
  outfile = fopen((fname+".bin").c_str(),"wb");
  if(!outfile) {
    cerr << "[ERROR] Could not open file for writing" << endl;
    return false; //better error handling?
  }
  size_t i = 0;
  size_t written;
  kmer_t * data_buf = new kmer_t[buf_size];
  for(auto km: kmer_set){
    data_buf[i] = km;
    //filled buffer
    if(i==buf_size-1){
      written = fwrite(data_buf,sizeof(kmer_t),buf_size,outfile);
      if(written != buf_size){ //better error handling
        fclose(outfile);
        delete data_buf;
        return false;
      }
      i = 0;
      continue;
    }
    i++;
  }
  written = fwrite(data_buf,sizeof(kmer_t),i,outfile);
  fclose(outfile);
  delete data_buf;
  return (written == i);
}

////////////////////////////////////////////////
//read in a kmer binary file and convert to kmer set
////////////////////////////////////////////////
bool readKmerFile(unordered_set<kmer_t> & kmer_set, const string & fname){
  FILE* infile;
  size_t buf_size = 1024*1024;
  infile = fopen(fname.c_str(),"rb");
  if(!infile) return false;
  kmer_t * data_buf = new kmer_t[buf_size];
  //get the number of elements
  fseek(infile,0,SEEK_END);
  long long fSize = ftell(infile);
  rewind(infile);
  long long num_kmers = fSize/sizeof(kmer_t);

  size_t numRead;
  while(num_kmers > 0){
    numRead = fread(data_buf,sizeof(kmer_t),buf_size,infile);
    if(numRead < buf_size && numRead < num_kmers){
      delete data_buf;
      fclose(infile);
      return false;
    }
    for(size_t i = 0; i < numRead; i++){
      kmer_set.insert(data_buf[i]);
    }
    num_kmers -= numRead;
  }
  delete data_buf;
  return true;
}

////////////////////////////////////////////////////////////////
// parse a fasta file and return vector of reads
////////////////////////////////////////////////////////////////
vector<string> parseFasta(string const & path) {
    vector<string> reads;
    FastaReader fr(path.c_str());
    kseq_t * seq;
    size_t cnt = 0;
    while ( (seq = fr.nextSequence() ) ) {
      // cerr << seq->seq.s << endl;
      reads.push_back(seq->seq.s);
      cnt++;
    }
    cerr << "(" << cnt << " reads) ";
    return reads;
}

////////////////////////////////////////////////////////////////
// take in reads and output the kmer set
////////////////////////////////////////////////////////////////
unordered_set<kmer_t> getKmers(vector<string> & reads, const int K) {
    unordered_set<kmer_t> kmers;
    for (auto r : reads) {

        if (r.size() < K) continue;
        for (int i = 0; i < r.size() - K + 1; i++) {
            // if (i+K > r.size()) {
            //   cerr << i + K << endl;
            // }
            kmer_t kmer_bin = mer_string_to_binary(&r[i], K);
            // kmers.insert( r.substr(i, K) );
            kmers.insert(kmer_bin);
        }
    }
    // cerr << endl;
    return kmers;
}

// take in reads and output both the kmer set and the edge kmer set
// edge kmers are the kmers that are within extend_len of either end of a read
void getKmersAndEdgeKmers(vector<string> & reads, const int K, const unsigned extend_len, unordered_set<kmer_t> & kmers, unordered_set<kmer_t> & edgeKmers){
  for (auto r : reads) {
    if (r.size() < K) continue;
    for (size_t i = 0; i < r.size() - K + 1; i++) {
      kmer_t kmer_bin = mer_string_to_binary(&r[i], K);
      kmers.insert(kmer_bin);
      if ( i < extend_len || i > r.size()-K-extend_len){
        edgeKmers.insert(kmer_bin);
      }
    }
  }
}

// same, but only taking every (N+1)th kmer
// ** this is to be used for sparse SBFs where the inputs
// are full length sequences, not different reads that could overlap **
void getNthKmersAndEdgeKmers(vector<string> & sequences, const int K, const unsigned skip_len, unordered_set<kmer_t> & kmers, unordered_set<kmer_t> & edgeKmers){
  for (auto seq : sequences) {
    if (seq.size() < K) continue;
    for (size_t i = 0; i < seq.size() - K + 1; i++) {
      kmer_t kmer_bin = mer_string_to_binary(&seq[i], K);
      if(!(i%(skip_len+1))) // N+1th kmer
        kmers.insert(kmer_bin);
      if (i == 0 || i > seq.size()-K-skip_len){
        edgeKmers.insert(kmer_bin);
      }
    }
  }
}


#endif
