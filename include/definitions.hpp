//
// Type definitions
//
//

#ifndef DEFINITIONS
#define DEFINITIONS

#include <utility>

// kmer type
typedef uint64_t kmer_t;

typedef kmer_t bin_kmer_t;

typedef uint genomic_coordinate_t;

typedef uint32_t reference_id_t;

typedef std::pair<reference_id_t,genomic_coordinate_t> seed_position_t;

#endif
