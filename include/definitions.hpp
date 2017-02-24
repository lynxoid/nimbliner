//
// Type definitions
//
//

#ifndef DEFINITIONS
#define DEFINITIONS

#include <utility>
#include <iostream>
#include <fstream>
#include <string>

// kmer type
typedef uint64_t kmer_t;

typedef kmer_t bin_kmer_t;

typedef uint genomic_coordinate_t;

typedef uint32_t reference_id_t;

typedef std::pair<reference_id_t,genomic_coordinate_t> seed_position_t;

inline bool can_read_or_quit(const std::ifstream & in, const std::string & name, const bool quit = true) {
    if (!in) {
        std::cerr << "[ERROR] Can not read from file " << name << std::endl;
        if (quit)
            exit(1);
        return false;
    }
    return true;
}

inline void can_write_or_quit(const std::ofstream & in) {

}


#endif
