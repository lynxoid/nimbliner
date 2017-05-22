#ifndef FOFN_READER_H
#define FOFN_READER_H

namespace nimble {

class FOFNReader {

    ifstream fofn_in;
    shared_ptr<FastaReader> current_fasta;

    shared_ptr<FastaReader> open_next_fasta() {
        string fasta_path;
        if ( !getline(fofn_in, fasta_path) ) {
            cerr << "Parsed all reference sequences" << endl;
            return nullptr;
        }
        cerr << "Parsing " << fasta_path << endl;
        current_fasta = shared_ptr<FastaReader>(new FastaReader( fasta_path.c_str() ) );
        if ( !current_fasta->is_open() ) {
            cerr << "[ERROR] Can not read from " << fasta_path << endl;
            current_fasta = nullptr;
        }
        return current_fasta;
    }

public:

    FOFNReader(const string & fofn_path) {
        fofn_in.open(fofn_path);
        can_read_or_quit(fofn_in, fofn_path, true);
    }

    ~FOFNReader() {
        current_fasta.reset();
        fofn_in.close();
    }

    kseq_t * getNextSequence() {
        kseq_t * seq = nullptr;
        if (current_fasta == nullptr) {
            // first time
            current_fasta = open_next_fasta();
            if (current_fasta == nullptr)
                // cant read anymore
                return nullptr;
        }
        seq = current_fasta->nextSequence();
        if (seq == nullptr) {
            current_fasta = open_next_fasta();
            if (current_fasta == nullptr)
                // cant read anymore
                return nullptr;
            seq = current_fasta->nextSequence();
        }
        return seq;
    }
};

}
#endif
