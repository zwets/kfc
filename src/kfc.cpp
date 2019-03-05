/* kfc.cpp - kmer frequency count
 * 
 * Copyright (C) 2019  Marco van Zwetselaar <io@zwets.it>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <fstream>
#include <iostream>
#include <string>

#include "kmercounter.h"
#include "seqreader.h"
#include "utils.h"

using namespace kfc;

static const int DEFAULT_KSIZE = 15;
static const int MAX_KSIZE = 32;

static const char USAGE[] = "\n"
"Usage: kfc [OPTIONS] [FILE ...]\n"
"\n"
"  Count the kmers in FILE or from standard input"
"\n"
"  OPTIONS\n"
"   -k KSIZE  k-mer size (default %d), must be odd unless option -s is present\n"
"   -s        consider input to be single stranded, do not canonicalise k-mers\n"
"   -n        output k-mers as bit-encoded numbers only, no DNA sequence output\n"
"   -z        include k-mers with a zero count in the output (default omitted)\n"
"   -i        include the invalid k-mer count in the output (default stderr)\n"
"   -q        suppress output headers\n"
"   -m MEM    constrain memory use to about MEM GB (default: all minus 2GB)\n"
"   -t NUM    use NUM threads (default: all system threads)\n"
"   -x        force k-mers in 32-bits (you normally won't need this option)\n"
"   -v        produce verbose output to stderr\n"
"\n"
"  Each FILE can be an (optionally gzipped) FASTA, FASTQ, or plain DNA file.\n"
"  If FILE is omitted or '-', it is read from stdin\n"
"\n"
"  Only k-mers consisting of proper bases (acgtACGT) are counted.  All k-mers\n"
"  containing other letters are counted as invalid.\n"
"\n"
"  Unless option -s is present, a k-mer and its reverse complement are counted\n"
"  as two occurrences of the same k-mer, whose canonical form is whichever of\n"
"  the two has A or C as its middle base.\n"
"\n"
"  If option -s is present, then k-mers are output as they occur in the input,\n"
"  reverse complements are counted separately, and KSIZE may be odd or even.\n"
"\n"
"  The output has three columns: k-mer dna sequence, k-mer number, count.  The\n"
"  DNA column can be suppressed with option -n.  K-mers with a 0 count are not\n"
"  output unless option -z is present.  Option -i includes the invalid k-mer\n"
"  count in the output (by default printed to standard error)\n"
"\n"
"  More information: http://io.zwets.it/kfc.\n"
"\n";


void
usage_exit()
{
    fprintf(stderr, USAGE, DEFAULT_KSIZE);
    std::exit(1);
}

int main (int, char *argv[]) 
{
    int ksize = DEFAULT_KSIZE;
    bool single_strand = false;
    int max_mem = 0;
    int kmer_32bit = false;
    int n_threads = 0;
    unsigned o_opts = output_opts::none;

    set_progname("kfc");

        // Parse arguments

    while (*++argv && **argv == '-' && (*argv)[1] != '\0')
    {
        char opt = (*argv)[1];

        if (opt == '-' && (*argv)[2] == '\0') {
            ++argv;            // double dash marks end of options
            break;
        }
        else if (opt == 'v') {
            set_verbose(true);
        }
        else if (opt == 's') {
            single_strand = true;
        }
        else if (opt == 'x') {
            kmer_32bit = true;
        }
        else if (opt == 'n') {
            o_opts |= output_opts::no_dna;
        }
        else if (opt == 'z') {
            o_opts |= output_opts::zeros;
        }
        else if (opt == 'i') {
            o_opts |= output_opts::invalids;
        }
        else if (opt == 'q') {
            o_opts |= output_opts::no_headers;
        }
        else if (!*++argv) {    // subsequent options require argument
            usage_exit();
        }
        else if (opt == 'k') {
            ksize = std::atoi(*argv);
            if (ksize < 1 || ksize > MAX_KSIZE) 
                raise_error("invalid k-mer size: %s", *argv);
        }
        else if (opt == 'm') {
            if ((max_mem = std::atoi(*argv)) < 1)
                raise_error("invalid memory size: %s", *argv);
        }
        else if (opt == 't') {
            if ((n_threads = std::atoi(*argv)) < 1)
                raise_error("invalid number of threads: %s", *argv);
        }
        else
            usage_exit();
    }

        // Create the kmer_counter

    // It seems the 32-bit count_t is faster than the fast32
    //
    //std::unique_ptr<kmer_counter_Q> counter(kmer_counter_Q::create(ksize, single_strand, max_mem, kmer_32bit, n_threads));
    std::unique_ptr<kmer_counter_S> counter(kmer_counter_S::create(ksize, single_strand, max_mem, kmer_32bit, n_threads));

        // Iterate over files

    std::string fname(*argv ? *argv++ : "-");

    do {
        verbose_emit("reading file: %s", fname.c_str());

        std::istream *is = &std::cin;

        std::ifstream in_file;
        if (fname != "-") {
            in_file.open(fname, std::ios_base::in|std::ios_base::binary);
            if (!in_file)
                raise_error("failed to open file: %s", fname.c_str());
            is = &in_file;
        }

        sequence_reader reader(*is);
        sequence seq;

        while (reader.next(seq))
            counter->process(seq.data);

        in_file.close();

        fname = *argv ? *argv++ : "";

    } while (!fname.empty());

        // Output kmer_counter results

    counter->write_results(std::cout, o_opts);

    return 0;
}

// vim: sts=4:sw=4:et:si:ai
