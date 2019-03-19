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
#ifndef NO_THREADS
#  include <thread>
#  include <mutex>
#endif
#include "implpicker.h"
#include "seqreader.h"
#include "utils.h"

using namespace kfc;

#ifndef NO_THREADS
static void process_thread(sequence_reader* reader, kmer_counter* counter)
{
    static std::mutex read_mutex;
    sequence seq;

    while (true) 
    {
        std::unique_lock<std::mutex> read_lock(read_mutex);

        if (!reader->next(seq))
            break;
        
        read_lock.unlock();

        counter->process(std::move(seq.data));
    }
}
#endif

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
"   -z        include k-mers with a zero count in the output (default: omit)\n"
"   -i        include the invalid k-mer count in the output (default: stderr)\n"
"   -n        output k-mers as encoded numbers only, no DNA sequences\n"
"   -q        suppress output headers, just show k-mers and counts\n"
"   -l MBASE  limit counting capacity to MBASE million bases (optimises speed)\n"
"   -m MEMGB  constrain memory use to about MEM GB (default: all minus 2GB)\n"
#ifndef NO_THREADS
"   -t NUM    use NUM threads (default: all system threads)\n"
#endif
"   -x l|v|m  override the implementation choice to be list, vector, or map\n"
"   -v        produce verbose output to stderr\n"
"\n"
"  Each FILE can be an (optionally gzipped) FASTA, FASTQ, or plain DNA file.  If\n"
"  FILE is omitted or '-', input is read from stdin.  With gzipped input, note\n"
"  that 'gunzip | kfc' is often fastest (due to multiprocessing).\n"
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
"  Use option -l for speed gains by telling kfc how much input to expect (in\n"
"  millions of bases).  E.g. for a bacterial assembly, '-l 10' will usually\n"
"  suffice, whereas for human use '-l 3200'.\n"
"\n"
"  The output has three columns: k-mer dna sequence, k-mer number, count.  The\n"
"  DNA column can be suppressed with option -n.  K-mers with a 0 count are not\n"
"  output unless option -z is present.  Option -i includes the invalid k-mer\n"
"  count in the output (with DNA sequence \"XXX..\").  By default the invalid\n"
"  count is printed to standard error).\n"
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
    unsigned max_mbp = 0.0;
    unsigned max_gb = 0;
    char force_impl = '\0';
#ifndef NO_THREADS
    int n_threads = 0;
#endif
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
        // subsequent options require an argument
        else if (!*++argv) {
            usage_exit();
        }
        else if (opt == 'k') {
            ksize = std::atoi(*argv);
            if (ksize < 1 || ksize > MAX_KSIZE) 
                raise_error("invalid k-mer size: %s", *argv);
        }
        else if (opt == 'l') {
            if ((max_mbp = std::atoi(*argv)) < 1)
                raise_error("invalid input size: %s", *argv);
        }
        else if (opt == 'm') {
            if ((max_gb = std::atoi(*argv)) < 1)
                raise_error("invalid memory size: %s", *argv);
        }
        else if (opt == 't') {
#ifndef NO_THREADS
            if ((n_threads = std::atoi(*argv)) < 1)
                raise_error("invalid number of threads: %s", *argv);
#else
            raise_error("not compiled with thread support");
#endif
        }
        else if (opt == 'x') {
            switch (force_impl = *argv[0]) {
                case 'l': case 'v': case 'm': break;
                default: raise_error("invalid implementation: %c", force_impl);
            }
        }
        else
            usage_exit();
    }

        // Create the kmer_counter via the pick_implementation method

    std::unique_ptr<kmer_counter> counter(pick_implementation(ksize, single_strand, max_mbp, max_gb, force_impl));

#ifndef NO_THREADS
    if (n_threads == 0) {
        n_threads = get_system_threads();
        verbose_emit("defaulting to %u system threads", n_threads);
    }
    else {
        verbose_emit("using %u threads", n_threads);
    }

    std::vector<std::thread> threads;
    threads.reserve(n_threads);
#endif
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

#ifndef NO_THREADS
        if (n_threads > 1) {
            verbose_emit("spawning %d threads", n_threads);
            for (int i = 0; i < n_threads; ++i)
                threads.push_back(std::thread(&process_thread, &reader, counter.get()));

            for (int i = 0; i < n_threads; ++i)
                threads[i].join();

            threads.clear();
        }
        else {
            sequence seq;

            while (reader.next(seq))
                counter->process(std::move(seq.data));
        }
#else
        sequence seq;

        while (reader.next(seq))
            counter->process(std::move(seq.data));
#endif
        in_file.close();

        fname = *argv ? *argv++ : "";

    } while (!fname.empty());

        // Output kmer_counter results

    counter->write_results(std::cout, o_opts);

    return 0;
}

// vim: sts=4:sw=4:et:si:ai
