/* kmercounter.h
 *
 * Copyright (C) 2018  Marco van Zwetselaar <io@zwets.it>
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
#ifndef kmercounter_h_INCLUDED
#define kmercounter_h_INCLUDED

#include <string>
#include <ostream>
#include <memory>
#include <algorithm>
#include "tallyman.h"
#include "kmerencoder.h"

namespace kfc {


// output_opts - bit field flag for kmer_counter<>::write_results()
//
enum output_opts : unsigned {
    none=0,         // default: three tab-separated columns: dna, kmer, count
    no_dna=1,       // omit the first (DNA string) column
    no_headers=2,   // omit the header line(s)
    zeros=4,        // include k-mers with zero counts
    invalids=8      // include a pseudo-k-mer with the invalid count
};


// kmer_counter
//
// Counts distinct kmers in any number of sequences of DNA.  Writes the counts
// as a table, sorted alphabetically on k-mer, to an output stream.  Its two
// core methods are process(sequence) and write_results().
//
// This class has three implementations: two based on tallying the kmers as
// they are processed, of which one uses a vector of tallies and one uses a
// map of kmer to tally, and one which does not tally but collects the list
// of kmers as is, then sorts this when results are requested.
//
// The count_t template parameter determines the size of the counts that can
// be kept, and has memory impact: factor 2 with k-space for the vector
// implementation, factor 2 with k-count for the map, none for the list.
//
// The kmer_t template parameter has impact factor 2 with k-count for both map
// and list, none for vector.  When k-size > 15, then it must be 64-bit (as
// each base in a k-mer takes up two bits and we reserve one for reporting
// invalid k-mers).  When k-size is 15 or below, std::uint32_t seems to be the
// more optimal choice than std::uint32_fast_t (which may be 64-bit).
//
// The tallyman component in the implementation classes encapsulates further
// optimisations for speed and memory consumption; see tallyman.h for details.
//
// The static "create" factory method attempts to return the optimal
// implementation class (vector, map, or list).  Actual performance is hard to
// predict. @TODO@ add guidelines.
//
class kmer_counter
{
    protected:
        int ksize_;
        bool s_strand_;

    public:
        kmer_counter(int ksize, bool s_strand);
        virtual ~kmer_counter() { }

        int ksize() const { return ksize_; }
        bool single_strand() const { return s_strand_; }

        virtual void process(const std::string& data) = 0;
        virtual void process(std::string &&data) = 0;
        virtual std::ostream& write_results(std::ostream& os, unsigned = output_opts::none) const = 0;
};


// kmer_counter_tally ----------------------------------------------------------
//
// Implements kmer_counter by keeping a tally for every kmer.  The tally counter
// has two possible implementations: a vector with an entry for every possible
// value of kmer_t, or a map whose keys are k-mers and values are counts.

template <typename kmer_t, typename count_t>
class kmer_counter_tally : public kmer_counter
{
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");
    static_assert(std::is_integral<count_t>::value || std::is_floating_point<count_t>::value,
            "template argument count_t must be a numerical type");

    public:
        constexpr static int max_ksize = kmer_encoder<kmer_t>::max_ksize;

    private:
        std::unique_ptr<tallyman<kmer_t,count_t>> tallyman_;
        kmer_encoder<kmer_t> encoder_;

    public:
	kmer_counter_tally(tallyman<kmer_t,count_t>*, int ksize, bool s_strand);
        kmer_counter_tally(const kmer_counter_tally<kmer_t,count_t>&) = delete;
        kmer_counter_tally& operator=(const kmer_counter_tally<kmer_t,count_t>&) = delete;
        virtual ~kmer_counter_tally() { }

        virtual void process(const std::string& data);
        virtual void process(std::string &&data);
        virtual std::ostream& write_results(std::ostream& os, unsigned = output_opts::none) const;

    private:
        void write_vec_results(std::ostream&, const count_t*, const count_t*, bool dna, bool zeros) const;
        void write_map_results(std::ostream&, bool dna, bool zeros) const;
};


// kmer_counter_list ------------------------------------------------------------
//
// This implementation does not keep tallies but instead keeps the list of kmers
// as they are coming in.  When write_results is called, the list is sorted and
// the counted kmers are output.

template <typename kmer_t>
class kmer_counter_list : public kmer_counter
{
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");

    public:
        constexpr static int max_ksize = kmer_encoder<kmer_t>::max_ksize;

    private:
        kmer_t *kmers_, *pkmers_cur_, *pkmers_end_;
        kmer_encoder<kmer_t> encoder_;

    public:
	kmer_counter_list(int ksize, bool s_strand, size_t max_count);
        kmer_counter_list(const kmer_counter_list<kmer_t>&) = delete;
        kmer_counter_list& operator=(const kmer_counter_list<kmer_t>&) = delete;
        virtual ~kmer_counter_list();

        virtual void process(const std::string& data);
        virtual void process(std::string &&data);
        virtual std::ostream& write_results(std::ostream& os, unsigned = output_opts::none) const;
};

// kmer_counter_tally methods --------------------------------------------------

template <typename kmer_t,typename count_t>
kmer_counter_tally<kmer_t,count_t>::kmer_counter_tally(
        tallyman<kmer_t,count_t>* tman, int ksize, bool s_strand)
    : kmer_counter(ksize, s_strand),
      tallyman_(tman),
      encoder_(ksize, s_strand)
{
    if (ksize > max_ksize)
        raise_error("k-mer size %d too large for this impl (max %d)", ksize, max_ksize);
}

template <typename kmer_t,typename count_t>
void
kmer_counter_tally<kmer_t,count_t>::process(const std::string& data)
{
    std::vector<kmer_t> kmers = encoder_.encode(data);
    tallyman_->tally(kmers);
}

template <typename kmer_t,typename count_t>
void
kmer_counter_tally<kmer_t,count_t>::process(std::string &&data)
{
    std::vector<kmer_t> kmers = encoder_.encode(data);
    tallyman_->tally(kmers);
}

template <typename kmer_t, typename count_t>
std::ostream&
kmer_counter_tally<kmer_t, count_t>::write_results(std::ostream &os, unsigned opts) const
{
    if (!os)
        return os;

    bool do_headers = (opts & output_opts::no_headers) == 0;
    bool do_dna = (opts & output_opts::no_dna) == 0;
    bool do_invalid = (opts & output_opts::invalids) != 0;
    bool do_zeros = (opts & output_opts::zeros) != 0;

    int k = kmer_counter::ksize_;
    bool s = kmer_counter::s_strand_;
    count_t n_invalid = tallyman_->invalid_count();

    if (do_headers) {
        // Line 1
        os << "# kfc " << k << "-mer counts "
            << (s ? "(single strand directional)": "(canonical, destranded)" );
        if (!do_invalid && n_invalid) // if !do_invalid then show in header
            os << "; excluding " << n_invalid << " invalid k-mers";
        if (!do_zeros)
            os << "; omitting zero counts";
        os << std::endl;
        // Line 2
        os << "#";
        if (do_dna) os << "k-mer\t";
        os << (s ? "s-code" : "c-code") << '\t' << "count" << std::endl;
    }

    if (tallyman_->is_vec()) {
        const count_t *data = tallyman_->get_results_vec();
        write_vec_results(os, data, data + tallyman_->max_value() + 1, do_dna, do_zeros);
    }
    else {
        write_map_results(os, do_dna, do_zeros);
    }

    if (do_invalid && (n_invalid || do_zeros)) {
        if (do_dna)
            os << "invalid\t";
        os << encoder_.max_kmer() + 1 << '\t' << n_invalid << std::endl;
    }

    if (n_invalid)
        verbose_emit("counted %lu invalid k-mers", static_cast<unsigned long>(n_invalid));

    return os;
}

template <typename kmer_t, typename count_t>
void
kmer_counter_tally<kmer_t, count_t>::write_vec_results(std::ostream &os, const count_t *pdata, const count_t *pend, bool dna, bool zeros) const
{
    const count_t *p = pdata - 1;
    kmer_t kmer = 0;
    while (++p != pend) {
        if (*p || zeros) {
            if (dna)
                os << encoder_.decode(kmer) << '\t';
            os << kmer << '\t' << *p << std::endl;
        }
        ++kmer;
    }
}

template <typename kmer_t, typename count_t>
void
kmer_counter_tally<kmer_t, count_t>::write_map_results(std::ostream &os, bool dna, bool zeros) const
{
    const std::map<kmer_t,count_t>& map = tallyman_->get_results_map();

    typename std::map<kmer_t,count_t>::const_iterator p = map.begin();
    typename std::map<kmer_t,count_t>::const_iterator pend = map.end();

    if (!zeros) {
        while (p != pend) {
            if (dna)
                os << encoder_.decode(p->first) << '\t';
            os << p->first << '\t' << p->second << std::endl;
            ++p;
        }
    }
    else {
        kmer_t kmer = 0;
        kmer_t done_kmer = tallyman_->max_value() + 1;
        kmer_t next_kmer = p == pend ? done_kmer : p->first;

        while (kmer != done_kmer) {

            while (kmer != next_kmer) {
                if (dna)
                    os << encoder_.decode(kmer) << '\t';
                os << kmer << "\t0" << std::endl;
                ++kmer;
            }

            if (kmer != done_kmer) { // so it is next_kmer and p->first
                if (dna)
                    os << encoder_.decode(kmer) << '\t';
                os << kmer << '\t' << p->second << std::endl;
                next_kmer = ++p == pend ? done_kmer : p->first;
                ++kmer;
            }
        }
    }
}


// kmer_counter_list methods --------------------------------------------------

template <typename kmer_t>
kmer_counter_list<kmer_t>::kmer_counter_list(int ksize, bool s_strand, size_t max_count)
    : kmer_counter(ksize, s_strand),
      kmers_(0), pkmers_cur_(0), pkmers_end_(0),
      encoder_(ksize, s_strand)
{
    if (ksize > max_ksize)
        raise_error("k-mer size %d too large for this impl (max %d)", ksize, max_ksize);

    //kmers_ = (kmer_t*) malloc(max_count * sizeof(kmer_t));
    kmers_ = ::new kmer_t[max_count];

    if (kmers_) {
        pkmers_cur_ = kmers_;
        pkmers_end_ = kmers_ + max_count;
    }
    else
        raise_error("failed to allocate memory (%uMB) for k-mer list",
                static_cast<unsigned>((max_count * sizeof(kmer_t)) >> 20));
}

template <typename kmer_t>
kmer_counter_list<kmer_t>::~kmer_counter_list()
{ 
    if (kmers_) 
    //    free (kmers_);
        delete [] kmers_;
}

template <typename kmer_t>
void
kmer_counter_list<kmer_t>::process(const std::string& data)
{
    size_t len = data.size() + 1;
    if (static_cast<size_t>(ksize_) < len)
        len -= ksize_;
    else
        return;

    kmer_t *new_pcur = pkmers_cur_ + len;

    if (new_pcur <= pkmers_end_) {
        // first bump the pcur, so later next thread can enter before encode
        kmer_t *encode_ptr = pkmers_cur_;
        pkmers_cur_ = new_pcur;
        encoder_.encode(data, encode_ptr);
    }
    else
        raise_error("k-mer list capacity (%uM k-mers) exhausted",
                static_cast<unsigned>((pkmers_end_ - kmers_) >> 20));
}

template <typename kmer_t>
void
kmer_counter_list<kmer_t>::process(std::string &&data)
{
    process(data);
}

template <typename kmer_t>
std::ostream&
kmer_counter_list<kmer_t>::write_results(std::ostream &os, unsigned opts) const
{
    // // shrink memory to be nice (can't here, we're const)
    //
    // size_t cur_count = pkmers_cur_ - kmers_;
    // kmers_ = (kmer_t*) std::realloc(kmers_, cur_count * sizeof(kmer_t));
    // pkmers_cur_ = pkmers_end_ = kmers_ + cur_count;

    if (!os)
        return os;

    bool do_headers = (opts & output_opts::no_headers) == 0;
    bool do_dna = (opts & output_opts::no_dna) == 0;
    bool do_invalid = (opts & output_opts::invalids) != 0;
    bool do_zeros = (opts & output_opts::zeros) != 0;

    const size_t k = kmer_counter::ksize_;
    const bool s = kmer_counter::s_strand_;

    if (do_headers) {
        // Line 1
        os << "# kfc " << k << "-mer counts "
            << (s ? "(single strand directional)": "(canonical, destranded)" );
        if (!do_invalid) // if !do_invalid then show in header
            os << "; excluding invalid k-mers";
        if (!do_zeros)
            os << "; omitting zero counts";
        os << std::endl;
        // Line 2
        os << "#";
        if (do_dna) os << "k-mer\t";
        os << (s ? "s-code" : "c-code") << '\t' << "count" << std::endl;
    }

    std::uint64_t n_invalid = 0;
    kmer_t *p = kmers_;

    if (p && p != pkmers_cur_) {

        std::sort(kmers_, pkmers_cur_);

        if (encoder_.is_invalid(*p)) {
            // we hit invalid right away, but since we're sorted
            // this means we are done (invalid is beyond max)
            n_invalid = pkmers_cur_ - p;
            if (do_zeros)
                for (kmer_t i = 0; i <= encoder_.max_kmer(); ++i) {
                    if (do_dna) os << encoder_.decode(i) << '\t';
                    os << i << '\t' << 0 << std::endl;
                }
        }
        else {

            kmer_t last = *p;
            std::uint64_t count = 1;

            // optionally generate zeros for kmers 0..last-1
            if (do_zeros)
                for (kmer_t i = 0; i < last; ++i) {
                    if (do_dna) os << encoder_.decode(i) << '\t';
                    os << i << '\t' << 0 << std::endl;
                }

            while (++p != pkmers_cur_) {
                if (*p == last) {   // continue run of last value
                    ++count;
                }
                else if (encoder_.is_invalid(*p)) {
                    n_invalid += pkmers_cur_ - p;
                    break;
                }
                else {
                    if (do_dna) os << encoder_.decode(last) << '\t';
                    os << last << '\t' << count << std::endl;

                    if (do_zeros)
                        while (++last != *p) {
                            if (do_dna) os << encoder_.decode(last) << '\t';
                            os << last << '\t' << 0 << std::endl;
                        }

                    last = *p;
                    count = 1;
                }
            }

            if (do_dna) os << encoder_.decode(last) << '\t';
            os << last << '\t' << count << std::endl;

            if (do_zeros)
                while (++last <= encoder_.max_kmer()) {
                    if (do_dna) os << encoder_.decode(last) << '\t';
                    os << last << '\t' << 0 << std::endl;
                }
        }
    }
    else if (do_zeros) {
        for (kmer_t i = 0; i <= encoder_.max_kmer(); ++i) {
            if (do_dna) os << encoder_.decode(i) << '\t';
            os << i << '\t' << 0 << std::endl;
        }
    }

    if (do_invalid && (n_invalid || do_zeros)) {
        if (do_dna) os << "invalid\t";
        os << encoder_.max_kmer() + 1 << '\t' << n_invalid << std::endl;
    }

    if (n_invalid) {
        verbose_emit("counted %lu k-mers, %lu invalid", 
                static_cast<unsigned long>(pkmers_cur_ - kmers_), static_cast<unsigned long>(n_invalid));
    }

    return os;
}


} // namespace kfc

#endif // kmercounter_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
