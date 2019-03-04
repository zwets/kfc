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
#include "tallyman.h"
#include "kmercodec.h"

namespace kfc {


// Forward declarations
template <typename count_t> class kmer_counter;
template <typename kmer_t, typename count_t> class kmer_counter_impl;


// kmer_counter_{L,Q,S,D} - convenience typedefs for kmer_counter with different
// count_t data types.  The count_t type tallies kmer hits, and must be small
// to save on memory, but large enough to not overflow.  It can be integral
// or floating point, though the latter is likely slower.
//
// - L: if there could ever be more than 4G occurrences of a single kmer
// - Q: if you value speed at the cost of possible larger memory consumption
// - S: if you value low memory consumption over possible faster speed
// - D: if you want to try using a double for keeping count
//
typedef kmer_counter<std::uint_fast64_t> kmer_counter_L;
typedef kmer_counter<std::uint_fast32_t> kmer_counter_Q;
typedef kmer_counter<std::uint32_t>      kmer_counter_S;
typedef kmer_counter<double>             kmer_counter_D;


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
// Perform any number of calls to process(seq), then invoke write_results()
// to output the detected kmers and their counts.
//
// This class has two implementations: one fits ksize up to 15 (32-bits), the
// other fits ksize up to 31 (64-bits).  The static create() factory method
// chooses between these depending on ksize and the k32_bit parameter.
//
// By default, if ksize is small, it picks std::uint_fast32_t which may mean
// it uses 64-bit.  With k32_bit = true, it forces 32-bit.  Note however that
// the bit-size of kmer_t isn't as impactful on memory as count_t.
//
// The tallyman component in the implementation classes encapsulates further
// optimisations for speed and memory consumption; see tallyman.h for details.
//
template <typename count_t>
class kmer_counter
{
    static_assert(std::is_integral<count_t>::value || std::is_floating_point<count_t>::value,
            "template argument count_t must be a numerical type");

    protected:
        int ksize_;
        bool s_strand_;
        int n_threads_;

    public:
        static kmer_counter<count_t>* create(int ksize, bool s_strand = false, int max_gb = 0, bool k32_bit = false, int n_threads = 0);

    public:
        kmer_counter(int ksize, bool s_strand, int n_threads);
        virtual ~kmer_counter() { }

        int ksize() const { return ksize_; }
        bool single_strand() const { return s_strand_; }

        virtual void process(const std::string& data) = 0;
        virtual std::ostream& write_results(std::ostream& os, unsigned = output_opts::none) const = 0;
};


// kmer_counter_impl ----------------------------------------------------------
//
// Implements kmer_counter for a given kmer_t.  You can instantiate this class
// directly if you know your ksize ahead of time that ksize will not exceed kmer_codec<kmer_t>::max_ksize.

template <typename kmer_t, typename count_t>
class kmer_counter_impl : public kmer_counter<count_t>
{
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");
    static_assert(std::is_integral<count_t>::value || std::is_floating_point<count_t>::value,
            "template argument count_t must be a numerical type");

    public:
        constexpr static int max_ksize = kmer_codec<kmer_t>::max_ksize;

    private:
        std::unique_ptr<tallyman<kmer_t,count_t> > tallyman_;
        kmer_codec<kmer_t> codec_;

    public:
	kmer_counter_impl(int ksize, bool s_strand, int mem_gb, int n_threads);
        kmer_counter_impl(const kmer_counter_impl<kmer_t,count_t>&) = delete;
        kmer_counter_impl& operator=(const kmer_counter_impl<kmer_t,count_t>&) = delete;
        virtual ~kmer_counter_impl() { }

        virtual void process(const std::string& data);
        virtual std::ostream& write_results(std::ostream& os, unsigned = output_opts::none) const;

    private:
        void write_vec_results(std::ostream&, const count_t*, const count_t*, bool dna, bool zeros) const;
};


// kmer_counter factory  --------------------------------------------------------

template <typename count_t>
kmer_counter<count_t>* kmer_counter<count_t>::create(int ksize, bool s_strand, int max_gb, bool k32_bit, int n_threads)
{
    kmer_counter<count_t> *ret = 0;

    if (ksize < 1)
        raise_error("invalid k-mer size: %d", ksize);
    else if (ksize <= kmer_counter_impl<std::uint32_t,count_t>::max_ksize)
        if (k32_bit)
            ret = new kmer_counter_impl<std::uint32_t,count_t>(ksize, s_strand, max_gb, n_threads);
        else
            ret = new kmer_counter_impl<std::uint_fast32_t,count_t>(ksize, s_strand, max_gb, n_threads);
    else if (ksize <= kmer_counter_impl<std::uint64_t,count_t>::max_ksize)
        ret = new kmer_counter_impl<std::uint64_t,count_t>(ksize, s_strand, max_gb, n_threads);
    else
        raise_error("k-mer size %d is not supported, maximum is %d", ksize, kmer_counter_impl<std::uint64_t,count_t>::max_ksize);

    return ret;
}

// kmer_counter methods -------------------------------------------------------

template <typename count_t>
kmer_counter<count_t>::kmer_counter(int ksize, bool s_strand, int n_threads)
    : ksize_(ksize), s_strand_(s_strand), n_threads_(n_threads)
{
    if (ksize < 1)
        raise_error("invalid k-mer size: %d", ksize);
}

// kmer_counter_impl methods --------------------------------------------------

template <typename kmer_t,typename count_t>
kmer_counter_impl<kmer_t,count_t>::kmer_counter_impl(int ksize, bool s_strand, int mem_gb, int n_threads)
    : kmer_counter<count_t>(ksize, s_strand, n_threads),
      tallyman_(tallyman<kmer_t,count_t>::create(2*ksize-(s_strand?0:1), mem_gb)),
      codec_(ksize, s_strand)
{
    if (ksize > max_ksize)
        raise_error("k-mer size %d too large for this impl (max %d)", ksize, max_ksize);
}

template <typename kmer_t,typename count_t>
void
kmer_counter_impl<kmer_t,count_t>::process(const std::string& data)
{
    std::vector<kmer_t> kmers = codec_.encode(data);
    tallyman_->tally(kmers);
}

template <typename kmer_t, typename count_t>
std::ostream&
kmer_counter_impl<kmer_t, count_t>::write_results(std::ostream &os, unsigned opts) const
{
    if (!os)
        return os;

    bool do_headers = (opts & output_opts::no_headers) == 0;
    bool do_dna = (opts & output_opts::no_dna) == 0;
    bool do_invalid = (opts & output_opts::invalids) != 0;
    bool do_zeros = (opts & output_opts::zeros) != 0;

    int k = kmer_counter<count_t>::ksize_;
    bool s = kmer_counter<count_t>::s_strand_;
    count_t n_invalid = tallyman_->invalid_count();

    if (n_invalid && !do_invalid)
        emit("counted %lu invalid k-mers", static_cast<unsigned long>(n_invalid));

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
        os << "... under construction ...\n";
    }

    if (do_invalid && (n_invalid || do_zeros)) {
        if (do_dna)
            os << "invalid\t";
        os << tallyman_->max_value() + 1 << '\t' << n_invalid << std::endl;
    }

/*
    const tallyman<kmer_t, count_t> *tallyman = counter.get_tallyman();

    std::cout << "invalid" << '\t' << tallyman->invalid_count() << std::endl;

    if (tallyman->is_vec()) {
    const count_t *p1 = p0 + tallyman->max_value() + 1;
	while (++p != p1)
	    if (*p)
		std::cout << *p << '\t' << p - p0 << std::endl;
    }
    else if (tallyman->is_map32()) {
	std::map<tallyman::val32_t,kfc::count_t>::const_iterator p0 = tallyman->get_results_map32().begin();
	std::map<tallyman::val32_t,kfc::count_t>::const_iterator p1 = tallyman->get_results_map32().end();
	for (std::map<tallyman::val32_t,kfc::count_t>::const_iterator p = p0; p != p1; ++p)
	    std::cout << p->first << '\t' << p->second << std::endl;
    }
    else { // is_map64
	std::map<tallyman::val64_t,kfc::count_t>::const_iterator p0 = tallyman->get_results_map64().begin();
	std::map<tallyman::val64_t,kfc::count_t>::const_iterator p1 = tallyman->get_results_map64().end();
	for (std::map<tallyman::val64_t,kfc::count_t>::const_iterator p = p0; p != p1; ++p)
	    std::cout << p->first << '\t' << p->second << std::endl;
    }
*/
    return os;
}


template <typename kmer_t, typename count_t>
void
kmer_counter_impl<kmer_t, count_t>::write_vec_results(std::ostream &os, const count_t *pdata, const count_t *pend, bool dna, bool zeros) const
{
    const count_t *p = pdata - 1;
    kmer_t kmer = 0;
    while (++p != pend) {
        if (*p || zeros) {
            if (dna)
                os << codec_.decode(kmer) << '\t';
            os << kmer << '\t' << *p << std::endl;
        }
        ++kmer;
    }
}


} // namespace kfc

#endif // kmercounter_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
