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
// or floating point, though the latter is likely to be slower.
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


// kmer_counter
//
// Perform any number of calls to process(seq), then invoke write_results()
// to output the detected kmers and their counts.
//
// This class has two implementations: one for ksize up to 15 (32-bits),
// and one for ksize up to 31 (64-bits).  The static create() factory method
// produces the most memory-efficient implementation for the given ksize.
//
// Note that the 32-bit (memory-efficient) implementation may be slower
// than the 64-bit implementation, so for small ksizes you may still want
// to instantiate the 64-bit implementation, or let the compiler choose
// the fastest implementation.  This is simply done by:
//
//     kmer_counter_impl<std::uint_fast32_t, std::uint_fast32_t> counter(...)
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
        static kmer_counter<count_t>* create(int ksize, bool s_strand = false, int max_gb = 0, int n_threads = 0);

    public:
        kmer_counter(int ksize, bool s_strand, int n_threads);
        virtual ~kmer_counter() { }

        int ksize() const { return ksize_; }
        bool single_strand() const { return s_strand_; }

        virtual void process(const std::string& data) = 0;
        virtual std::ostream& write_results(std::ostream& os) const = 0;
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
        virtual std::ostream& write_results(std::ostream& os) const;
};


// kmer_counter factory  --------------------------------------------------------

template <typename count_t>
kmer_counter<count_t>* kmer_counter<count_t>::create(int ksize, bool s_strand, int max_gb, int n_threads)
{
    kmer_counter<count_t> *ret;

    if (ksize < 1)
        raise_error("invalid k-mer size: %d", ksize);
    else if (ksize <= kmer_counter_impl<std::uint32_t,count_t>::max_ksize)
        ret = new kmer_counter_impl<std::uint32_t,count_t>(ksize, s_strand, max_gb, n_threads);
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
      tallyman_(tallyman<kmer_t,count_t>::create(2*ksize + (s_strand ? 1 : 0), mem_gb)),
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
kmer_counter_impl<kmer_t, count_t>::write_results(std::ostream &os) const
{
    /*
    if (os) {
        if (impl64) {
            if (impl64->is_vec) {
                impl64->get_results_vec();
            }
            else {
                impl64->get_results_map();
            }
        }
        else {
            const tallyman<std::uint32_t,count_t> *t = impl32->tallyman();
        }
    }

    const tallyman<kmer_t, count_t> *tallyman = counter.get_tallyman();

    std::cout << "invalid" << '\t' << tallyman->invalid_count() << std::endl;

    if (tallyman->is_vec()) {
	const kfc::count_t *p0 = tallyman->get_results_vec().get();
	const kfc::count_t *p1 = p0 + tallyman->max_value() + 1;
	const kfc::count_t *p = p0 - 1;
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


} // namespace kfc

#endif // kmercounter_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
