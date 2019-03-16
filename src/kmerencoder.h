/* kmerencoder.h
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
#ifndef kmerencoder_h_INCLUDED
#define kmerencoder_h_INCLUDED

#include <string>
#include <vector>
#include "bitfiddle.h"
#include "utils.h"

namespace kfc {


// kmer_encoder - encodes and decodes between DNA and k-mers
//
// A kmer_encoder is constructed for a fixed kmer size and encoding method.
// The encoding method is double (default) or single stranded; see the
// README.md for explanation.
//
// Template parameter kmer_t must be an unsigned integral type that can hold
// encoded k-mers of the given ksize plus a reserved bit for invalid k-mers.
// Constant member kmer_encoder<kmer_t>::max_ksize gives the maximum ksize.
//
// The encode() members encode a string of DNA to a vector of kmer values.
// Any kmer that contains an invalid base (not acgtACGT) is returned as an
// as an arbitrary number with its high bit set.  The is_invalid() method
// tests for this.
//
// The decode(kmer_t) member decodes a kmer to a string of DNA.  It decodes
// invalid kmers to strings of "X" bases.
//
// The encode() and decode() members delegate to private function pointers
// which are set in the init_() private member function.  This enables
// swapping out the back-end implementation.  The default implementation
// is in kmerencoder.cpp.
//
template <typename kmer_t>
class kmer_encoder {
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");

    typedef void (*encode_fn)(const char*, const char*, kmer_t*);
    typedef std::string (*decode_fn)(kmer_t, bool);

    public:
        constexpr static unsigned max_ksize = (bitsize<kmer_t> -1 ) / 2;
        constexpr static bool is_invalid(kmer_t kmer) { return kmer & high_bit<kmer_t>; }

    private:
        const unsigned ksize_;
        const bool sstrand_;
        const kmer_t max_kmer_;
        encode_fn encode_;
        decode_fn decode_;

    public:
        explicit kmer_encoder(unsigned ksize, bool sstrand = false);

        kmer_t max_kmer() const { return max_kmer_; }

        void encode(const char*, const char*, kmer_t*) const;
        void encode(const std::string&, kmer_t*) const;
        void encode(std::string&&, kmer_t*) const;

        kmer_t encode_kmer(const char*) const;
        std::vector<kmer_t> encode(const std::string&) const;
        std::vector<kmer_t> encode(std::string&&) const;

        std::string decode(kmer_t, bool rc = false) const;

    private:
        void init_(unsigned, bool);
};


// implementation ------------------------------------------------------------

template <typename kmer_t>
kmer_encoder<kmer_t>::kmer_encoder(unsigned ksize, bool sstrand)
    : ksize_(ksize), sstrand_(sstrand), 
      max_kmer_((((kmer_t)1)<<(2*ksize-(sstrand?0:1)))-1),
      encode_(0), decode_(0)
{
    if (ksize < 1)
        raise_error("invalid k-mer size: %d", ksize);

    if (ksize > max_ksize)
        raise_error("k-mer size %d too large for datatype (%d)", ksize, max_ksize);

    if (!sstrand && (ksize & 0x1) != 0x1)
        raise_error("k-mer size must be odd for double-stranded encoding");

    init_(ksize, sstrand);
}

template <typename kmer_t>
void
kmer_encoder<kmer_t>::encode(const char *pbeg, const char *pend, kmer_t* t) const
{
    if (pbeg + ksize_ <= pend)
        encode_(pbeg, pend, t);
}

template <typename kmer_t>
void
kmer_encoder<kmer_t>::encode(const std::string& s, kmer_t* t) const
{
    if (s.size() >= ksize_)
        encode_(s.data(), s.data() + s.size(), t);
}

template <typename kmer_t>
void
kmer_encoder<kmer_t>::encode(std::string &&s, kmer_t* t) const
{
    if (s.size() >= ksize_)
        encode_(s.data(), s.data() + s.size(), t);
}

template <typename kmer_t>
std::vector<kmer_t>
kmer_encoder<kmer_t>::encode(const std::string& s) const
{
    std::vector<kmer_t> v;

    if (s.size() >= ksize_) {
        v.resize(s.size() - ksize_ + 1);
        encode_(s.data(), s.data() + s.size(), v.data());
    }

    return v;
}

template <typename kmer_t>
std::vector<kmer_t>
kmer_encoder<kmer_t>::encode(std::string&& s) const
{
    return encode((const std::string&)s);
}

template <typename kmer_t>
kmer_t
kmer_encoder<kmer_t>::encode_kmer(const char *s) const
{
    kmer_t kmer;
    encode_(s, s+ksize_, &kmer);
    return kmer;
}

template <typename kmer_t>
std::string
kmer_encoder<kmer_t>::decode(kmer_t kmer, bool rc) const
{
    return decode_(kmer, rc);
}


} // namespace kfc

#endif // kmerencoder_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
