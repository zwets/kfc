/* kmercodec.h
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
#ifndef kmercodec_h_INCLUDED
#define kmercodec_h_INCLUDED

#include <string>
#include <vector>
#include "utils.h"

namespace kfc {


// kmer_codec - encodes and decodes between DNA and k-mers
//
// A kmer_codec is constructed for a fixed kmer size and encoding method.
// The encoding method is canonical (default) or single stranded; see the
// README.md for explanation.
//
// Template parameter kmer_t must be an unsigned integral type capable of
// holding encoded k-mers of the given ksize, plus one reserved bit so that
// invalid k-mers (having non-ACGT content) can be reported.
//
// If constructor argument ksize is greater than kmer_codec<kmer_t>::max_ksize,
// the an exception will be thrown.
//
// The encode(string) member encodes a string of DNA to a vector of kmer
// values.  Any kmer in string that contains an invalid base (not acgtACGT)
// is reported as value kmer_codec<kmer_t>::invalid_kmer.
//
// The decode(kmer_t) member decodes a kmer to a string of DNA.
//
template <typename kmer_t>
class kmer_codec {
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");

    public:
        constexpr static kmer_t invalid_kmer = static_cast<kmer_t>(-1);
        constexpr static int max_ksize = 4*sizeof(kmer_t) - 1;

    private:
        int ksize_;
        bool sstrand_;
        kmer_t max_kmer_;

    public:
        explicit kmer_codec(int ksize, bool sstrand = false);

        kmer_t encode_base(char) const;
        kmer_t encode_kmer(std::string::const_iterator) const;
        std::vector<kmer_t> encode(const std::string&) const;

        std::string decode(kmer_t) const;
};


// implementation ------------------------------------------------------------

template <typename kmer_t>
kmer_codec<kmer_t>::kmer_codec(int ksize, bool sstrand)
    : ksize_(ksize), sstrand_(sstrand), max_kmer_((static_cast<kmer_t>(1)<<(2*max_ksize))-1)
{
    if (ksize < 1)
        raise_error("invalid k-mer size: %d", ksize);

    if (ksize > max_ksize)
        raise_error("k-mer size %d too large for datatype (%d)", ksize, max_ksize);

    if (!sstrand && (ksize & 0x1) != 0x1)
        raise_error("k-mer size must be odd for canonical encoding");
}

template <typename kmer_t>
kmer_t
kmer_codec<kmer_t>::encode_base(char c) const
{
    constexpr static kmer_t X = kmer_codec<kmer_t>::invalid_kmer;
    constexpr static kmer_t base_values[] = { 0, X, 1, X, X, X, 2, X, X, X, X, X, X, X, X, X, X, X, X, 3, X, X, X, X, X, X };

    size_t o = c - 'A';

    if (o < 0 || o > 19) {      // outside 'A'..'T'
        o = c - 'a';            // try 'a'..'t' smallcaps
        if (o < 0 || o > 19)
            return -1;
    }

    return base_values[o];
}

template <typename kmer_t>
kmer_t
kmer_codec<kmer_t>::encode_kmer(const std::string::const_iterator pcur) const
{
    kmer_t res = 0;

    std::string::const_iterator pmid = pcur + (ksize_ / 2);
    std::string::const_iterator pend = pcur + ksize_;

    if (sstrand_) {
        std::string::const_iterator p = pcur - 1;

        while (++p != pend)
            res = (res<<2) | encode_base(*p);
    }
    else {
        if (!(encode_base(*pmid) & 2))  {    // middle base is a or c, encode forward
            std::string::const_iterator p = pcur - 1;

            while (++p != pmid)
                res = (res<<2) | encode_base(*p);

            res = (res<<1) | encode_base(*p); // central base encoded as 1 bit: a->0, c->1

            while (++p != pend)
                res = (res<<2) | encode_base(*p);
        }
        else {
            std::string::const_iterator p = pend;

            while (--p != pmid)
                res = (res<<2) | (encode_base(*p) ^ 3); // xor with 3 is complementary base

            res = (res<<1) | (encode_base(*p) ^ 3); // t->a->0, g->c->1

            while (p-- != pcur)
                res = (res<<2) | (encode_base(*p) ^ 3);
        }
    }

    return res > max_kmer_ ? invalid_kmer : res;
}

template <typename kmer_t>
std::vector<kmer_t>
kmer_codec<kmer_t>::encode(const std::string& s) const
{
    std::vector<kmer_t> result;

    size_t n = s.size() + 1;
    if (static_cast<size_t>(ksize_) < n)
    {
        n -= ksize_;

        result.reserve(n);

        std::string::const_iterator pcur = s.cbegin();
        std::string::const_iterator pend = pcur + n;

        while (pcur != pend) {
            kmer_t kmer = encode_kmer(pcur);
            result.emplace_back(kmer);
            ++pcur;
        }
    }

    return result;
}

template <typename kmer_t>
std::string
kmer_codec<kmer_t>::decode(kmer_t kmer) const
{
    static const char CHARS[4] = { 'a', 'c', 'g', 't' };

    if (kmer == invalid_kmer)
        return "invalid";

    std::string result;
    result.reserve(ksize_);

    if (sstrand_)
        for (int s = 2*ksize_-2; s >= 0; s -= 2)
            result.push_back(CHARS[0x3 & (kmer>>s)]);
    else {
        int s = 2*ksize_ - 3;
        int m = (ksize_ - 1) / 2;

        while (s > m) {
            result.push_back(CHARS[0x3 & (kmer>>s)]);
            s -= 2;
        }

        result.push_back(CHARS[0x1 & (kmer>>m)]);

        --s;
        while (s >= 0) {
            result.push_back(CHARS[0x3 & (kmer>>s)]);
            s -= 2;
        }
    }

    return result;
}


} // namespace kfc

#endif // kmercodec_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
