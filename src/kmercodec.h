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
// Template parameter kmer_t must be an unsigned integral type that can hold
// encoded k-mers of the given ksize, plus one reserved bit so that invalid
// k-mers can be reported.  Constant member kmer_codec<kmer_t>::max_ksize
// gives the maximum ksize.  The constructor will throw an exception if ksize
// is larger.
//
// The encode(string) member encodes a string of DNA to a vector of kmer
// values.  Any kmer in string that contains an invalid base (not acgtACGT)
// is reported as value kmer_codec<kmer_t>::invalid_kmer.
//
// The decode(kmer_t) member decodes a kmer to a string of DNA.  It decodes
// the invalid_kmer to the string "invalid".
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

        kmer_t max_kmer() const { return max_kmer_; }

        kmer_t encode_base(char) const;
        kmer_t encode_kmer(std::string::const_iterator) const;

        void encode(const std::string&, kmer_t*) const;
        void encode(std::string&&, kmer_t*) const;

        std::vector<kmer_t> encode(const std::string&) const;
        std::vector<kmer_t> encode(std::string&&) const;

        std::string decode(kmer_t, bool rc = false) const;
};


// implementation ------------------------------------------------------------

template <typename kmer_t>
kmer_codec<kmer_t>::kmer_codec(int ksize, bool sstrand)
    : ksize_(ksize), sstrand_(sstrand), max_kmer_(0)
{
    //std::cerr << "kmer_codec<" << (sizeof(kmer_t)*8) << ">(" << ksize << "," << sstrand << ")" << std::endl;

    if (ksize < 1)
        raise_error("invalid k-mer size: %d", ksize);

    if (ksize > max_ksize)
        raise_error("k-mer size %d too large for datatype (%d)", ksize, max_ksize);

    if (!sstrand && (ksize & 0x1) != 0x1)
        raise_error("k-mer size must be odd for canonical encoding");

    max_kmer_ = (static_cast<kmer_t>(1) << (2*ksize-(sstrand?0:1))) - 1;
}

template <typename kmer_t>
kmer_t
kmer_codec<kmer_t>::encode_base(char c) const
{
    constexpr static kmer_t X = kmer_codec<kmer_t>::invalid_kmer;
    constexpr static kmer_t base_values[] = { 0, X, 1, X, X, X, 2, X, X, X, X, X, X, X, X, X, X, X, X, 3 };

    unsigned o = c - 'A';

    if (o > 19) {      // outside 'A'..'T'
        o = c - 'a';            // try 'a'..'t' smallcaps
        if (o > 19)
            return X;
    }

    return base_values[o];
}

template <typename kmer_t>
kmer_t
kmer_codec<kmer_t>::encode_kmer(const std::string::const_iterator pcur) const
{
    kmer_t res = 0;

    std::string::const_iterator pend = pcur + ksize_;

    if (sstrand_) {     // single strand: encode every base as two bits
        std::string::const_iterator p = pcur - 1;

        while (++p != pend)
            res = (res<<2) | encode_base(*p);
    }
    else {              // canonical: middle base as one bit, determines direction
        std::string::const_iterator pmid = pcur + (ksize_ / 2);

        if (!(encode_base(*pmid) & 2))  {     // middle base is a or c, encode forward
            std::string::const_iterator p = pcur - 1;

            while (++p != pmid)
                res = (res<<2) | encode_base(*p);

            res = (res<<1) | encode_base(*p); // central base encoded as 1 bit: a->0, c->1

            while (++p != pend)
                res = (res<<2) | encode_base(*p);
        }
        else {                                // middle base is g or t, encode reverse
            std::string::const_iterator p = pend;

            while (--p != pmid)
                res = (res<<2) | (encode_base(*p) ^ 3); // xor with 3 is complementary base

            res = (res<<1) | (encode_base(*p) ^ 3); // t->a->0, g->c->1

            while (p-- != pcur)
                res = (res<<2) | (encode_base(*p) ^ 3);
        }
    }

    // If any base was invalid, then the bitwise or of its all-ones representation
    // will have set the high bit in res, and we know we must return invalid_kmer.

    return res > max_kmer_ ? invalid_kmer : res;
}

template <typename kmer_t>
void
kmer_codec<kmer_t>::encode(const std::string& s, kmer_t* t) const
{
    std::string::const_iterator pcur = s.cbegin();
    std::string::const_iterator pend = pcur + s.size() + 1 - ksize_;

    while (pcur < pend)
        *t++ = encode_kmer(pcur++);
}

template <typename kmer_t>
void
kmer_codec<kmer_t>::encode(std::string &&s, kmer_t* t) const
{
    return encode((const std::string&)s, t);
}

template <typename kmer_t>
std::vector<kmer_t>
kmer_codec<kmer_t>::encode(const std::string& s) const
{
    std::vector<kmer_t> result;

    // a string of size n produces n+1-k kmers
    size_t n = s.size() + 1;

    if (static_cast<size_t>(ksize_) < n)
    {
        n -= ksize_;

        result.reserve(n);

        std::string::const_iterator pcur = s.cbegin();
        std::string::const_iterator pend = pcur + n;

        while (pcur != pend)
            result.push_back(encode_kmer(pcur++));
    }

    return result;
}

template <typename kmer_t>
std::vector<kmer_t>
kmer_codec<kmer_t>::encode(std::string &&s) const
{
    return encode(s);
}

template <typename kmer_t>
std::string
kmer_codec<kmer_t>::decode(kmer_t kmer, bool rc) const
{
    static const char CHARS[4] =  { 'a', 'c', 'g', 't' };
    static const char RCHARS[4] = { 't', 'g', 'c', 'a' };

    std::string result;

    if (kmer > max_kmer_) {
        result = "invalid";
    }
    else {
        result.resize(ksize_);

        kmer_t mer = kmer;

        if (!rc) {  // forward

            int n = ksize_;

            if (sstrand_) {
                while (n != 0) {
                    result[--n] = CHARS[0x3 & mer];
                    mer >>= 2;
                }
            }
            else {
                int m = ksize_/2 + 1;

                while (n != m) {
                    result[--n] = CHARS[0x3 & mer];
                    mer >>= 2;
                }

                result[--n] = CHARS[0x1 & mer];
                mer >>= 1;

                while (n != 0) {
                    result[--n] = CHARS[0x3 & mer];
                    mer >>= 2;
                }
            }
        }
        else {      // reverse complement
            int n = 0;

            if (sstrand_) {
                while (n != ksize_) {
                    result[n++] = RCHARS[0x3 & mer];
                    mer >>= 2;
                }
            }
            else {
                int m = ksize_/2;

                while (n != m) {
                    result[n++] = RCHARS[0x3 & mer];
                    mer >>= 2;
                }

                result[n++] = RCHARS[0x1 & mer];
                mer >>= 1;

                while (n != ksize_) {
                    result[n++] = RCHARS[0x3 & mer];
                    mer >>= 2;
                }
            }
        }
    }

    return result;
}


} // namespace kfc

#endif // kmercodec_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
