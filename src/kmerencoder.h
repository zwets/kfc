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
#include "utils.h"

namespace kfc {


// kmer_encoder - encodes and decodes between DNA and k-mers
//
// A kmer_encoder is constructed for a fixed kmer size and encoding method.
// The encoding method is canonical (default) or single stranded; see the
// README.md for explanation.
//
// Template parameter kmer_t must be an unsigned integral type that can hold
// encoded k-mers of the given ksize, plus one reserved bit so that invalid
// k-mers can be reported.  Constant member kmer_encoder<kmer_t>::max_ksize
// gives the maximum ksize.  The constructor will throw an exception if ksize
// is larger.
//
// The encode() members encode a string of DNA to a vector of kmer values.
// Any kmer in string that contains an invalid base (not acgtACGT) is reported
// as an arbitrary number with its high bit set.
//
// The decode(kmer_t) member decodes a kmer to a string of DNA.  It decodes
// invalid kmers to the string "invalid".
//
template <typename kmer_t>
class kmer_encoder {
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");

    public:
        constexpr static int max_ksize = 4*sizeof(kmer_t) - 1;

        constexpr static bool is_invalid(kmer_t kmer) { return kmer & high_bit; }

    private:
        const int ksize_;
        const bool sstrand_;
        const kmer_t max_kmer_;
        const kmer_t invalid_mask_;

    public:
        explicit kmer_encoder(int ksize, bool sstrand = false);

        kmer_t max_kmer() const { return max_kmer_; }

        kmer_t encode_kmer(std::string::const_iterator) const;

        void encode(const std::string&, kmer_t*) const;
        void encode(std::string&&, kmer_t*) const;

        std::vector<kmer_t> encode(const std::string&) const;
        std::vector<kmer_t> encode(std::string&&) const;

        std::string decode(kmer_t, bool rc = false) const;

    private:
        constexpr static kmer_t high_bit = ((kmer_t)1) << (8*sizeof(kmer_t)-1);
        constexpr static kmer_t all_ones = ~((kmer_t)0);
        constexpr static kmer_t low_bits(int n) { return ((((kmer_t)1)<<n)-1); }
        constexpr static kmer_t high_bits(int n) { return ((((kmer_t)1)<<n)-1) << (8*sizeof(kmer_t)-n); }

    private:
        void rolling_encode(std::string::const_iterator, std::string::const_iterator, kmer_t*) const;
        kmer_t chomp_invalids(std::string::const_iterator&, std::string::const_iterator, kmer_t*) const;
};


// implementation ------------------------------------------------------------

template <typename kmer_t>
kmer_encoder<kmer_t>::kmer_encoder(int ksize, bool sstrand)
    : ksize_(ksize), sstrand_(sstrand), 
      max_kmer_(low_bits(2*ksize-(sstrand?0:1))),
      invalid_mask_(high_bits(2*ksize))
{
    //std::cerr << "kmer_encoder<" << (sizeof(kmer_t)*8) << ">(" << ksize << "," << sstrand << ")" << std::endl;
    std::cerr << "ksize " << ksize_ << "; max_kmer " << max_kmer_ << "; invalid_mask " << invalid_mask_ << std::endl;

    if (ksize < 1)
        raise_error("invalid k-mer size: %d", ksize);

    if (ksize > max_ksize)
        raise_error("k-mer size %d too large for datatype (%d)", ksize, max_ksize);

    if (!sstrand && (ksize & 0x1) != 0x1)
        raise_error("k-mer size must be odd for canonical encoding");

}

template <typename kmer_t>
void
kmer_encoder<kmer_t>::encode(const std::string& s, kmer_t* t) const
{
    std::string::const_iterator pcur = s.cbegin();
    std::string::const_iterator pend = s.cend();

    if (sstrand_)
        encode_kmers_ss<kmer_t>(s.cbegin(), s.cend(), t);
    else
        encode_kmers_ds(s.cbegin(), s.cend(), t);

    if (sstrand_)
        rolling_encode(pcur, pend, t);
    else
        while (pcur < pend)
            *t++ = encode_kmer(pcur++);
}

template <typename kmer_t>
void
kmer_encoder<kmer_t>::encode(std::string &&s, kmer_t* t) const
{
    return encode((const std::string&)s, t);
}
/*
    kmer_t kmer = all_ones;

    while (kmer == invalid_kmer && p != pend) {

        std::string::const_iterator pchk = p + ksize_;

        while (encode_base(*--pchk) != invalid_kmer)
            if (pchk == p) {        // we made it without invalids, so p is good
                kmer = encode_kmer(p++);  // encode it, add it, and move p forward
                *t++ = kmer;
                break;              // so end with p at next starting point
            }

        if (kmer == invalid_kmer) { // pck has poison base, fill invalids up to it
            if (pchk >= pend)
                pchk = pend - 1;    // the last possible invalid is one before pend
            do {                    // add an invalid for p upto and including pchk
                *t++ = kmer;
            } while (p++ != pchk);  // and end with p++ just beyond pchk
        }
    }
*/
template <typename kmer_t>
void
kmer_encoder<kmer_t>::rolling_encode(std::string::const_iterator pbegin, std::string::const_iterator pend, kmer_t *t) const
    // Note: - pend is one past the last kmer starting position, not one past the end of string (which is kmer chars further)
    //       - in other words, we must encode exactly one kmer for each of [pbegin,pend)
    //       - this method is private, we may assume that it is called with pbegin < pend
{
    std::string::const_iterator p = pbegin;
    std::string::const_iterator pk = pbegin + ksize_;
   
    kmer_t kmer = 0;

    do { // fill the first
        kmer <<= 2;
        kmer |= encode_base(*p);
    } while (++p != pk);

    *t++ = kmer;    // but blank out when high bit set!
    //kmer | (signed)kmer>>all but one;

    while (p != pend) {
        
        kmer_t new_base = encode_base(p[ksize_-1]);

        if (is_invalid(new_base)) {
            kmer = -1 ; //chomp_invalids(p, pend, t);
        }
        else { // new base is good
            kmer <<= 2;        
            kmer |= new_base;
            kmer &= max_kmer_;
            *t++ = kmer;
            ++p;
        }
    }
}

template <typename kmer_t>
std::vector<kmer_t>
kmer_encoder<kmer_t>::encode(const std::string& s) const
{
    std::vector<kmer_t> result;

    // a string of length s produces n=s+1-k kmers
    size_t n = s.size() + 1;

    if (static_cast<size_t>(ksize_) < n)
    {
        n -= ksize_;

        result.reserve(n);

        std::string::const_iterator pcur = s.cbegin(); // first
        std::string::const_iterator pend = pcur + n;   // one past

        while (pcur != pend)
            result.push_back(encode_kmer(pcur++));
    }

    return result;
}

template <typename kmer_t>
std::vector<kmer_t>
kmer_encoder<kmer_t>::encode(std::string &&s) const
{
    return encode(s);
}

template <typename kmer_t>
std::string
kmer_encoder<kmer_t>::decode(kmer_t kmer, bool rc) const
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

#endif // kmerencoder_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
