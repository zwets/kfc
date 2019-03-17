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

#include "basecodec.h"
#include "bitfiddle.h"

namespace kfc {


// ss_encode_one - single-strand encode single kmer
//
// - p must point at the array of ksize bases to be encoded in the kmer
// - returns the encoded kmer or an arbitrary value with high bit set,
//   when there was an invalid base in the kmer
//
template <typename kmer_t, unsigned ksize>
kmer_t
ss_encode_one(const char *p)
{
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");
    static_assert(0 < ksize && ksize < 4*sizeof(kmer_t),
            "template argument ksize must be in range [1,bitsize/2)");

    constexpr static kmer_t invalid_value = high_bits<kmer_t,2*ksize>;

    kmer_t kmer = 0;

    for (unsigned i = 0; i < ksize; ++i) {
        kmer <<= 2;
        // add the new base or mask high bits without risk of shifting out
        kmer |= encode_base<kmer_t,invalid_value>(p[i]);
    }
    
    return kmer;
}


// ss_encode - single-strand encode dna to sequence of kmers
//
// - p0 must point at the first base of n>0 kmers to encode
// - p1 must point one beyond end of string, and at least at p0+ksize
// - t must point at an array of n (= p1-p0-ksize+1) kmer_t, to receive
//   encoded valid kmers, or values with high bit set to mark invalid
//
template <typename kmer_t, unsigned ksize>
void
ss_encode(const char *p0, const char *p1, kmer_t *t)
{
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");
    static_assert(0 < ksize && ksize < 4*sizeof(kmer_t),
            "template argument ksize must be in range [1,bitsize-2]");

    constexpr static kmer_t invalid_value = high_bits<kmer_t,2*ksize>;

    kmer_t kmer = 0;

    const char *p = p0;

    // we need to fill the kmer before we can roll new bases into it

    for (unsigned i = 0; i < ksize; ++i) {
        kmer_t new_base = encode_base<kmer_t,invalid_value>(*p++);
        // make place for the new base
        kmer <<= 2;
        // clear all when the new base has high bit set (invalid)
        kmer &= ~flood_hibit<kmer_t>(new_base);
        // or in the new base or set high bits for next ksize kmers
        kmer |= new_base;
    }
    
    *t++ = kmer;

    // p now points at the next base to roll into the kmer

    while (p != p1) {
        // clear bits of the outgoing base if current kmer is good
        kmer &= signed_shr<kmer_t>(kmer|~high_bit<kmer_t>, bitsize<kmer_t>-2*ksize+1);
        // make place for the new base
        kmer <<= 2;
        // retrieve the new base
        kmer_t new_base = encode_base<kmer_t,invalid_value>(*p++);
        // clear all when new base has high bit set (invalid kmer)
        kmer &= ~flood_hibit<kmer_t>(new_base);
        // or in the new base or invalid value with ksize high bits
        kmer |= new_base;
        // write it
        *t++ = kmer;
    }
}


// ss_revcomp - reverse complement single stranded kmer
//
// - returns the reverse complement of single-strand encoded kmer,
//   or an arbitrary number with high bit set if the input had it too
//
template <typename kmer_t, unsigned ksize>
kmer_t
ss_revcomp(kmer_t in)
{
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");
    static_assert(0 < ksize && ksize < 4*sizeof(kmer_t),
            "template argument ksize must be in range [1,bitsize/2)");
    static_assert(ksize & 1,
            "template argument ksize must be odd for double strand encoding");

    kmer_t out = 0;

    in ^= low_bits<kmer_t,2*ksize>;

    // reverse by shifting bitpairs from in to out
    for (unsigned i = 0; i < ksize; ++i) {
        out <<= 2;
        out |= (in & 0x3);
        in = signed_shr(in, 2);  // need to retain the high bit
    }

    // return result but set its error bit if it was set on in 
    return out | (in & high_bit<kmer_t>);
}


// ss_decode - decode kmer to dna string
//
template <typename kmer_t, unsigned ksize>
std::string
ss_decode(kmer_t kmer, bool rc = false)
{
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");
    static_assert(0 < ksize && ksize < 4*sizeof(kmer_t),
            "template argument ksize must be in range [1,bitsize-2]");

    std::string result;
    result.resize(ksize, 'X');

    if (!(kmer & high_bit<kmer_t>)) {
        if (rc)
            for (unsigned n = 0; n != ksize; ++n) {
                result[n] = decode_comp_base(kmer);
                kmer >>= 2;
            }
        else
            for (unsigned n = ksize; n != 0; ) {
                result[--n] = decode_base(kmer);
                kmer >>= 2;
            }
    }

    return result;
}


// ss_to_ds - convert ss-encoded kmer to ds-encoded kmer,
//            reverse complementing if needed
//
template <typename kmer_t, unsigned ksize>
kmer_t
ss_to_ds(kmer_t kmer)
{
    static_assert(ksize & 1,
            "template argument ksize must be odd to convert ss to ds");

    constexpr kmer_t half_mask = low_bits<kmer_t,ksize>;

    kmer_t out = ((kmer>>ksize) & 0x1) ? ss_revcomp<kmer_t,ksize>(kmer) : kmer;

    return ((out & (half_mask << ksize)) >> 1) | (out & (high_bit<kmer_t> | half_mask));
}


// --- ds encode ---------------------------------------------------------


// ds_encode_one - double-strand encode single kmer
//
// - p must point at the array of ksize bases to be encoded in the kmer
// - returns the encoded kmer or an arbitrary value with high bit set,
//   when there was an invalid base in the kmer
//
template <typename kmer_t, unsigned ksize>
kmer_t
ds_encode_one(const char *p0)
{
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");
    static_assert(0 < ksize && ksize < 4*sizeof(kmer_t),
            "template argument ksize must be in range [1,bitsize/2)");
    static_assert(ksize & 1,
            "template argument ksize must be odd for double strand encoding");

    constexpr static kmer_t invalid_value = high_bits<kmer_t,2*ksize>;

    kmer_t kmer = 0;

    // based on the middle base we encode forward or reverse kmer

    const char *pmid = p0 + (ksize / 2);
    kmer_t bmid = encode_base<kmer_t,invalid_value>(*pmid);

    if (!(bmid & 2)) { // middle base is a or c, encode forward
        const char *p = p0;

        while (p != pmid) {
            kmer <<= 2; 
            kmer |= encode_base<kmer_t,invalid_value>(*p++);
        }

        kmer <<= 1;
        kmer |= bmid;

        while (++p != p0 + ksize) {
            kmer <<= 2;
            kmer |= encode_base<kmer_t,invalid_value>(*p);
        }
    }
    else {  // middle base is g or t, encode reverse complement
        const char *q = p0 + ksize;

        while (--q != pmid) {
            kmer <<= 2;
            kmer |= encode_base<kmer_t, invalid_value>(*q);
        }

        kmer <<= 1;
        kmer |= (bmid & 1);

        while (q-- != p0) {
            kmer <<= 2;
            kmer |= encode_base<kmer_t, invalid_value>(*q);
        }

        // complement
        kmer ^= low_bits<kmer_t,2*ksize-1>;
    }

    return kmer;
}


// ds_encode - double-strand encode dna to sequence of kmers
//
// - p0 must point at the first base of n>0 kmers to encode
// - p1 must point one beyond end of string, and at least at p0+ksize
// - t must point at an array of n (= p1-p0-ksize+1) kmer_t, to receive
//   encoded valid kmers, or values with high bit set to mark invalid
//
template <typename kmer_t, unsigned ksize>
void
ds_encode(const char *p0, const char *p1, kmer_t *t)
{
    static_assert(std::is_unsigned<kmer_t>::value,
            "template argument kmer_t must be unsigned integral");

    static_assert(0 < ksize && ksize < 4*sizeof(kmer_t),
            "template argument ksize must be in range [1,kbits-1]");

    static_assert(ksize & 1,
            "template argument ksize must be odd");

#if 0
    // IMPLEMENTATION 1: use ss_encode, then ss_to_ds
    ss_encode<kmer_t,ksize>(p0, p1, t);

    for (kmer_t *pt = t; pt != t + (p1 - ksize + 1 - p0); ++pt)
        *pt = ss_to_ds<kmer_t,ksize>(*pt);
#else
    // IMPLEMENTATION 2: use encode_kmer_ds in turn
    const char *p = p0;
    while (p != p1 - ksize + 1)
        *t++ = ds_encode_one<kmer_t,ksize>(p++);
// #else
    // IMPLEMENTATION 3: write rolling encode
#endif
}


// ds_decode ----------------------------------------------------------
//
template <typename kmer_t, unsigned ksize>
std::string
ds_decode(kmer_t kmer, bool rc = false)
{
    std::string result;
    result.resize(ksize, 'X');

    if (!(kmer & high_bit<kmer_t>)) {

        if (!rc) {
            constexpr unsigned m = ksize/2 + 1;
            int n = ksize;

            while (n != m) {
                result[--n] = decode_base(kmer);
                kmer >>= 2;
            }

            result[--n] = decode_base(kmer & 1);
            kmer >>= 1;

            while (n != 0) {
                result[--n] = decode_base(kmer);
                kmer >>= 2;
            }
        }
        else {
            constexpr unsigned m = ksize/2;
            int n = 0;

            while (n != m) {
                result[n++] = decode_comp_base(kmer);
                kmer >>= 2;
            }

            result[n++] = decode_comp_base(kmer & 1);
            kmer >>= 1;

            while (n != ksize) {
                result[n++] = decode_comp_base(kmer);
                kmer >>= 2;
            }
        }
    }

    return result;
}


} // namespace kfc

#endif // kmercodec_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
