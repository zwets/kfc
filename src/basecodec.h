/* basecodec.h
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
#ifndef basecodec_h_INCLUDED
#define basecodec_h_INCLUDED

//
// basecodec.h - functions for encoding and decoding bases
//


namespace kfc {


// encode_base - encode base character to kmer_t value
//
// template parameters
// - kmer_t - the type to encode into
// - X      - the value to return for an invalid base
//
template <typename kmer_t, kmer_t X>
inline kmer_t encode_base(unsigned char c)
{
    static_assert(sizeof(c) == 1, "out of bounds risk: unsigned char not one byte");

    constexpr static kmer_t base_values[] = { 
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0x00
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0x10
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0x20
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0x30
        X, 0, X, 1, X, X, X, 2, X, X, X, X, X, X, X, X,  // 0x40
        X, X, X, X, 3, X, X, X, X, X, X, X, X, X, X, X,  // 0x50
        X, 0, X, 1, X, X, X, 2, X, X, X, X, X, X, X, X,  // 0x60
        X, X, X, X, 3, X, X, X, X, X, X, X, X, X, X, X,  // 0x70
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0x80
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0x90
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0xA0
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0xB0
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0xC0
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0xD0
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X,  // 0xE0
        X, X, X, X, X, X, X, X, X, X, X, X, X, X, X, X   // 0xF0
    };

    return base_values[c];
}


// decode_base - decode bottom two bits of kmer to base
//
template <typename kmer_t>
inline char decode_base(kmer_t kmer)
{
    constexpr static char kmer_bases[] = { 'a', 'c', 'g', 't' };
    return kmer_bases[ kmer & 0x3 ];
}


// decode_comp_base - decode bottom two bits of kmer to complement of base
//
template <typename kmer_t>
inline char decode_comp_base(kmer_t kmer)
{
    constexpr static char kmer_bases[] = { 't', 'g', 'c', 'a' };
    return kmer_bases[ kmer & 0x3 ];
}


} // namespace kfc

#endif // basecodec_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
