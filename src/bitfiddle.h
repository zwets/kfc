/* bitfiddle.h
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
#ifndef bitfiddle_h_INCLUDED
#define bitfiddle_h_INCLUDED

namespace kfc {


// bitsize - constant number of bits in T
//
template <typename T>
constexpr unsigned bitsize = 8 * sizeof(T);

// low_bits - constant integral T with the N low bits set (N < bits in T)
//
template <typename T, unsigned N>
constexpr T low_bits = (((T)1)<<N)-1;

// high_bit - constant integral T just the high bit set
//
template <typename T>
constexpr T high_bit = ((T)1) << (bitsize<T>-1);

// high_bits - constant integral T with the N high bits set (N>0)
//
template <typename T, unsigned N>
constexpr T high_bits = typename std::make_signed<T>::type(high_bit<T>) >> (N-1); 

// signed_shr - perform signed right shift, even when T is unsigned
//
template <typename T>
constexpr T signed_shr(T t, unsigned n) {
    return typename std::make_signed<T>::type(t) >> n;
}

// flush_hibit - flush the value with the value of the high bit
//
template <typename T>
constexpr T flush_hibit(T t) {
    return signed_shr(t, bitsize<T>-1);
}


} // namespace kfc

#endif // bitfiddle_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
