/* kmerencoder.cpp
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

#include <cstdlib>
#include "kmerencoder.h"
#include "kmercodec.h"

namespace kfc {
/*
typedef void (*encode_fn32)(const char*, const char*, std::uint32_t*);
typedef void (*encode_fn64)(const char*, const char*, std::uint64_t*);
typedef std::uint32_t (*decode_fn32)(std::uint32_t);
typedef std::uint64_t (*decode_fn64)(std::uint64_t);

struct fn32 {
    encode_fn32 encode;
    decode_fn32 decode;
};

struct fn64 {
    encode_fn64 encode;
    decode_fn64 decode;
};

static fn32 ss_fn32[] = {
    { 0, 0 }, 
    { ss_encode<std::uint32_t,1>, ss_decode<std::uint32_t,1> },
    { ss_encode<std::uint32_t,2>, ss_decode<std::uint32_t,2> },
    { ss_encode<std::uint32_t,3>, ss_decode<std::uint32_t,3> },
    { ss_encode<std::uint32_t,4>, ss_decode<std::uint32_t,4> },
    { ss_encode<std::uint32_t,5>, ss_decode<std::uint32_t,5> },
    { ss_encode<std::uint32_t,6>, ss_decode<std::uint32_t,6> },
    { ss_encode<std::uint32_t,7>, ss_decode<std::uint32_t,7> },
    { ss_encode<std::uint32_t,8>, ss_decode<std::uint32_t,8> },
    { ss_encode<std::uint32_t,9>, ss_decode<std::uint32_t,9> },
    { ss_encode<std::uint32_t,10>, ss_decode<std::uint32_t,10> },
    { ss_encode<std::uint32_t,11>, ss_decode<std::uint32_t,11> },
    { ss_encode<std::uint32_t,12>, ss_decode<std::uint32_t,12> },
    { ss_encode<std::uint32_t,13>, ss_decode<std::uint32_t,13> },
    { ss_encode<std::uint32_t,13>, ss_decode<std::uint32_t,13> },
    { ss_encode<std::uint32_t,14>, ss_decode<std::uint32_t,14> },
    { ss_encode<std::uint32_t,15>, ss_decode<std::uint32_t,15> },
}

static fn32 ds_fn32 = {
}

static fn64 ss_fn64 = {
}

static fn64 ds_fn64 = {
}
*/
} // namespace kfc

// vim: sts=4:sw=4:ai:si:et
