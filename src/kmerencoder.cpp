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
#include <string>
#include "kmerencoder.h"
#include "kmercodec.h"

namespace kfc {

typedef void (*encode_fn32)(const char*, const char*, std::uint32_t*);
typedef void (*encode_fn64)(const char*, const char*, std::uint64_t*);

typedef std::string (*decode_fn32)(std::uint32_t, bool);
typedef std::string (*decode_fn64)(std::uint64_t, bool);

// encode functions ----------------------------------------------------

static encode_fn32 ss_enc32[16] = {
    0,
    ss_encode<std::uint32_t,1>,
    ss_encode<std::uint32_t,2>,
    ss_encode<std::uint32_t,3>,
    ss_encode<std::uint32_t,4>,
    ss_encode<std::uint32_t,5>,
    ss_encode<std::uint32_t,6>,
    ss_encode<std::uint32_t,7>,
    ss_encode<std::uint32_t,8>,
    ss_encode<std::uint32_t,9>,
    ss_encode<std::uint32_t,10>,
    ss_encode<std::uint32_t,11>,
    ss_encode<std::uint32_t,12>,
    ss_encode<std::uint32_t,13>,
    ss_encode<std::uint32_t,14>,
    ss_encode<std::uint32_t,15>
};

static encode_fn32 ds_enc32[16] = {
    0, ds_encode<std::uint32_t,1>,
    0, ds_encode<std::uint32_t,3>,
    0, ds_encode<std::uint32_t,5>,
    0, ds_encode<std::uint32_t,7>,
    0, ds_encode<std::uint32_t,9>,
    0, ds_encode<std::uint32_t,11>,
    0, ds_encode<std::uint32_t,13>,
    0, ds_encode<std::uint32_t,15>
};

static encode_fn64 ss_enc64[32] = {
    0,
    ss_encode<std::uint64_t,1>,
    ss_encode<std::uint64_t,2>,
    ss_encode<std::uint64_t,3>,
    ss_encode<std::uint64_t,4>,
    ss_encode<std::uint64_t,5>,
    ss_encode<std::uint64_t,6>,
    ss_encode<std::uint64_t,7>,
    ss_encode<std::uint64_t,8>,
    ss_encode<std::uint64_t,9>,
    ss_encode<std::uint64_t,10>,
    ss_encode<std::uint64_t,11>,
    ss_encode<std::uint64_t,12>,
    ss_encode<std::uint64_t,13>,
    ss_encode<std::uint64_t,14>,
    ss_encode<std::uint64_t,15>,
    ss_encode<std::uint64_t,16>,
    ss_encode<std::uint64_t,17>,
    ss_encode<std::uint64_t,18>,
    ss_encode<std::uint64_t,19>,
    ss_encode<std::uint64_t,20>,
    ss_encode<std::uint64_t,21>,
    ss_encode<std::uint64_t,22>,
    ss_encode<std::uint64_t,23>,
    ss_encode<std::uint64_t,24>,
    ss_encode<std::uint64_t,25>,
    ss_encode<std::uint64_t,26>,
    ss_encode<std::uint64_t,27>,
    ss_encode<std::uint64_t,28>,
    ss_encode<std::uint64_t,29>,
    ss_encode<std::uint64_t,30>,
    ss_encode<std::uint64_t,31>
};

static encode_fn64 ds_enc64[32] = {
    0, ds_encode<std::uint64_t,1>,
    0, ds_encode<std::uint64_t,3>,
    0, ds_encode<std::uint64_t,5>,
    0, ds_encode<std::uint64_t,7>,
    0, ds_encode<std::uint64_t,9>,
    0, ds_encode<std::uint64_t,11>,
    0, ds_encode<std::uint64_t,13>,
    0, ds_encode<std::uint64_t,15>,
    0, ds_encode<std::uint64_t,17>,
    0, ds_encode<std::uint64_t,19>,
    0, ds_encode<std::uint64_t,21>,
    0, ds_encode<std::uint64_t,23>,
    0, ds_encode<std::uint64_t,25>,
    0, ds_encode<std::uint64_t,27>,
    0, ds_encode<std::uint64_t,29>,
    0, ds_encode<std::uint64_t,31>
};

// decode functions ---------------------------------------------------

static decode_fn32 ss_dec32[16] = {
    0,
    ss_decode<std::uint32_t,1>,
    ss_decode<std::uint32_t,2>,
    ss_decode<std::uint32_t,3>,
    ss_decode<std::uint32_t,4>,
    ss_decode<std::uint32_t,5>,
    ss_decode<std::uint32_t,6>,
    ss_decode<std::uint32_t,7>,
    ss_decode<std::uint32_t,8>,
    ss_decode<std::uint32_t,9>,
    ss_decode<std::uint32_t,10>,
    ss_decode<std::uint32_t,11>,
    ss_decode<std::uint32_t,12>,
    ss_decode<std::uint32_t,13>,
    ss_decode<std::uint32_t,14>,
    ss_decode<std::uint32_t,15>
};

static decode_fn32 ds_dec32[16] = {
    0, ds_decode<std::uint32_t,1>,
    0, ds_decode<std::uint32_t,3>,
    0, ds_decode<std::uint32_t,5>,
    0, ds_decode<std::uint32_t,7>,
    0, ds_decode<std::uint32_t,9>,
    0, ds_decode<std::uint32_t,11>,
    0, ds_decode<std::uint32_t,13>,
    0, ds_decode<std::uint32_t,15>
};

static decode_fn64 ss_dec64[32] = {
    0,
    ss_decode<std::uint64_t,1>,
    ss_decode<std::uint64_t,2>,
    ss_decode<std::uint64_t,3>,
    ss_decode<std::uint64_t,4>,
    ss_decode<std::uint64_t,5>,
    ss_decode<std::uint64_t,6>,
    ss_decode<std::uint64_t,7>,
    ss_decode<std::uint64_t,8>,
    ss_decode<std::uint64_t,9>,
    ss_decode<std::uint64_t,10>,
    ss_decode<std::uint64_t,11>,
    ss_decode<std::uint64_t,12>,
    ss_decode<std::uint64_t,13>,
    ss_decode<std::uint64_t,14>,
    ss_decode<std::uint64_t,15>,
    ss_decode<std::uint64_t,16>,
    ss_decode<std::uint64_t,17>,
    ss_decode<std::uint64_t,18>,
    ss_decode<std::uint64_t,19>,
    ss_decode<std::uint64_t,20>,
    ss_decode<std::uint64_t,21>,
    ss_decode<std::uint64_t,22>,
    ss_decode<std::uint64_t,23>,
    ss_decode<std::uint64_t,24>,
    ss_decode<std::uint64_t,25>,
    ss_decode<std::uint64_t,26>,
    ss_decode<std::uint64_t,27>,
    ss_decode<std::uint64_t,28>,
    ss_decode<std::uint64_t,29>,
    ss_decode<std::uint64_t,30>,
    ss_decode<std::uint64_t,31>
};

static decode_fn64 ds_dec64[32] = {
    0, ds_decode<std::uint64_t,1>,
    0, ds_decode<std::uint64_t,3>,
    0, ds_decode<std::uint64_t,5>,
    0, ds_decode<std::uint64_t,7>,
    0, ds_decode<std::uint64_t,9>,
    0, ds_decode<std::uint64_t,11>,
    0, ds_decode<std::uint64_t,13>,
    0, ds_decode<std::uint64_t,15>,
    0, ds_decode<std::uint64_t,17>,
    0, ds_decode<std::uint64_t,19>,
    0, ds_decode<std::uint64_t,21>,
    0, ds_decode<std::uint64_t,23>,
    0, ds_decode<std::uint64_t,25>,
    0, ds_decode<std::uint64_t,27>,
    0, ds_decode<std::uint64_t,29>,
    0, ds_decode<std::uint64_t,31>
};

// define the init_ member ---------------------------------------------

template<>
void
kmer_encoder<std::uint32_t>::init_(unsigned ksize, bool sstrand)
{
    encode_ = sstrand ? ss_enc32[ksize] : ds_enc32[ksize];
    decode_ = sstrand ? ss_dec32[ksize] : ds_dec32[ksize];
}

template<>
void
kmer_encoder<std::uint64_t>::init_(unsigned ksize, bool sstrand)
{
    encode_ = sstrand ? ss_enc64[ksize] : ds_enc64[ksize];
    decode_ = sstrand ? ss_dec64[ksize] : ds_dec64[ksize];
}


} // namespace kfc

// vim: sts=4:sw=4:ai:si:et
