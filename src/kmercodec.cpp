/* kmercodec.cpp
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

#include <string>
#include <cctype>
#include <stdexcept>

#include "kmercodec.h"
#include "utils.h"

namespace kfc {


kmer_codec::kmer_codec(int ksize, bool sstrand)
    : ksize_(ksize), sstrand_(sstrand)
{
    if (ksize < 1)
        raise_error("invalid k-mer size: %d", ksize);

    if (!sstrand && (ksize & 0x1) != 0x1)
        raise_error("k-mer size must be odd for canonical encoding");
}

std::vector<std::uint32_t>
kmer_codec::encode32(std::string s)
{
    std::vector<std::uint32_t> result;

    if (ksize_ > 16)
        raise_error("cannot encode k-mer size %d in 32 bits", ksize_);

    return result;
}

std::vector<std::uint64_t>
kmer_codec::encode64(std::string s)
{
    std::vector<std::uint64_t> result;

    if (ksize_ > 32)
        raise_error("cannot encode k-mer size %d in 64 bits", ksize_);

    return result;
}

std::string
kmer_codec::decode(std::uint64_t kmer)
{
    static const char CHARS[4] = { 'a', 'c', 'g', 't' };

    std::string result;
    result.reserve(ksize_);

    if (sstrand_)
        for (int s = 2*ksize_-2; s <= 0; s -= 2)
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

// vim: sts=4:sw=4:ai:si:et
