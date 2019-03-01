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


// kmer_codec - encodes and decodes between dna and k-mers
//
template <typename kmer_t>
class kmer_codec {

    public:
        constexpr static kmer_t invalid_kmer = 4*sizeof(kmer_t);
       // static_cast<kmer_t>(-1);
        constexpr static int max_ksize = 4*sizeof(kmer_t) - 1;

    private:
        int ksize_;
        bool sstrand_;

    public:
        explicit kmer_codec(int ksize, bool sstrand = false);

        std::vector<kmer_t> encode(const std::string&);
        std::string decode(kmer_t);
};


// implementation ------------------------------------------------------------

template <typename kmer_t>
kmer_codec<kmer_t>::kmer_codec(int ksize, bool sstrand)
    : ksize_(ksize), sstrand_(sstrand)
{
    if (ksize < 1)
        raise_error("invalid k-mer size: %d", ksize);

    if (ksize > max_ksize)
        raise_error("k-mer size %d too large for datatype (%d)", ksize, max_ksize);

    if (!sstrand && (ksize & 0x1) != 0x1)
        raise_error("k-mer size must be odd for canonical encoding");
}

template <typename kmer_t>
std::vector<kmer_t>
kmer_codec<kmer_t>::encode(const std::string&)
{
    std::vector<kmer_t> result;
    raise_error("TODO: implement kmercodec::encode");
    return result;
}

template <typename kmer_t>
std::string
kmer_codec<kmer_t>::decode(kmer_t kmer)
{
    static const char CHARS[4] = { 'a', 'c', 'g', 't' };

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

// vim: sts=4:sw=4:ai:si:et
#endif // kmercodec_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
