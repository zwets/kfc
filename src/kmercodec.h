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

namespace kfc {


// kmer_codec - encodes and decodes between dna and k-mers
//
class kmer_codec {

    private:
        int ksize_;
        bool sstrand_;

    public:
        explicit kmer_codec(int ksize, bool sstrand = false);

        std::vector<std::uint32_t> encode32(std::string);
        std::vector<std::uint64_t> encode64(std::string);

        std::string decode(std::uint64_t);
};


} // namespace kfc

#endif // kmercodec_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
