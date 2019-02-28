/* kmercounter.h
 * 
 * Copyright (C) 2018  Marco van Zwetselaar <io@zwets.it>
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
#ifndef kmercounter_h_INCLUDED
#define kmercounter_h_INCLUDED

#include <memory>
#include <vector>
#include "tallyman.h"

namespace kfc {


// kmer_counter - counter for kmers of certain size and strandness
//
// Processes a series of strings, returns the counts of each k-mer seen.
//
class kmer_counter
{
    private:
        std::unique_ptr<tallyman> tallyman_;
        int ksize_;
        int s_strand_;
        int n_threads_;

    public:
        kmer_counter(int ksize, bool single_strand = false, int max_gb = 0, int n_threads = 0);
        void process(std::string data);
        void reset();
        const tallyman* get_tallyman() const { return tallyman_.get(); }
};


} // namespace kfc

#endif // kmercounter_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
