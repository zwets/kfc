/* kmercounter.cpp
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

#include "kmercounter.h"


namespace kfc {


kmer_counter::kmer_counter(int ksize, bool s_strand, int mem_gb, int n_threads)
    : ksize_(ksize), s_strand_(s_strand), n_threads_(n_threads)
{
    tallyman_ = tallyman::create(2 * ksize + (s_strand ? 1 : 0), mem_gb);
}

//void
//kmer_counter::process(std::string data)
//{
//}


} // namespace kfc

// vim: sts=4:sw=4:ai:si:et
