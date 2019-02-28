/* tallyman.cpp
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

#include <cstring>
#include "tallyman.h"

#include "utils.h"

namespace kfc {


std::unique_ptr<tallyman>
tallyman::create(int nbits, int max_gb)
{
    tallyman *ret;

    const int max_bits32 = 8*sizeof(val32_t);
    const int max_bits64 = 8*sizeof(val64_t);

    if (nbits < 1)
        raise_error("number of bits must be one or more: %d", nbits);

    if (nbits > max_bits64)
        raise_error("%d bits requested exceeds limit %d", nbits, max_bits64);

    std::uintmax_t max_mb = static_cast<std::uintmax_t>(max_gb) << 10;
    std::uintmax_t vec_mcount = static_cast<std::uintmax_t>(1) << (nbits > 20 ? nbits - 20 : 0);
    std::uintmax_t vec_mb = vec_mcount * sizeof(count_t);
    
    if (max_gb == 0)
    {
        unsigned long phy_mb = get_system_memory() >> 20;
        max_mb = phy_mb > 2048 ? phy_mb - 2048 : phy_mb;

        verbose_emit("defaulting max memory to all%s physical memory: %luG",
                phy_mb > 2048 ? " but 2G" : "",
                static_cast<unsigned long>(max_mb >> 10));
    }

    if (vec_mb > max_mb)
    {
        verbose_emit("vector memory (%luG) would exceed %luG: tallying using map",
                static_cast<unsigned long>(vec_mb >> 10),
                static_cast<unsigned long>(max_mb >> 10));

        if (nbits <= max_bits32)
            ret = new tallyman_map<val32_t>(nbits);
        else
            ret = new tallyman_map<val64_t>(nbits);
    }
    else
    {
        verbose_emit("vector memory (%luG) fits %luG: tallying using linear vector",
                static_cast<unsigned long>(vec_mb >> 10),
                static_cast<unsigned long>(max_mb >> 10));

        ret = new tallyman_vec(nbits);
    }
 
    return std::unique_ptr<tallyman>(ret);
}

// tallyman -----------------------------------------------------------------

void
tallyman::tally(tallyman::val32_t)
{
    raise_error("64-bit tallyman invoked with 32-bit argument");
}

void
tallyman::tally(tallyman::val64_t)
{
    raise_error("32-bit tallyman invoked with 64-bit argument");
}

void 
tallyman::tally(std::vector<tallyman::val32_t>)
{
    raise_error("64-bit tallyman invoked with vector of 32-bit arguments");
}

void
tallyman::tally(std::vector<tallyman::val64_t>)
{
    raise_error("32-bit tallyman invoked with vector of 64-bit arguments");
}

bool
tallyman::is_vec() const
{
    return false;
}

bool
tallyman::is_map32() const
{
    return false;
}

const std::unique_ptr<kfc::count_t[]>&
tallyman::get_results_vec() const
{
    static std::unique_ptr<kfc::count_t[]> dummy;
    raise_error("get_results_vec invoked on non-vector implementation ");
    return dummy;
}


const std::map<tallyman::val32_t,kfc::count_t>&
tallyman::get_results_map32() const
{
    static std::map<val32_t,count_t> dummy;
    raise_error("get_results_map32 invoked on non-map32 implementation ");
    return dummy;
}

const std::map<tallyman::val64_t,kfc::count_t>&
tallyman::get_results_map64() const
{
    static std::map<val64_t,count_t> dummy;
    raise_error("get_results_map64 invoked on non-map64 implementation ");
    return dummy;
}

// tallyman_vec -------------------------------------------------------------

tallyman_vec::tallyman_vec(int nbits)
    : tallyman(nbits), vec_(nullptr)
{
    const int max_bits = 8*sizeof(val64_t);
    if (nbits > max_bits)
        raise_error("number of bits (%d) exceeds maximum %d", nbits, max_bits);

    size_t alloc_n = max_value_ + 1;
    size_t alloc_size = alloc_n * sizeof(count_t);

    count_t *data = new count_t[alloc_n];
    if (!data)
        raise_error("failed to allocate memory (%luMB) for tally vector",
                static_cast<unsigned long>(alloc_size >> 20));

    std::memset(data, 0, alloc_size);

    vec_.reset(data);
}

void
tallyman_vec::tally(tallyman::val32_t i)
{
    if (i > max_value_)
        ++n_invalid_;
    else
        ++vec_[i];
}

void
tallyman_vec::tally(tallyman::val64_t i)
{
    if (i > max_value_)
        ++n_invalid_;
    else
        ++vec_[i];
}

void
tallyman_vec::tally(std::vector<tallyman::val32_t> ii)
{
    for (auto i : ii)
        if (i > max_value_)
            ++n_invalid_;
        else
            ++vec_[i];
}

void
tallyman_vec::tally(std::vector<tallyman::val64_t> ii)
{
    for (auto i : ii)
        if (i > max_value_)
            ++n_invalid_;
        else
            ++vec_[i];
}

bool
tallyman_vec::is_vec() const
{
    return true;
}

const std::unique_ptr<kfc::count_t[]>&
tallyman_vec::get_results_vec() const
{
    return vec_;
}

// mapXX --------------------------------------------------------------------

template<typename value_t>
tallyman_map<value_t>::tallyman_map(int nbits)
    : tallyman(nbits)
{
    const int max_bits = 8*sizeof(value_t);
    if (nbits > max_bits)
        raise_error("number of bits (%d) exceeds maximum %d", nbits, max_bits);
}


template<typename value_t>
void
tallyman_map<value_t>::tally(value_t i)
{
    if (i > max_value_)
        ++n_invalid_;
    else {
        iterator p = map_.lower_bound(i);

        if (p == map_.end() || i != p->first)
            map_.insert(p, std::make_pair(i,1));
        else
            ++(p->second);
    }
}

template<typename value_t>
void
tallyman_map<value_t>::tally(std::vector<value_t> ii)
{
    for (value_t i : ii)
        tally(i);
}

// map32 --------------------------------------------------------------------

template<>
bool
tallyman_map<tallyman::val32_t>::is_map32() const
{
    return true;
}

template<>
const std::map<tallyman::val32_t,kfc::count_t>&
tallyman_map<tallyman::val32_t>::get_results_map32() const
{
    return map_;
}

template<>
const std::map<tallyman::val64_t,kfc::count_t>&
tallyman_map<tallyman::val32_t>::get_results_map64() const
{
    static std::map<val64_t,count_t> dummy;
    raise_error("get_results_map64 invoked on non-map64 implementation ");
    return dummy;
}

// map64 --------------------------------------------------------------------

template<>
bool
tallyman_map<tallyman::val64_t>::is_map32() const
{
    return false;
}

template<>
const std::map<tallyman::val32_t,kfc::count_t>&
tallyman_map<tallyman::val64_t>::get_results_map32() const
{
    static std::map<val32_t,count_t> dummy;
    raise_error("get_results_map32 invoked on non-map32 implementation ");
    return dummy;
}

template<>
const std::map<tallyman::val64_t,kfc::count_t>&
tallyman_map<tallyman::val64_t>::get_results_map64() const
{
    return map_;
}


} // namespace kfc

// vim: sts=4:sw=4:ai:si:et
