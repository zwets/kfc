/* tallyman.h
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
#ifndef tallyman_h_INCLUDED
#define tallyman_h_INCLUDED

#include <cctype>
#include <memory>
#include <vector>
#include <cstring>
#include <map>
#include "utils.h"

namespace kfc {


// tallyman - keeps counts of encoded kmers (or any integer type)
//
// Core operation is tally(n), which increments the count of integer n by 1,
// or increments the invalid_count if n is larger than nbits can contain.
//
// This class has four implementations with different space and time
// characteristics.  Given B the number of bits in the integers to be tallied,
// N the number of integers tallied, and C the size of count_t,
// - tallyman_vec uses a linear array with O(1) lookup and O(2^B) memory,
//   with C contributing a factor 2 if its type is 64 rather than 32-bit
// - tallyman_map uses a map of O(N) storage and O(log N) lookup, with C
//   contributing additional 25-50% if it is 64-bit
//
// The tallyman::create() factory method picks the fastest implementation
// that fits a given approximate memory limit.
//
// The get_results_X() members return the tallied counts.  For performance
// reasons, these members do not shield from the underlying implementation.
// Use the is_vec() and is_map() selectors to find out which get_results_X()
// is applicable to the actual class.
//
template <typename value_t, typename count_t>
class tallyman {

    protected:
        value_t max_value_;
        count_t n_invalid_;

    public:
        static tallyman<value_t,count_t>* create(int nbits, int max_gb = 0);

    protected:
        tallyman(int nbits);

    public:
        virtual ~tallyman() { }

	virtual void tally(const std::vector<value_t>&) = 0;

        virtual bool is_vec() const { return false; }
        virtual bool is_map() const { return false; }

        virtual const std::unique_ptr<count_t[]>& get_results_vec() const = 0;
        virtual const std::map<value_t,count_t>& get_results_map() const = 0;

        count_t invalid_count() { return n_invalid_; }
};

template <typename value_t, typename count_t>
class tallyman_vec : public tallyman<value_t,count_t>
{
    private:
        std::unique_ptr<count_t[]> vec_;

        void tally(value_t i);

    public:
        tallyman_vec<value_t,count_t>(int nbits);
        tallyman_vec<value_t,count_t>(const tallyman_vec&) = delete;
        tallyman_vec<value_t,count_t>& operator=(const tallyman_vec&) = delete;

	virtual void tally(const std::vector<value_t>& ii) { for (auto i : ii) tally(i); }

        virtual bool is_vec() const { return true; }

        virtual const std::unique_ptr<count_t[]>& get_results_vec() const { return vec_; }
        virtual const std::map<value_t,count_t>& get_results_map() const;
};

template<typename value_t, typename count_t>
class tallyman_map : public tallyman<value_t,count_t>
{
    private:
        typedef typename std::map<value_t,count_t>::iterator iterator;
        std::map<value_t,count_t> map_;

        void tally(value_t i);

    public:
        tallyman_map<value_t,count_t>(int nbits);
        tallyman_map<value_t,count_t>(const tallyman_map<value_t,count_t>&) = delete;
        tallyman_map<value_t,count_t>& operator=(const tallyman_map<value_t,count_t>&) = delete;

        virtual void tally(const std::vector<value_t>& ii) { for (auto i : ii) tally(i); }

        virtual bool is_map() const { return true; }

        virtual const std::unique_ptr<count_t[]>& get_results_vec() const;
        virtual const std::map<value_t,count_t>& get_results_map() const { return map_; }
};


// static create -------------------------------------------------------------

template<typename value_t, typename count_t>
tallyman<value_t,count_t>*
tallyman<value_t,count_t>::create(int nbits, int max_gb)
{
    const int max_bits = 8*sizeof(value_t);

    if (nbits < 1)
        raise_error("number of bits must be one or more: %d", nbits);

    if (nbits > max_bits)
        raise_error("%d bits requested exceeds limit %d", nbits, max_bits);

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

        return new tallyman_map<value_t,count_t>(nbits);
    }
    else
    {
        verbose_emit("vector memory (%luG) fits %luG: tallying using linear vector",
                static_cast<unsigned long>(vec_mb >> 10),
                static_cast<unsigned long>(max_mb >> 10));

        return new tallyman_vec<value_t,count_t>(nbits);
    }
}

// constructors --------------------------------------------------------------

template<typename value_t, typename count_t>
tallyman<value_t,count_t>::tallyman(int nbits) 
    : max_value_((1<<nbits)-1), n_invalid_(0)
{
    const int max_bits = 8*sizeof(value_t);
    if (nbits > 8*(int)sizeof(value_t))
        raise_error("number of bits (%d) exceeds maximum %d", nbits, max_bits);
}

template<typename value_t, typename count_t>
tallyman_vec<value_t,count_t>::tallyman_vec(int nbits)
    : tallyman<value_t,count_t>(nbits), vec_(nullptr)
{
    size_t alloc_n = tallyman<value_t,count_t>::max_value_ + 1;
    size_t alloc_size = alloc_n * sizeof(count_t);

    count_t *data = new count_t[alloc_n];
    if (!data)
        raise_error("failed to allocate memory (%luMB) for tally vector",
                static_cast<unsigned long>(alloc_size >> 20));

    std::memset(data, 0, alloc_size);

    vec_.reset(data);
}

template<typename value_t, typename count_t>
tallyman_map<value_t,count_t>::tallyman_map(int nbits) 
    : tallyman<value_t,count_t>(nbits)
{
}

template<typename value_t, typename count_t>
inline void
tallyman_map<value_t,count_t>::tally(value_t i)
{
    if (i > tallyman<value_t,count_t>::max_value_)
        ++tallyman<value_t,count_t>::n_invalid_;
    else {
        tallyman_map<value_t,count_t>::iterator p = map_.lower_bound(i);
        if (p == map_.end() || i != p->first)
            map_.insert(p, std::make_pair(i,1));
        else
            ++(p->second);
    }
}

template<typename value_t, typename count_t>
inline void
tallyman_vec<value_t,count_t>::tally(value_t i)
{
    if (i > tallyman<value_t,count_t>::max_value_) 
        ++tallyman<value_t,count_t>::n_invalid_; 
    else 
        ++vec_[i]; 
}

template<typename value_t, typename count_t>
const std::unique_ptr<count_t[]>&
tallyman_map<value_t,count_t>::get_results_vec() const
{
    static std::unique_ptr<count_t[]> dummy;
    raise_error("invalid invocation: get_results_vec on map implementation");
    return dummy;
}

template<typename value_t, typename count_t>
const std::map<value_t,count_t>&
tallyman_vec<value_t,count_t>::get_results_map() const
{ 
    static std::map<value_t,count_t> dummy;
    raise_error("invalid invocation: get_results_map on vec implementation");
    return dummy;
}


} // namespace kfc

#endif // tallyman_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et