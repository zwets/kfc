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


// tallyman - keeps counts of encoded kmers, or generally of any
//            unsigned integral type of a specified maximum number of bits
//
// Tallyman has implementations with different space and time characteristics.
// Given B=nbits the bit size of the items to be tallied, N the number of
// values tallied, and C the size of count_t, then:
// - tallyman_vec uses a linear array, with O(1) lookup and C*2^B memory;
// - tallyman_map uses a map, with O(log N) lookup and O(N) storage
//
// The tallyman::create(nbits,max_gb) factory method returns the vector
// implementation if it fits in max_gb memory, or else the map implementation.
//
// The core operation is tally(items), which tallies each i in items by either
// incrementing its item count, or incrementing the invalid_count if i exceeds
// max_value, the largest possible nbit number.
//
// Template parameter value_t must be an unsigned integral type of at least
// nbits bits (or an exception is thrown).  If its bit size equals nbits, then
// by definition invalid_count will remain 0, as no item can exceed max_value.
//
// Template parameter count_t can be any numeric type, though integral types
// will likely perform better than floating point.  It must be wide enough
// to hold the maximum tally of any element (and of the invalid count).
// Roll-over of integral types is silent; this is C++.
//
// The get_results_X() members return the tallied counts.  For performance
// reasons, these members do not shield from the underlying implementation
// (vector or map).  Use the is_vec() and is_map() selectors to find out
// which get_results_X() is applicable to the actual class.
//
template <typename value_t, typename count_t>
class tallyman {
    static_assert(std::is_unsigned<value_t>::value, 
            "template argument value_t must be an unsigned integral type");
    static_assert(std::is_integral<count_t>::value || std::is_floating_point<count_t>::value, 
            "template argument count_t must be a numerical type");

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
    static_assert(std::is_unsigned<value_t>::value,
            "template argument value_t must be unsigned integral");
    static_assert(std::is_integral<count_t>::value || std::is_floating_point<count_t>::value, 
            "template argument count_t must be a numerical type");

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
    static_assert(std::is_unsigned<value_t>::value,
            "template argument value_t must be unsigned integral");
    static_assert(std::is_integral<count_t>::value || std::is_floating_point<count_t>::value, 
            "template argument count_t must be a numerical type");

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
    constexpr int max_bits = 8*sizeof(value_t);

    if (nbits < 1)
        raise_error("number of bits must be one or more: %d", nbits);

    if (nbits > max_bits)
        raise_error("%d bits requested exceeds limit %d", nbits, max_bits);

    unsigned long max_mb = static_cast<unsigned long>(max_gb) << 10;
    unsigned long vec_mb = sizeof(count_t) * (1UL << std::max(0, nbits-20));

    if (max_gb == 0)
    {
        unsigned long phy_mb = get_system_memory() >> 20;
        max_mb = phy_mb > 2048 ? phy_mb - 2048 : phy_mb;

        verbose_emit("defaulting max memory to all%s physical memory: %luG",
                phy_mb > 2048 ? " but 2G" : "", (max_mb>>10));
    }

    if (vec_mb > max_mb)
    {
        verbose_emit("vector memory (%luG) would exceed %luG: tallying using map",
                (vec_mb>> 10), (max_mb>>10));
        return new tallyman_map<value_t,count_t>(nbits);
    }
    else
    {
        verbose_emit("vector memory (%luG) fits %luG: tallying using linear vector",
                (vec_mb>>10), (max_mb>>10));
        return new tallyman_vec<value_t,count_t>(nbits);
    }
}

// constructors --------------------------------------------------------------

template<typename value_t, typename count_t>
tallyman<value_t,count_t>::tallyman(int nbits) 
    : max_value_((1<<nbits)-1), n_invalid_(0)
{
    constexpr int max_bits = 8*sizeof(value_t);
    if (nbits > max_bits)
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
