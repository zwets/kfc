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
#ifndef NO_THREADS
#  include <mutex>
#endif
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
// The core operation is tally(items), which tallies each i in items by either
// incrementing its item count, or incrementing the invalid_count if i exceeds
// max_value, the largest possible nbit number.
//
// Template parameter value_t must be an unsigned integral type of at least
// nbits bits (or the program exits).  If its bit size equals nbits, then
// by definition invalid_count will remain 0, as no item can exceed max_value.
//
// Template parameter count_t can be any numeric type, though integral types
// will likely perform better than floating point.  It must be wide enough
// to hold the maximum tally of any element (and of the invalid count).
// Roll-over of integral types is silent.
//
// The get_results_X() members return the tallied counts.  For performance
// reasons, these members do not shield from the underlying implementation
// (vector or map).  Use the is_vec() and is_map() selectors to find out
// whether get_results_vec() or get_results_map() should be called.
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

    protected:
        tallyman(int nbits);

    public:
        virtual ~tallyman() { }

	virtual void tally(std::vector<value_t> &&) = 0;
	virtual void tally(const std::vector<value_t>&) = 0;

        virtual bool is_vec() const { return false; }
        virtual bool is_map() const { return false; }

        virtual const count_t *get_results_vec() const = 0;
        virtual const std::map<value_t,count_t>& get_results_map() const = 0;

        value_t max_value() const { return max_value_; }
        count_t invalid_count() const { return n_invalid_; }
};

template <typename value_t, typename count_t>
class tallyman_vec : public tallyman<value_t,count_t>
{
    static_assert(std::is_unsigned<value_t>::value,
            "template argument value_t must be unsigned integral");
    static_assert(std::is_integral<count_t>::value || std::is_floating_point<count_t>::value,
            "template argument count_t must be a numerical type");

    private:
        count_t *vec_;
#ifndef NO_THREADS
        std::mutex tally_mutex_;
#endif
        void tally(value_t i);

    public:
        tallyman_vec<value_t,count_t>(int nbits);
        tallyman_vec<value_t,count_t>(const tallyman_vec&) = delete;
        tallyman_vec<value_t,count_t>& operator=(const tallyman_vec&) = delete;
        virtual ~tallyman_vec<value_t,count_t>() { if (vec_) free(vec_); }

	virtual void tally(std::vector<value_t> &&);
	virtual void tally(const std::vector<value_t>&);

        virtual bool is_vec() const { return true; }

        virtual const count_t *get_results_vec() const { return vec_; }
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
#ifndef NO_THREADS
        std::mutex tally_mutex_;
#endif
        void tally(value_t i);

    public:
        tallyman_map<value_t,count_t>(int nbits);
        tallyman_map<value_t,count_t>(const tallyman_map<value_t,count_t>&) = delete;
        tallyman_map<value_t,count_t>& operator=(const tallyman_map<value_t,count_t>&) = delete;

        virtual void tally(std::vector<value_t>&&);
        virtual void tally(const std::vector<value_t>&);

        virtual bool is_map() const { return true; }

        virtual const count_t *get_results_vec() const;
        virtual const std::map<value_t,count_t>& get_results_map() const { return map_; }
};

// constructors --------------------------------------------------------------

template<typename value_t, typename count_t>
tallyman<value_t,count_t>::tallyman(int nbits)
    : max_value_((static_cast<value_t>(1)<<nbits)-1), n_invalid_(0)
{
    constexpr int max_bits = 8*sizeof(value_t);
    if (nbits < 1)
        raise_error("invalid number of bits: %d", nbits);
    if (nbits > max_bits)
        raise_error("number of bits (%d) exceeds maximum %d", nbits, max_bits);
}

template<typename value_t, typename count_t>
tallyman_vec<value_t,count_t>::tallyman_vec(int nbits)
    : tallyman<value_t,count_t>(nbits), vec_(0)
{
    size_t alloc_n = tallyman<value_t,count_t>::max_value_ + 1;
    size_t alloc_size = alloc_n * sizeof(count_t);

    vec_ = (count_t*) std::calloc(alloc_n, sizeof(count_t));
    if (!vec_)
        raise_error("failed to allocate memory (%luMB) for tally vector",
                static_cast<unsigned long>(alloc_size >> 20));
}

template<typename value_t, typename count_t>
tallyman_map<value_t,count_t>::tallyman_map(int nbits)
    : tallyman<value_t,count_t>(nbits)
{
}

// tallyman_vec --------------------------------------------------------------

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
inline void
tallyman_vec<value_t,count_t>::tally(std::vector<value_t> &&ii)
{
#ifndef NO_THREADS
    std::unique_lock<std::mutex> tally_lock(tally_mutex_);
#endif
    for (auto i : ii)
        tally(i);
}

template<typename value_t, typename count_t>
inline void
tallyman_vec<value_t,count_t>::tally(const std::vector<value_t>& ii)
{
#ifndef NO_THREADS
    std::unique_lock<std::mutex> tally_lock(tally_mutex_);
#endif
    for (auto i : ii)
        tally(i);
}

template<typename value_t, typename count_t>
const std::map<value_t,count_t>&
tallyman_vec<value_t,count_t>::get_results_map() const
{
    static std::map<value_t,count_t> dummy;
    raise_error("invalid invocation: get_results_map on vec implementation");
    return dummy;
}

// tallyman_map --------------------------------------------------------------

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
tallyman_map<value_t,count_t>::tally(const std::vector<value_t>& ii)
{
#ifndef NO_THREADS
    std::unique_lock<std::mutex> tally_lock(tally_mutex_);
#endif
    for (auto i : ii)
        tally(i);
}

template<typename value_t, typename count_t>
inline void
tallyman_map<value_t,count_t>::tally(std::vector<value_t> &&ii)
{
#ifndef NO_THREADS
    std::unique_lock<std::mutex> tally_lock(tally_mutex_);
#endif
    for (auto i : ii)
        tally(i);
}

template<typename value_t, typename count_t>
const count_t*
tallyman_map<value_t,count_t>::get_results_vec() const
{
    raise_error("invalid invocation: get_results_vec on map implementation");
    return 0;
}


} // namespace kfc

#endif // tallyman_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
