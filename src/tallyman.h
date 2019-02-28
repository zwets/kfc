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
#include <map>

namespace kfc {


// count_t - the datatype holding tallies
//
// If there could ever be more than 4G occurrences of a single kmer, pick
// typedef std::uint_fast64_t count_t;  // up to 16E counts
// Else, if you value speed at the cost of possible memory consumption, pick
// typedef std::uint_fast32_t count_t;  // up to 4G counts, maybe more, speed
// Else, if you value memory consumption over possible faster speed, pick
// typedef std::uint32_t count_t;       // up to 4G counts, save on memory
//
typedef std::uint_fast32_t count_t;


// tallyman - keeps counts of encoded kmers (or any uint32 or uint64 integers).
//
// Core operation is tally(n), which increments the count of integer n by 1,
// or increments the invalid_count if n is larger than nbits can contain.
//
// This class has multiple implementations with differing space and time
// characteristics.  Given B the number of bits in the integers to be tallied,
// N the number of integers tallied, then
// - tallyman_vec uses a linear array with O(1) lookup and O(2^B) memory;
//   the count_t defined above is a (possible) factor 2 in memory consumption
// - tallyman_map uses a map of O(N) storage and O(log N) lookup;
//   the count_t defined above is a possible factor N/2 in memory consumption
//
// The tallyman::create() factory method picks the fastest implementation that
// fits a given approximate memory limit.
//
// The get_results_X() members return the tallied counts.  For performance
// reasons, these members do not shield from the underlying implementation.
// Use the is_vec() and is_map32() selectors to find out which get_results_X()
// is applicable to the actual class.
//
class tallyman {

    public:
        typedef std::uint64_t val64_t;
        typedef std::uint32_t val32_t;

    protected:
        val64_t max_value_;
        count_t n_invalid_;

    public:
        static std::unique_ptr<tallyman> create(int nbits, int max_gb = 0);

    protected:
        tallyman(int nbits) : max_value_((1 << nbits) - 1), n_invalid_(0) { }

    public:
	virtual void tally(val32_t);
	virtual void tally(val64_t);
	virtual void tally(std::vector<val32_t>);
	virtual void tally(std::vector<val64_t>);

        virtual bool is_vec() const;
        virtual bool is_map32() const;

        virtual const std::unique_ptr<count_t[]>& get_results_vec() const;
        virtual const std::map<val32_t,count_t>& get_results_map32() const;
        virtual const std::map<val64_t,count_t>& get_results_map64() const;

        val64_t max_value() const { return max_value_; }
        count_t invalid_count() const { return n_invalid_; }
};

class tallyman_vec : public tallyman
{
    private:
        std::unique_ptr<count_t[]> vec_;

    public:
        tallyman_vec(int nbits);
        tallyman_vec(const tallyman_vec&) = delete;
        tallyman_vec& operator=(const tallyman_vec&) = delete;

	virtual void tally(val32_t);
	virtual void tally(val64_t);
	virtual void tally(std::vector<val32_t>);
	virtual void tally(std::vector<val64_t>);

        virtual bool is_vec() const;

        virtual const std::unique_ptr<count_t[]>& get_results_vec() const;
};

template<typename value_t>
class tallyman_map : public tallyman
{
    private:
        typedef typename std::map<value_t,count_t>::iterator iterator;
        std::map<value_t,count_t> map_;

    public:
        tallyman_map<value_t>(int nbits);
        tallyman_map<value_t>(const tallyman_map<value_t>&) = delete;
        tallyman_map<value_t>& operator=(const tallyman_map<value_t>&) = delete;

        virtual void tally(value_t);
        virtual void tally(std::vector<value_t>);

        virtual bool is_map32() const;

        virtual const std::map<val32_t,count_t>& get_results_map32() const;
        virtual const std::map<val64_t,count_t>& get_results_map64() const;
};


} // namespace kfc

#endif // tallyman_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
