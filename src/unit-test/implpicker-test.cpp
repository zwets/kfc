/* implpicker-test.cpp
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

#include <gtest/gtest.h>
#include "implpicker.h"
#include "utils.h"

using namespace kfc;

namespace {

typedef std::uint32_t u32;
typedef std::uint64_t u64;

std::unique_ptr<kmer_counter>
pick_impl_wrap(int ks, bool ss = false, unsigned mc = 0, unsigned mg = 0, char fi = '\0', unsigned nt = 0)
{
    return std::unique_ptr<kmer_counter>(pick_implementation(ks, ss, mc, mg, fi, nt));
}

bool is_tally3232(kmer_counter *p) {
    return dynamic_cast<kmer_counter_tally<u32,u32>*>(p) != nullptr;
}

bool is_tally6432(kmer_counter *p) {
    return dynamic_cast<kmer_counter_tally<u64,u32>*>(p) != nullptr;
}

bool is_tally3264(kmer_counter *p) {
    return dynamic_cast<kmer_counter_tally<u32,u64>*>(p) != nullptr;
}

//bool is_tally6464(kmer_counter *p) {
//    return dynamic_cast<kmer_counter_tally<u64,u64>*>(p) != nullptr;
//}

bool is_list32(kmer_counter *p) {
    return dynamic_cast<kmer_counter_list<u32>*>(p) != nullptr;
}

bool is_list64(kmer_counter *p) {
    return dynamic_cast<kmer_counter_list<u64>*>(p) != nullptr;
}

// sizes and limits -----------------------------------------------------

TEST(implpicker_test, no_ksize_zero) {
    EXPECT_DEATH(pick_impl_wrap(0), ".*");
}

TEST(implpicker_test, ksize_1) {
    pick_impl_wrap(1);
    pick_impl_wrap(1, true);
}

TEST(implpicker_test, no_ksize_16) {
    EXPECT_DEATH(pick_impl_wrap(16), ".*");
}

TEST(implpicker_test, ksize_15) {
    pick_impl_wrap(15,false,1);
    pick_impl_wrap(15,true,1);
}

TEST(implpicker_test, no_ksize_32) {
    EXPECT_DEATH(pick_impl_wrap(32), ".*");
    EXPECT_DEATH(pick_impl_wrap(32), ".*");
}

TEST(implpicker_test, no_ksize_32_ss) {
    EXPECT_DEATH(pick_impl_wrap(32,true), ".*");
    EXPECT_DEATH(pick_impl_wrap(32,true), ".*");
}

TEST(implpicker_test, ksize_31) {
    pick_impl_wrap(31,false,1);
    pick_impl_wrap(31,true,1);
}

TEST(implpicker_test, no_ksize_even) {
    EXPECT_DEATH(pick_impl_wrap(6), ".*");
}

TEST(implpicker_test, ksize_even_ss) {
    pick_impl_wrap(6,true,1);
}

// types returned --------------------------------------------------------

TEST(implpicker_test, small_k_is_tally) {
    EXPECT_TRUE(is_tally3232(pick_impl_wrap(7,false,1).get()));
}

TEST(implpicker_test, large_count_is_tally) {
    EXPECT_TRUE(is_tally3264(pick_impl_wrap(7,false,1<<13).get()));
}

TEST(implpicker_test, small_count_is_list) {
    EXPECT_TRUE(is_list32(pick_impl_wrap(15,true,2).get()));
}

TEST(implpicker_test, big_ksize_small_count_is_list64) {
    EXPECT_TRUE(is_list64(pick_impl_wrap(24,true,2).get()));
}

TEST(implpicker_test, big_ksize_big_count_is_list64) {
    EXPECT_TRUE(is_list64(pick_impl_wrap(17,false,1<<13).get()));
}

// errors for user specified ----------------------------------------------

TEST(implpicker_test, exceed_1g) {
    EXPECT_DEATH(pick_impl_wrap(15,true,1<<10,1), ".*");
}

TEST(implpicker_test, not_valid_impl) {
    EXPECT_DEATH(pick_impl_wrap(5,true,1,1,'x'), ".*");
}

TEST(implpicker_test, force_vec_impl64_too_big) {
    EXPECT_DEATH(is_tally6432(pick_impl_wrap(30,true,2,1,'v').get()), ".*");
}

// picking implementations ------------------------------------------------

TEST(implpicker_test, force_map_impl) {
    EXPECT_TRUE(is_tally3232(pick_impl_wrap(15,true,2,0,'m').get()));
}

TEST(implpicker_test, force_vec_impl) {
    EXPECT_TRUE(is_tally3232(pick_impl_wrap(15,true,2,0,'v').get()));
}

TEST(implpicker_test, force_lst_impl) {
    EXPECT_TRUE(is_list32(pick_impl_wrap(15,true,2,0,'l').get()));
}

TEST(implpicker_test, force_lst_impl64) {
    EXPECT_TRUE(is_list64(pick_impl_wrap(17,false,2,0,'l').get()));
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
