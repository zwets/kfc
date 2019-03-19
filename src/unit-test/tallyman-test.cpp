/* tallyman-test.cpp
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

#include <list>
#include <gtest/gtest.h>
#include "tallyman.h"
#include "utils.h"

using namespace kfc;

namespace {

typedef tallyman<std::uint32_t,std::uint32_t> tman3232;
typedef tallyman<std::uint64_t,std::uint32_t> tman6432;
typedef tallyman<std::uint32_t,std::uint64_t> tman3264;
typedef tallyman<std::uint64_t,std::uint64_t> tman6464;
typedef tallyman<std::uint32_t,double> tman32dd;

typedef std::unique_ptr<tman3232> uptr3232;
typedef std::unique_ptr<tman6432> uptr6432;
typedef std::unique_ptr<tman3264> uptr3264;
typedef std::unique_ptr<tman6464> uptr6464;
typedef std::unique_ptr<tman32dd> uptr32dd;

typedef tallyman_vec<std::uint32_t,std::uint32_t> tvec3232;
typedef tallyman_vec<std::uint64_t,std::uint32_t> tvec6432;
typedef tallyman_vec<std::uint32_t,std::uint64_t> tvec3264;
typedef tallyman_vec<std::uint64_t,std::uint64_t> tvec6464;

typedef tallyman_map<std::uint32_t,std::uint32_t> tmap3232;
typedef tallyman_map<std::uint32_t,std::uint64_t> tmap3264;
typedef tallyman_map<std::uint64_t,std::uint32_t> tmap6432;
typedef tallyman_map<std::uint64_t,std::uint64_t> tmap6464;

typedef std::map<std::uint32_t,std::uint32_t> map3232;
typedef std::map<std::uint64_t,std::uint32_t> map6432;

static uptr3232 create_uptr3232(int nbits, bool map = false) {
    return uptr3232(map ? (tman3232*) new tmap3232(nbits) : (tman3232*) new tvec3232(nbits));
}

static uptr3264 create_uptr3264(int nbits, bool map = false) {
    return uptr3264(map ? (tman3264*) new tmap3264(nbits) : (tman3264*) new tvec3264(nbits));
}

static uptr6432 create_uptr6432(int nbits, bool map = false) {
    return uptr6432(map ? (tman6432*) new tmap6432(nbits) : (tman6432*) new tvec6432(nbits));
}

static uptr6464 create_uptr6464(int nbits, bool map = false) {
    return uptr6464(map ? (tman6464*) new tmap6464(nbits) : (tman6464*) new tvec6464(nbits));
}

static uptr32dd create_uptr32dd(int nbits) {
    return uptr32dd(new tallyman_vec<std::uint32_t, double>(nbits));
}

TEST(tallyman_test, no_bits_zero) {
    EXPECT_DEATH(new tvec3232(0), ".*");
}

TEST(tallyman_test, double_counter) {
    auto p = create_uptr32dd(4);
    p->tally({0,1,1,2,2,2});
    EXPECT_EQ(p->get_results_vec()[0],1.0);
}

TEST(tallyman_test, no_bits_33) {
    EXPECT_DEATH(new tmap3232(33), ".*");
    EXPECT_DEATH(new tvec3232(33), ".*");
    EXPECT_DEATH(new tmap3264(33), ".*");
    EXPECT_DEATH(new tvec3264(33), ".*");
}

TEST(tallyman_test, no_bits_65) {
    EXPECT_DEATH(new tmap6432(65), ".*");
    EXPECT_DEATH(new tvec6432(65), ".*");
    EXPECT_DEATH(new tmap6464(65), ".*");
    EXPECT_DEATH(new tvec6464(65), ".*");
}

TEST(tallyman_test, bits_1_vec) {
    uptr3232 r1 = create_uptr3232(1);
    uptr6432 r2 = create_uptr6432(1);
    uptr3264 r3 = create_uptr3264(1);
    uptr6464 r4 = create_uptr6464(1);
}

TEST(tallyman_test, bits_1_map) {
    uptr3232 r1 = create_uptr3232(1, true);
    uptr6432 r2 = create_uptr6432(1, true);
    uptr3264 r3 = create_uptr3264(1, true);
    uptr6464 r4 = create_uptr6464(1, true);
}

TEST(tallyman_test, bits_64) {
    uptr6432 r1 = create_uptr6432(64);
    uptr6464 r2 = create_uptr6464(64);
}

TEST(tallyman_test, bits_64_map) {
    uptr6432 r1 = create_uptr6432(64,true);
    uptr6464 r2 = create_uptr6464(64,true);
}

TEST(tallyman_test, store_none) {
    uptr3232 r = create_uptr3232(2);
    EXPECT_TRUE(r->is_vec());
    const std::uint32_t* v = r->get_results_vec();
    EXPECT_EQ(v[0],0);
}

TEST(tallyman_test, store_one) {
    uptr3264 r = create_uptr3264(2);
    r->tally({1});
    EXPECT_TRUE(r->is_vec());
    const uint64_t* v = r->get_results_vec();
    EXPECT_EQ(v[1], 1);
}

TEST(tallyman_test, store_two_ones) {
    uptr6432 r = create_uptr6432(2);
    r->tally({1,1});
    EXPECT_TRUE(r->is_vec());
    const uint32_t* v = r->get_results_vec();
    EXPECT_EQ(v[1], 2);
}

TEST(tallyman_test, store_invalid) {
    uptr3232 r = create_uptr3232(2);
    r->tally({4,3});
    EXPECT_EQ(r->invalid_count(),1);
    EXPECT_TRUE(r->is_vec());
    const uint32_t* v = r->get_results_vec();
    EXPECT_EQ(v[3], 1);
}

TEST(tallyman_test, store_invalid_64) {
    uptr6432 r = create_uptr6432(2);
    r->tally({4,3});
    EXPECT_EQ(r->invalid_count(),1);
    EXPECT_TRUE(r->is_vec());
    const uint32_t* v = r->get_results_vec();
    EXPECT_EQ(v[3], 1);
}

TEST(tallyman_test, store_map_none) {
    uptr3232 r = create_uptr3232(29,true);
    EXPECT_TRUE(r->is_map());
    const map3232& m = r->get_results_map();
    EXPECT_EQ(m.begin(),m.end());
}

TEST(tallyman_test, store_map_one) {
    uptr3232 r = create_uptr3232(29,true);
    r->tally({1234567});
    EXPECT_TRUE(r->is_map());
    const map3232& m = r->get_results_map();
    map3232::const_iterator i = m.cbegin();
    EXPECT_EQ(i->first, 1234567);
    EXPECT_EQ(i->second, 1);
    EXPECT_EQ(++i, m.cend());
}

TEST(tallyman_test, store_map_two_ones) {
    uptr3232 r = create_uptr3232(29,true);
    r->tally({7654321,7654321});
    EXPECT_TRUE(r->is_map());
    const map3232& m = r->get_results_map();
    map3232::const_iterator i = m.cbegin();
    EXPECT_EQ(i->first, 7654321);
    EXPECT_EQ(i->second, 2);
    EXPECT_EQ(++i, m.cend());
}

TEST(tallyman_test, store_map_invalid) {
    uptr3232 r = create_uptr3232(29,true);
    r->tally({std::uint32_t(1)<<29,(std::uint32_t(1)<<29)-1});
    EXPECT_EQ(r->invalid_count(),1);
    EXPECT_TRUE(r->is_map());
    const map3232& m = r->get_results_map();
    map3232::const_iterator i = m.cbegin();
    EXPECT_EQ(i->first, (std::uint32_t(1)<<29)-1);
    EXPECT_EQ(i->second, 1);
    EXPECT_EQ(++i, m.cend());
}

TEST(tallyman_test, store_map_invalid_64) {
    uptr6432 r = create_uptr6432(29,true);
    r->tally({std::uint64_t(1)<<29,(std::uint64_t(1)<<29)-1});
    EXPECT_EQ(r->invalid_count(),1);
    EXPECT_TRUE(r->is_map());
    const map6432& m = r->get_results_map();
    map6432::const_iterator i = m.cbegin();
    EXPECT_EQ(i->first, (std::uint64_t(1)<<29)-1);
    EXPECT_EQ(i->second, 1);
    EXPECT_EQ(++i, m.cend());
}

} // namespace
// vim: sts=4:sw=4:ai:si:et
