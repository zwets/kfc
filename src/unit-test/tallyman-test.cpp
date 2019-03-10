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

typedef tallyman_map<std::uint32_t,std::uint32_t> tmap3232;
typedef tallyman_map<std::uint32_t,std::uint64_t> tmap3264;
typedef tallyman_map<std::uint64_t,std::uint32_t> tmap6432;

typedef std::map<std::uint32_t,std::uint32_t> map3232;

static uptr3232 create_uptr3232(int nbits, int max_gb = 0) {
    return uptr3232(tman3232::create(nbits, max_gb));
}

static uptr3264 create_uptr3264(int nbits, int max_gb = 0) {
    return uptr3264(tman3264::create(nbits, max_gb));
}

static uptr6432 create_uptr6432(int nbits, int max_gb = 0) {
    return uptr6432(tman6432::create(nbits, max_gb));
}

static uptr32dd create_uptr32dd(int nbits, int max_gb = 0) {
    return uptr32dd(tman32dd::create(nbits, max_gb));
}

TEST(tallyman_test, no_bits_zero) {
    EXPECT_DEATH(create_uptr3232(0), ".*");
}

TEST(tallyman_test, double_counter) {
    auto p = create_uptr32dd(4);
    p->tally({0,1,1,2,2,2});
    EXPECT_EQ(p->get_results_vec()[0],1.0);
}

TEST(tallyman_test, no_bits_33) {
    EXPECT_DEATH(tman3232::create(33), ".*");
}

TEST(tallyman_test, no_bits_65) {
    EXPECT_DEATH(tman6432::create(65), ".*");
}

TEST(tallyman_test, no_bits_65_either) {
    EXPECT_DEATH(tman6464::create(65), ".*");
}

TEST(tallyman_test, bits_1) {
    uptr3232 r = create_uptr3232(1);
}

TEST(tallyman_test, bits_64) {
    uptr6432 r = create_uptr6432(64);
}

TEST(tallyman_test, bits_1_is_vec) {
    uptr3232 r = create_uptr3232(1);
    tman3232 *p = dynamic_cast<tvec3232*>(r.get());
    EXPECT_NE(p, (tman3232*)0);
}

TEST(tallyman_test, bits_64_is_map) {
    uptr6432 r = create_uptr6432(64);
    tman6432 *p = dynamic_cast<tmap6432*>(r.get());
    EXPECT_NE(p, (tman6432*)0);
}

TEST(tallyman_test, bits_32_map) {
    uptr3232 r = create_uptr3232(32,1);
    tman3232 *p = dynamic_cast<tmap3232*>(r.get());
    EXPECT_NE(p, (tman3232*)0);
}

TEST(tallyman_test, bits_33_map) {
    uptr6432 r = create_uptr6432(33,1);
    tman6432 *p = dynamic_cast<tmap6432*>(r.get());
    EXPECT_NE(p, (tman6432*)0);
}

TEST(tallyman_test, bits_27_1g_vec) {
    uptr3232 r = create_uptr3232(27,1);
    tman3232 *p = dynamic_cast<tvec3232*>(r.get());
    EXPECT_NE(p, (tman3232*)0);
}

TEST(tallyman_test, bits_28_1g_depends) {
    uptr3232 r0 = create_uptr3232(28,1);
    tman3232 *p0 = dynamic_cast<tvec3232*>(r0.get());
    EXPECT_NE(p0, (tman3232*)0);
    uptr3264 r1 = create_uptr3264(28,1);
    tman3264 *p1 = dynamic_cast<tmap3264*>(r1.get());
    EXPECT_NE(p1, (tman3264*)0);
}

TEST(tallyman_test, bits_29_1g_map) {
    uptr3232 r = create_uptr3232(29,1);
    tman3232 *p = dynamic_cast<tmap3232*>(r.get());
    EXPECT_NE(p, (tman3232*)0);
}

TEST(tallyman_test, store_none) {
    uptr3232 r = create_uptr3232(2);
    EXPECT_TRUE(r->is_vec());
    const std::uint32_t* v = r->get_results_vec();
    EXPECT_EQ(v[0],0);
}

TEST(tallyman_test, store_one) {
    uptr3232 r = create_uptr3232(2);
    r->tally({1});
    EXPECT_TRUE(r->is_vec());
    const uint32_t* v = r->get_results_vec();
    EXPECT_EQ(v[1], 1);
}

TEST(tallyman_test, store_two_ones) {
    uptr3232 r = create_uptr3232(2);
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

TEST(tallyman_test, store_map_none) {
    uptr3232 r = create_uptr3232(29,1);
    EXPECT_TRUE(r->is_map());
    const map3232& m = r->get_results_map();
    EXPECT_EQ(m.begin(),m.end());
}

TEST(tallyman_test, store_map_one) {
    uptr3232 r = create_uptr3232(29,1);
    r->tally({1234567});
    EXPECT_TRUE(r->is_map());
    const map3232& m = r->get_results_map();
    map3232::const_iterator i = m.cbegin();
    EXPECT_EQ(i->first, 1234567);
    EXPECT_EQ(i->second, 1);
    EXPECT_EQ(++i, m.cend());
}

TEST(tallyman_test, store_map_two_ones) {
    uptr3232 r = create_uptr3232(29,1);
    r->tally({7654321,7654321});
    EXPECT_TRUE(r->is_map());
    const map3232& m = r->get_results_map();
    map3232::const_iterator i = m.cbegin();
    EXPECT_EQ(i->first, 7654321);
    EXPECT_EQ(i->second, 2);
    EXPECT_EQ(++i, m.cend());
}

TEST(tallyman_test, store_map_invalid) {
    uptr3232 r = create_uptr3232(29,1);
    r->tally({std::uint32_t(1)<<29,(std::uint32_t(1)<<29)-1});
    EXPECT_EQ(r->invalid_count(),1);
    EXPECT_TRUE(r->is_map());
    const map3232& m = r->get_results_map();
    map3232::const_iterator i = m.cbegin();
    EXPECT_EQ(i->first, (std::uint32_t(1)<<29)-1);
    EXPECT_EQ(i->second, 1);
    EXPECT_EQ(++i, m.cend());
}

} // namespace
// vim: sts=4:sw=4:ai:si:et
