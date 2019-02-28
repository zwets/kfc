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

#include <stdexcept>
#include <list>
#include <gtest/gtest.h>
#include "tallyman.h"
#include "utils.h"

using namespace kfc;

namespace {

TEST(tallyman_test, no_bits_zero) {
    std::unique_ptr<tallyman> r;
    //EXPECT_THROW(r = new tallyman(0), std::runtime_error);
    EXPECT_DEATH(r = tallyman::create(0), ".*");
}

TEST(tallyman_test, no_bits_65) {
    std::unique_ptr<tallyman> r;
    //EXPECT_THROW(r = new tallyman(tallyman::max_ksize+1), std::runtime_error);
    EXPECT_DEATH(r = tallyman::create(65), ".*");
}

TEST(tallyman_test, bits_1) {
    std::unique_ptr<tallyman> r = tallyman::create(1);
}

TEST(tallyman_test, bits_64) {
    std::unique_ptr<tallyman> r = tallyman::create(64);
}

TEST(tallyman_test, bits_1_is_vec) {
    std::unique_ptr<tallyman> r = tallyman::create(1);
    tallyman *p = dynamic_cast<tallyman_vec*>(r.get());
    EXPECT_NE(p, (tallyman*)0);
}

TEST(tallyman_test, bits_64_is_map) {
    std::unique_ptr<tallyman> r = tallyman::create(64);
    tallyman *p = dynamic_cast<tallyman_map<tallyman::val64_t>*>(r.get());
    EXPECT_NE(p, (tallyman*)0);
}

TEST(tallyman_test, bits_32_map) {
    std::unique_ptr<tallyman> r = tallyman::create(32, 1);
    tallyman *p = dynamic_cast<tallyman_map<tallyman::val32_t>*>(r.get());
    EXPECT_NE(p, (tallyman*)0);
}

TEST(tallyman_test, bits_33_map) {
    std::unique_ptr<tallyman> r = tallyman::create(32, 1);
    tallyman *p = dynamic_cast<tallyman_map<tallyman::val32_t>*>(r.get());
    EXPECT_NE(p, (tallyman*)0);
}

TEST(tallyman_test, bits_27_1g_vec) {
    std::unique_ptr<tallyman> r = tallyman::create(27, 1);
    tallyman *p0 = r.get();
    EXPECT_NE(p0, (tallyman*)0);
    tallyman *p1 = dynamic_cast<tallyman_vec*>(p0);
    EXPECT_NE(p1, (tallyman*)0);
}

TEST(tallyman_test, bits_28_1g_which) {
    std::unique_ptr<tallyman> r = tallyman::create(28, 1);
    tallyman *p0 = r.get();
    EXPECT_NE(p0, (tallyman*)0);
    tallyman *p1 = sizeof(kfc::count_t) == 4 
        ? (tallyman*)dynamic_cast<tallyman_vec*>(p0)
        : (tallyman*)dynamic_cast<tallyman_map<tallyman::val32_t>*>(p0);
    EXPECT_NE(p1, (tallyman*)0);
}

TEST(tallyman_test, bits_29_1g_map) {
    std::unique_ptr<tallyman> r = tallyman::create(29, 1);
    tallyman *p0 = r.get();
    EXPECT_NE(p0, (tallyman*)0);
    tallyman *p1 = dynamic_cast<tallyman_map<tallyman::val32_t>*>(p0);
    EXPECT_NE(p1, (tallyman*)0);
}

TEST(tallyman_test, store_none) {
    std::unique_ptr<tallyman> r = tallyman::create(2);
    EXPECT_TRUE(r->is_vec());
    const std::unique_ptr<kfc::count_t[]> &v = r->get_results_vec();
    EXPECT_EQ(v[0],0);
}

TEST(tallyman_test, store_one) {
    std::unique_ptr<tallyman> r = tallyman::create(2);
    r->tally(tallyman::val32_t(1));

    EXPECT_TRUE(r->is_vec());
    const std::unique_ptr<kfc::count_t[]> &v = r->get_results_vec();
    EXPECT_EQ(v[1], 1);
}

TEST(tallyman_test, store_two_ones) {
    std::unique_ptr<tallyman> r = tallyman::create(2);
    r->tally(tallyman::val32_t(1));
    r->tally(tallyman::val32_t(1));

    EXPECT_TRUE(r->is_vec());
    const std::unique_ptr<kfc::count_t[]> &v = r->get_results_vec();
    EXPECT_EQ(v[1], 2);
}

TEST(tallyman_test, store_invalid) {
    std::unique_ptr<tallyman> r = tallyman::create(2);

    r->tally(tallyman::val32_t(4));
    r->tally(tallyman::val32_t(3));

    EXPECT_EQ(r->invalid_count(),1);
    EXPECT_TRUE(r->is_vec());
    const std::unique_ptr<kfc::count_t[]> &v = r->get_results_vec();
    EXPECT_EQ(v[3], 1);
}

TEST(tallyman_test, store_map_none) {
    std::unique_ptr<tallyman> r = tallyman::create(29,1);

    EXPECT_TRUE(r->is_map32());
    const std::map<tallyman::val32_t,kfc::count_t> &m = r->get_results_map32();
    EXPECT_EQ(m.begin(),m.end());
}

TEST(tallyman_test, store_map_one) {
    std::unique_ptr<tallyman> r = tallyman::create(29,1);
    r->tally(tallyman::val32_t(1234567));

    EXPECT_TRUE(r->is_map32());
    std::map<tallyman::val32_t,kfc::count_t>::const_iterator i = r->get_results_map32().begin();
    EXPECT_EQ(i->first, 1234567);
    EXPECT_EQ(i->second, 1);
    EXPECT_EQ(++i, r->get_results_map32().end());
}

TEST(tallyman_test, store_map_two_ones) {
    std::unique_ptr<tallyman> r = tallyman::create(29,1);

    r->tally(tallyman::val32_t(7654321));
    r->tally(tallyman::val32_t(7654321));

    EXPECT_TRUE(r->is_map32());
    std::map<tallyman::val32_t,kfc::count_t> m = r->get_results_map32();
    std::map<tallyman::val32_t,kfc::count_t>::const_iterator i = m.begin();
    EXPECT_EQ(i->first, 7654321);
    EXPECT_EQ(i->second, 2);
    EXPECT_EQ(++i, m.end());
}

TEST(tallyman_test, store_map_invalid) {
    std::unique_ptr<tallyman> r = tallyman::create(29,1);

    r->tally(tallyman::val32_t(1<<29));
    r->tally(tallyman::val32_t((1<<29)-1));

    EXPECT_EQ(r->invalid_count(),1);
    EXPECT_TRUE(r->is_map32());
    std::map<tallyman::val32_t,kfc::count_t> m = r->get_results_map32();
    std::map<tallyman::val32_t,kfc::count_t>::const_iterator i = m.begin();
    EXPECT_EQ(i->first, (1<<29)-1);
    EXPECT_EQ(i->second, 1);
    EXPECT_EQ(++i, m.end());
}


} // namespace
// vim: sts=4:sw=4:ai:si:et
