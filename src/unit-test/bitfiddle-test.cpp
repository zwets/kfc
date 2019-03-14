/* bitfiddle-test.cpp
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
#include "bitfiddle.h"

using namespace kfc;

namespace {

// decode_base --------------------------------------------------------

TEST(bitfiddle_test, test_bitsize) {
    EXPECT_EQ(bitsize<std::uint32_t>, 32);
    EXPECT_EQ(bitsize<std::uint64_t>, 64);
}

TEST(bitfiddle_test, test_low_bits32) {
    const uint32_t low_bits32_0 = low_bits<std::uint32_t,0>;
    const uint32_t low_bits32_1 = low_bits<std::uint32_t,1>;
    const uint32_t low_bits32_2 = low_bits<std::uint32_t,2>;
    const uint32_t low_bits32_3 = low_bits<std::uint32_t,3>;
    const uint32_t low_bits32_4 = low_bits<std::uint32_t,4>;
    const uint32_t low_bits32_29 = low_bits<std::uint32_t,29>;
    const uint32_t low_bits32_30 = low_bits<std::uint32_t,30>;
    const uint32_t low_bits32_31 = low_bits<std::uint32_t,31>;
    // doesn't compile: const uint32_t low_bits32_32 = low_bits<std::uint32_t,32>;

    EXPECT_EQ(low_bits32_0, 0x00000000);
    EXPECT_EQ(low_bits32_1, 0x00000001);
    EXPECT_EQ(low_bits32_2, 0x00000003);
    EXPECT_EQ(low_bits32_3, 0x00000007);
    EXPECT_EQ(low_bits32_4, 0x0000000F);
    EXPECT_EQ(low_bits32_29, 0x1FFFFFFF);
    EXPECT_EQ(low_bits32_30, 0x3FFFFFFF);
    EXPECT_EQ(low_bits32_31, 0x7FFFFFFF);
    //EXPECT_EQ(low_bits32_32, 0xFFFFFFFF);
}

TEST(bitfiddle_test, test_high_bits32) {

    // doesn't compile: const uint32_t high_bits32_0 = high_bits<std::uint32_t,0>;
    const uint32_t high_bits32_1 = high_bits<std::uint32_t,1>;
    const uint32_t high_bits32_2 = high_bits<std::uint32_t,2>;
    const uint32_t high_bits32_3 = high_bits<std::uint32_t,3>;
    const uint32_t high_bits32_4 = high_bits<std::uint32_t,4>;
    const uint32_t high_bits32_28 = high_bits<std::uint32_t,28>;
    const uint32_t high_bits32_29 = high_bits<std::uint32_t,29>;
    const uint32_t high_bits32_30 = high_bits<std::uint32_t,30>;
    const uint32_t high_bits32_31 = high_bits<std::uint32_t,31>;
    const uint32_t high_bits32_32 = high_bits<std::uint32_t,32>;

    // doesn't compile: EXPECT_EQ(high_bits32_0, 0x00000000);
    EXPECT_EQ(high_bits32_1, 0x80000000);
    EXPECT_EQ(high_bits32_2, 0xC0000000);
    EXPECT_EQ(high_bits32_3, 0xE0000000);
    EXPECT_EQ(high_bits32_4, 0xF0000000);
    EXPECT_EQ(high_bits32_28, 0xFFFFFFF0);
    EXPECT_EQ(high_bits32_29, 0xFFFFFFF8);
    EXPECT_EQ(high_bits32_30, 0xFFFFFFFC);
    EXPECT_EQ(high_bits32_31, 0xFFFFFFFE);
    EXPECT_EQ(high_bits32_32, 0xFFFFFFFF);
}

TEST(bitfiddle_test, test_low_bits64) {
    const uint64_t low_bits64_0 = low_bits<std::uint64_t,0>;
    const uint64_t low_bits64_1 = low_bits<std::uint64_t,1>;
    const uint64_t low_bits64_2 = low_bits<std::uint64_t,2>;
    const uint64_t low_bits64_3 = low_bits<std::uint64_t,3>;
    const uint64_t low_bits64_4 = low_bits<std::uint64_t,4>;
    const uint64_t low_bits64_61 = low_bits<std::uint64_t,61>;
    const uint64_t low_bits64_62 = low_bits<std::uint64_t,62>;
    const uint64_t low_bits64_63 = low_bits<std::uint64_t,63>;
    // doesn't compile: const uint64_t low_bits64_64 = low_bits<std::uint64_t,64>;

    EXPECT_EQ(low_bits64_0, 0x0000000000000000L);
    EXPECT_EQ(low_bits64_1, 0x0000000000000001L);
    EXPECT_EQ(low_bits64_2, 0x0000000000000003L);
    EXPECT_EQ(low_bits64_3, 0x0000000000000007L);
    EXPECT_EQ(low_bits64_4, 0x000000000000000FL);
    EXPECT_EQ(low_bits64_61, 0x1FFFFFFFFFFFFFFFL);
    EXPECT_EQ(low_bits64_62, 0x3FFFFFFFFFFFFFFFL);
    EXPECT_EQ(low_bits64_63, 0x7FFFFFFFFFFFFFFFL);
    //EXPECT_EQ(low_bits64_64, 0xFFFFFFFFFFFFFFFFL);
}

TEST(bitfiddle_test, test_high_bits64) {

    // doesn't compile: const uint64_t high_bits64_0 = high_bits<std::uint64_t,0>;
    const uint64_t high_bits64_1 = high_bits<std::uint64_t,1>;
    const uint64_t high_bits64_2 = high_bits<std::uint64_t,2>;
    const uint64_t high_bits64_3 = high_bits<std::uint64_t,3>;
    const uint64_t high_bits64_4 = high_bits<std::uint64_t,4>;
    const uint64_t high_bits64_61 = high_bits<std::uint64_t,61>;
    const uint64_t high_bits64_62 = high_bits<std::uint64_t,62>;
    const uint64_t high_bits64_63 = high_bits<std::uint64_t,63>;
    const uint64_t high_bits64_64 = high_bits<std::uint64_t,64>;

    // doesn't compile: EXPECT_EQ(high_bits64_0, 0x00000000);
    EXPECT_EQ(high_bits64_1, 0x8000000000000000);
    EXPECT_EQ(high_bits64_2, 0xC000000000000000);
    EXPECT_EQ(high_bits64_3, 0xE000000000000000);
    EXPECT_EQ(high_bits64_4, 0xF000000000000000);
    EXPECT_EQ(high_bits64_61, 0xFFFFFFFFFFFFFFF8);
    EXPECT_EQ(high_bits64_62, 0xFFFFFFFFFFFFFFFC);
    EXPECT_EQ(high_bits64_63, 0xFFFFFFFFFFFFFFFE);
    EXPECT_EQ(high_bits64_64, 0xFFFFFFFFFFFFFFFF);
}

TEST(bitfiddle_test, test_signed_shr) {
    EXPECT_EQ(signed_shr(1, 1), 0);
    EXPECT_EQ(signed_shr(-1, 1), -1);
    EXPECT_EQ(signed_shr(-1L, 1), -1L);
    EXPECT_EQ(signed_shr(0x80000000, 1), 0xC0000000);
    EXPECT_EQ(signed_shr(0x80000000L, 1), 0x40000000);
    EXPECT_EQ(signed_shr(0x8000000000000000L, 1), 0xC000000000000000L);
    EXPECT_EQ(signed_shr(0x8888888800000000L, 1), 0xC444444400000000L);
    EXPECT_EQ(signed_shr(0x4000000000000000L, 1), 0x2000000000000000L);
}

TEST(bitfiddle_test, test_flush_hibit) {
    EXPECT_EQ(flush_hibit(0), 0);
    EXPECT_EQ(flush_hibit(0xF0), 0xF0);
    EXPECT_EQ(flush_hibit(0x80000000), 0xFFFFFFFF);
    EXPECT_EQ(flush_hibit(0xF0000000L), 0xF0000000L);
    EXPECT_EQ(flush_hibit(0xA000000000000000L), 0xFFFFFFFFFFFFFFFFL);
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
