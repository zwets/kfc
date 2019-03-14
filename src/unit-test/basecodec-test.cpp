/* basecodec-test.cpp
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
#include "basecodec.h"

using namespace kfc;

namespace {

const unsigned A_VAL = 0, C_VAL = 1, G_VAL = 2, T_VAL = 3;
const unsigned INVAL = -1;

static std::uint32_t (*encode_base32)(unsigned char) = encode_base<std::uint32_t,INVAL>;
static std::uint64_t (*encode_base64)(unsigned char) = encode_base<std::uint64_t,INVAL>;
static char (*decode_base32)(std::uint32_t) = decode_base<std::uint32_t>;
static char (*decode_base64)(std::uint64_t) = decode_base<std::uint64_t>;

// decode_base --------------------------------------------------------

TEST(basecodec_test, test_encode_a) {
    EXPECT_EQ(encode_base32('a'), A_VAL);
    EXPECT_EQ(encode_base32('A'), A_VAL);
    EXPECT_EQ(encode_base64('a'), A_VAL);
    EXPECT_EQ(encode_base64('A'), A_VAL);
}

TEST(basecodec_test, test_encode_c) {
    EXPECT_EQ(encode_base32('c'), C_VAL);
    EXPECT_EQ(encode_base32('C'), C_VAL);
    EXPECT_EQ(encode_base64('c'), C_VAL);
    EXPECT_EQ(encode_base64('C'), C_VAL);
}

TEST(basecodec_test, test_encode_g) {
    EXPECT_EQ(encode_base32('g'), G_VAL);
    EXPECT_EQ(encode_base32('G'), G_VAL);
    EXPECT_EQ(encode_base64('g'), G_VAL);
    EXPECT_EQ(encode_base64('G'), G_VAL);
}

TEST(basecodec_test, test_encode_t) {
    EXPECT_EQ(encode_base32('t'), T_VAL);
    EXPECT_EQ(encode_base32('T'), T_VAL);
    EXPECT_EQ(encode_base64('t'), T_VAL);
    EXPECT_EQ(encode_base64('T'), T_VAL);
}

TEST(basecodec_test, test_encode_non_base) {
    unsigned char c = 0;
    do {
        switch (c) { 
            case 'a': case 'c': case 'g': case 't':
            case 'A': case 'C': case 'G': case 'T':
                break;
            default:
                EXPECT_EQ(encode_base32(c), INVAL);
                EXPECT_EQ(encode_base64(c), INVAL);
        }
    } while (++c != 0);
}

TEST(basecodec_test, test_encode_out_of_bounds) {
    EXPECT_EQ(encode_base32(int(-2)), INVAL);
    EXPECT_EQ(encode_base64(int(-2)), INVAL);
}

TEST(basecodec_test, test_encode_signed_char) {
    EXPECT_EQ(encode_base32(int(-2)), INVAL);
    EXPECT_EQ(encode_base64(int(-2)), INVAL);
}

// decode_base --------------------------------------------------------

TEST(basecodec_test, test_decode_a) {
    EXPECT_EQ(decode_base32(A_VAL), 'a');
    EXPECT_EQ(decode_base64(A_VAL), 'a');
}

TEST(basecodec_test, test_decode_c) {
    EXPECT_EQ(decode_base32(C_VAL), 'c');
    EXPECT_EQ(decode_base64(C_VAL), 'c');
}

TEST(basecodec_test, test_decode_g) {
    EXPECT_EQ(decode_base32(G_VAL), 'g');
    EXPECT_EQ(decode_base64(G_VAL), 'g');
}

TEST(basecodec_test, test_decode_t) {
    EXPECT_EQ(decode_base32(T_VAL), 't');
    EXPECT_EQ(decode_base64(T_VAL), 't');
}

TEST(basecodec_test, test_decode_oob_a) {
    EXPECT_EQ(decode_base32(4 + A_VAL), 'a');
    EXPECT_EQ(decode_base64(4 + A_VAL), 'a');
}

TEST(basecodec_test, test_decode_oob_c) {
    EXPECT_EQ(decode_base32(4 + C_VAL), 'c');
    EXPECT_EQ(decode_base64(4 + C_VAL), 'c');
}

TEST(basecodec_test, test_decode_oob_g) {
    EXPECT_EQ(decode_base32(4 + G_VAL), 'g');
    EXPECT_EQ(decode_base64(4 + G_VAL), 'g');
}

TEST(basecodec_test, test_decode_oob_t) {
    EXPECT_EQ(decode_base32(4 + T_VAL), 't');
    EXPECT_EQ(decode_base64(4 + T_VAL), 't');
}

const std::uint32_t HIBIT_32 = ((uint32_t)1) << 31;
const std::uint64_t HIBIT_64 = ((uint64_t)1) << 63;

TEST(basecodec_test, test_decode_neg_a) {
    EXPECT_EQ(decode_base32(HIBIT_32|A_VAL), 'a');
    EXPECT_EQ(decode_base64(HIBIT_64|A_VAL), 'a');
}

TEST(basecodec_test, test_decode_neg_c) {
    EXPECT_EQ(decode_base32(HIBIT_32|C_VAL), 'c');
    EXPECT_EQ(decode_base64(HIBIT_64|C_VAL), 'c');
}

TEST(basecodec_test, test_decode_neg_g) {
    EXPECT_EQ(decode_base32(HIBIT_32|G_VAL), 'g');
    EXPECT_EQ(decode_base64(HIBIT_64|G_VAL), 'g');
}

TEST(basecodec_test, test_decode_neg_t) {
    EXPECT_EQ(decode_base32(HIBIT_32|T_VAL), 't');
    EXPECT_EQ(decode_base64(HIBIT_64|T_VAL), 't');
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
