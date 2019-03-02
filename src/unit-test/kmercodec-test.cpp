/* kmercodec-test.cpp
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

#include "kmercodec.h"

using namespace kfc;

namespace {

const unsigned A_VAL = 0, C_VAL = 1, G_VAL = 2, T_VAL = 3;

typedef kmer_codec<std::uint32_t> codec32;
static std::uint32_t invalid32 = kmer_codec<std::uint32_t>::invalid_kmer;

typedef kmer_codec<std::uint64_t> codec64;
static std::uint64_t invalid64 = kmer_codec<std::uint64_t>::invalid_kmer;

// sizes and limits -----------------------------------------------------

TEST(kmercodec_test, no_ksize_zero) {
    EXPECT_DEATH(codec32(0), ".*");
}

TEST(kmercodec_test, no_ksize_16) {
    EXPECT_DEATH(codec32(16), ".*");
}

TEST(kmercodec_test, no_ksize_16_ss) {
    EXPECT_DEATH(codec32(16, true), ".*");
}

TEST(kmercodec_test, ksize_15) {
    codec32(15);
    codec32(15, true);
}

TEST(kmercodec_test, no_ksize_32) {
    EXPECT_DEATH(codec64(32), ".*");
}

TEST(kmercodec_test, no_ksize_32_ss) {
    EXPECT_DEATH(codec64(32, true), ".*");
}

TEST(kmercodec_test, ksize_31) {
    codec64(31);
    codec64(31, true);
}

TEST(kmercodec_test, no_ksize_even) {
    EXPECT_DEATH(codec32(6), ".*");
}

TEST(kmercodec_test, ksize_even_ss) {
    codec32(6, true);
}

// decoding -------------------------------------------------------------

TEST(kmercodec_test, ksize_dec1_a) {
    codec32 c(1);
    EXPECT_EQ(c.decode(A_VAL), "a");
}

TEST(kmercodec_test, ksize_dec1_c) {
    codec32 c(1);
    EXPECT_EQ(c.decode(C_VAL), "c");
}

TEST(kmercodec_test, ksize_dec1_g_is_a) {
    codec32 c(1);
    EXPECT_EQ(c.decode(G_VAL), "a");
}

TEST(kmercodec_test, ksize_dec1_t_is_c) {
    codec32 c(1);
    EXPECT_EQ(c.decode(T_VAL), "c");
}

TEST(kmercodec_test, ksize_dec1_a_ss) {
    codec32 c(1,true);
    EXPECT_EQ(c.decode(A_VAL), "a");
}

TEST(kmercodec_test, ksize_dec1_c_ss) {
    codec32 c(1,true);
    EXPECT_EQ(c.decode(C_VAL), "c");
}

TEST(kmercodec_test, ksize_dec1_g_ss) {
    codec32 c(1,true);
    EXPECT_EQ(c.decode(G_VAL), "g");
}

TEST(kmercodec_test, ksize_dec1_t_ss) {
    codec32 c(1,true);
    EXPECT_EQ(c.decode(T_VAL), "t");
}

TEST(kmercodec_test, ksize_dec3_acg) {
    codec32 c(3);
    EXPECT_EQ(c.decode(A_VAL<<3|C_VAL<<2|G_VAL), "acg");
}

TEST(kmercodec_test, ksize_dec3_acg_ss) {
    codec32 c(3,true);
    EXPECT_EQ(c.decode(A_VAL<<4|C_VAL<<2|G_VAL), "acg");
}

TEST(kmercodec_test, ksize_dec3_cgt_ss) {
    codec32 c(3,true);
    EXPECT_EQ(c.decode(C_VAL<<4|G_VAL<<2|T_VAL), "cgt");
}

// encoding -------------------------------------------------------------

TEST(kmercodec_test, reverse_long) {
    const int k = 7;
    codec32 c(k);
    char seq[] = "acgattagcgatagggt";
    char rev[] = "accctatcgctaatcgt";

    std::vector<std::uint32_t> r1 = c.encode(seq);
    std::vector<std::uint32_t> r2 = c.encode(rev);

    EXPECT_EQ(r1.size(), sizeof(seq)-7);
    EXPECT_EQ(r2.size(), sizeof(rev)-7);

    for (size_t i = 0; i != r1.size(); ++i)
        EXPECT_EQ(r1[i], r2[r2.size()-1-i]);
}

TEST(kmercodec_test, ksize_3) {
    codec32 c(3);
    char seq[] = "acgtca";
    std::vector<std::uint32_t> r1 = c.encode(seq);
    EXPECT_EQ(r1.size(),4); // acg -> 00110
    EXPECT_EQ(6,r1[0]); // acg -> 00110
    EXPECT_EQ(6,r1[1]); // cgt -> acg -> 00110
    EXPECT_EQ(17,r1[2]); // gtc -> gac -> 10001
    EXPECT_EQ(28,r1[3]); // tca -> 11100
}

// invalid input --------------------------------------------------------

TEST(kmercodec_test, invalid32) {
    codec32 c(3);
    char seq[] = "cgn";
    std::vector<std::uint32_t> r1 = c.encode(seq);
    EXPECT_EQ(r1.size(),1);
    EXPECT_EQ(r1[0], invalid32);
}


TEST(kmercodec_test, invalid64) {
    codec64 c(3);
    char seq[] = "cgn";
    std::vector<std::uint64_t> r1 = c.encode(seq);
    EXPECT_EQ(r1.size(),1);
    EXPECT_EQ(r1[0], invalid64);
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
