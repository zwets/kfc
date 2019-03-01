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

typedef kmer_codec<std::uint32_t> kmer_codec32;
typedef kmer_codec<std::uint64_t> kmer_codec64;

TEST(kmercodec_test, no_ksize_zero) {
    EXPECT_DEATH(kmer_codec32(0), ".*");
}

TEST(kmercodec_test, no_ksize_16) {
    EXPECT_DEATH(kmer_codec32(16), ".*");
}

TEST(kmercodec_test, no_ksize_16_ss) {
    EXPECT_DEATH(kmer_codec32(16, true), ".*");
}

TEST(kmercodec_test, ksize_15) {
    kmer_codec32(15);
    kmer_codec32(15, true);
}

TEST(kmercodec_test, no_ksize_32) {
    EXPECT_DEATH(kmer_codec64(32), ".*");
}

TEST(kmercodec_test, no_ksize_32_ss) {
    EXPECT_DEATH(kmer_codec64(32, true), ".*");
}

TEST(kmercodec_test, ksize_31) {
    kmer_codec64(31);
    kmer_codec64(31, true);
}

TEST(kmercodec_test, no_ksize_even) {
    EXPECT_DEATH(kmer_codec32(6), ".*");
}

TEST(kmercodec_test, ksize_even_ss) {
    kmer_codec32(6, true);
}

TEST(kmercodec_test, ksize_1_a) {
    kmer_codec32 c(1);
    EXPECT_EQ(c.decode(A_VAL), "a");
}

TEST(kmercodec_test, ksize_1_c) {
    kmer_codec32 c(1);
    EXPECT_EQ(c.decode(C_VAL), "c");
}

TEST(kmercodec_test, ksize_1_g) {
    kmer_codec32 c(1);
    EXPECT_EQ(c.decode(G_VAL), "a");
}

TEST(kmercodec_test, ksize_1_t) {
    kmer_codec32 c(1);
    EXPECT_EQ(c.decode(T_VAL), "c");
}

TEST(kmercodec_test, ksize_1_a_ss) {
    kmer_codec32 c(1,true);
    EXPECT_EQ(c.decode(A_VAL), "a");
}

TEST(kmercodec_test, ksize_1_c_ss) {
    kmer_codec32 c(1,true);
    EXPECT_EQ(c.decode(C_VAL), "c");
}

TEST(kmercodec_test, ksize_1_g_ss) {
    kmer_codec32 c(1,true);
    EXPECT_EQ(c.decode(G_VAL), "g");
}

TEST(kmercodec_test, ksize_1_t_ss) {
    kmer_codec32 c(1,true);
    EXPECT_EQ(c.decode(T_VAL), "t");
}

TEST(kmercodec_test, ksize_3_acg) {
    kmer_codec32 c(3);
    EXPECT_EQ(c.decode(A_VAL<<3|C_VAL<<2|G_VAL), "acg");
}

TEST(kmercodec_test, ksize_3ss_acg) {
    kmer_codec32 c(3,true);
    EXPECT_EQ(c.decode(A_VAL<<4|C_VAL<<2|G_VAL), "acg");
}

TEST(kmercodec_test, ksize_3ss_cgt) {
    kmer_codec32 c(3,true);
    EXPECT_EQ(c.decode(C_VAL<<4|G_VAL<<2|T_VAL), "cgt");
}

TEST(kmercodec_test, reverse_long) {
    kmer_codec32 c(7);
    char seq[] = "acgattagcgatagggt";
    char rev[] = "accctatcgctaatcgt";

    std::vector<std::uint32_t> r1 = c.encode(seq);
    std::vector<std::uint32_t> r2 = c.encode(rev);

    EXPECT_EQ(r1, r2);
}

TEST(kmercodec_test, ksize_3) {
    kmer_codec32 c(3);
    char seq[] = "acgtca";
    std::vector<std::uint32_t> r1 = c.encode(seq);
    EXPECT_EQ(r1.size(),4); // acg -> 00110
    EXPECT_EQ(6,r1[0]); // acg -> 00110
    EXPECT_EQ(6,r1[1]); // cgt -> acg -> 00110
    EXPECT_EQ(17,r1[2]); // gtc -> gac -> 10001
    EXPECT_EQ(28,r1[3]); // tca -> 11100
}

TEST(kmercodec_test, invalids) {
    kmer_codec32 c(3);
    char seq[] = "cgn";
    std::vector<std::uint32_t> r1 = c.encode(seq);
    EXPECT_EQ(r1.size(),1); // acg -> 00110
    EXPECT_EQ(r1[0], kmer_codec32::invalid_kmer);
}

#if 0
TEST(kmercodec_test, skip_n) {
    kmeriser r(3, true);
    char seq[] = "cgnaaa";
    r.set(seq,seq+strlen(seq));
    EXPECT_TRUE(r.next());
    EXPECT_EQ(0,r.knum());
    EXPECT_FALSE(r.next());
}

TEST(kmercodec_test, skip_all_degen) {
    kmeriser r(3, true);
    char seq[] = "cgbcgdcghcgvcgscgwcgrcgyaaanttt";
    r.set(seq,seq+strlen(seq));
    EXPECT_TRUE(r.next());
    EXPECT_EQ(0,r.knum());
    EXPECT_TRUE(r.next());
    EXPECT_EQ(0,r.knum());
    EXPECT_FALSE(r.next());
}

TEST(kmercodec_test, no_skip_bad) {
    kmeriser r(3, true);
    char seq[] = "cgxaaa";
    r.set(seq,seq+strlen(seq));
    EXPECT_TRUE(r.next());
    EXPECT_DEATH(r.knum(), ".*");
}

TEST(kmercodec_test, same_as_ator) {
    kmeriser ki(5); kmerator ka(5, 1);
    char seq[] = "acgtaaccggttagacatgtacgggattaatag";
    ki.set(seq, seq+strlen(seq));
    ka.set(seq, seq+strlen(seq));
    while (ki.next() && ka.next())
        EXPECT_EQ(ki.knum(), ka.knum());
    EXPECT_FALSE(ki.next());
    EXPECT_FALSE(ka.next());
}

#endif

} // namespace
  // vim: sts=4:sw=4:ai:si:et
