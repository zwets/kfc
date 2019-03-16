/* kmerencoder-test.cpp
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
#include "kmerencoder.h"
#include "utils.h"

using namespace kfc;

namespace {

const unsigned A_VAL = 0, C_VAL = 1, G_VAL = 2, T_VAL = 3;

typedef kmer_encoder<std::uint32_t> encoder32;
typedef std::vector<std::uint32_t> vector32;
typedef std::vector<std::uint32_t>::const_iterator v32iter;
typedef std::vector<std::uint32_t>::const_reverse_iterator r32iter;

typedef kmer_encoder<std::uint64_t> encoder64;
typedef std::vector<std::uint64_t> vector64;
typedef std::vector<std::uint64_t>::const_iterator v64iter;
typedef std::vector<std::uint64_t>::const_reverse_iterator r64iter;

// sizes and limits -----------------------------------------------------

TEST(kmerencoder_test, no_ksize_zero) {
    EXPECT_DEATH({ encoder32 c(0); }, ".*");
}

TEST(kmerencoder_test, no_ksize_16) {
    EXPECT_DEATH(encoder32(16), ".*");
}

TEST(kmerencoder_test, no_ksize_16_ss) {
    EXPECT_DEATH(encoder32(16, true), ".*");
}

TEST(kmerencoder_test, ksize_15) {
    encoder32(15);
    encoder32(15, true);
}

TEST(kmerencoder_test, no_ksize_32) {
    EXPECT_DEATH(encoder64(32), ".*");
}

TEST(kmerencoder_test, no_ksize_32_ss) {
    EXPECT_DEATH(encoder64(32, true), ".*");
}

TEST(kmerencoder_test, ksize_31) {
    encoder64(31);
    encoder64(31, true);
}

TEST(kmerencoder_test, no_ksize_even) {
    EXPECT_DEATH(encoder32(6), ".*");
}

TEST(kmerencoder_test, ksize_even_ss) {
    encoder32(6, true);
}

// decoding -------------------------------------------------------------

TEST(kmerencoder_test, ksize_dec1_a) {
    encoder32 c(1);
    EXPECT_EQ(c.decode(A_VAL), "a");
    EXPECT_EQ(c.decode(A_VAL,true), "t");
}

TEST(kmerencoder_test, ksize_dec1_c) {
    encoder32 c(1);
    EXPECT_EQ(c.decode(C_VAL), "c");
    EXPECT_EQ(c.decode(C_VAL,true), "g");
}

TEST(kmerencoder_test, ksize_dec1_g_is_a) {
    encoder32 c(1);
    EXPECT_EQ(c.decode(G_VAL), "a");
    EXPECT_EQ(c.decode(G_VAL,true), "t");
}

TEST(kmerencoder_test, ksize_dec1_t_is_c) {
    encoder32 c(1);
    EXPECT_EQ(c.decode(T_VAL), "c");
    EXPECT_EQ(c.decode(T_VAL,true), "g");
}

TEST(kmerencoder_test, ksize_dec1_a_ss) {
    encoder32 c(1,true);
    EXPECT_EQ(c.decode(A_VAL), "a");
    EXPECT_EQ(c.decode(A_VAL,true), "t");
}

TEST(kmerencoder_test, ksize_dec1_c_ss) {
    encoder32 c(1,true);
    EXPECT_EQ(c.decode(C_VAL), "c");
    EXPECT_EQ(c.decode(C_VAL,true), "g");
}

TEST(kmerencoder_test, ksize_dec1_g_ss) {
    encoder32 c(1,true);
    EXPECT_EQ(c.decode(G_VAL), "g");
    EXPECT_EQ(c.decode(G_VAL,true), "c");
}

TEST(kmerencoder_test, ksize_dec1_t_ss) {
    encoder32 c(1,true);
    EXPECT_EQ(c.decode(T_VAL), "t");
    EXPECT_EQ(c.decode(T_VAL,true), "a");
}

TEST(kmerencoder_test, ksize_dec3_acg) {
    encoder32 c(3);
    EXPECT_EQ(c.decode(A_VAL<<3|C_VAL<<2|G_VAL), "acg");
    EXPECT_EQ(c.decode(A_VAL<<3|C_VAL<<2|G_VAL,true), "cgt");
}

TEST(kmerencoder_test, ksize_dec3_acg_ss) {
    encoder32 c(3,true);
    EXPECT_EQ(c.decode(A_VAL<<4|C_VAL<<2|G_VAL), "acg");
    EXPECT_EQ(c.decode(A_VAL<<4|C_VAL<<2|G_VAL,true), "cgt");
}

TEST(kmerencoder_test, ksize_dec3_cgt_ss) {
    encoder32 c(3,true);
    EXPECT_EQ(c.decode(C_VAL<<4|G_VAL<<2|T_VAL), "cgt");
    EXPECT_EQ(c.decode(C_VAL<<4|G_VAL<<2|T_VAL,true), "acg");
}

// encoding -------------------------------------------------------------

TEST(kmerencoder_test, ksize_is_string) {
    encoder32 c(7);
    vector32 r1 = c.encode("acgtaaa");
    EXPECT_EQ(r1.size(), 1);
}

TEST(kmerencoder_test, ksize_longer_than_string) {
    encoder32 c(7);
    vector32 r1 = c.encode("acgtaa");
    EXPECT_EQ(r1.size(), 0);
}

TEST(kmerencoder_test, reverse_long) {
    const int k = 7;
    encoder32 c(k);

    char seq[] = "acgattagcgatagggt";
    char rev[] = "accctatcgctaatcgt";

    vector32 r1 = c.encode(seq);
    vector32 r2 = c.encode(rev);

    ASSERT_EQ(r1.size(), sizeof(seq)-k);
    ASSERT_EQ(r2.size(), sizeof(rev)-k);
    v32iter p1 = r1.begin();
    r32iter p2 = r2.rbegin();
    for (; p1 != r1.end() && p2 != r2.rend(); ++p1, ++p2)
        EXPECT_EQ(*p1, *p2);
}

TEST(kmerencoder_test, encode_ksize_3) {
    encoder32 c(3);
    char seq[] = "acgtca";
    vector32 r1 = c.encode(seq);
    ASSERT_EQ(r1.size(),4); // acg -> 00110
    EXPECT_EQ(6,r1[0]); // acg -> 00110
    EXPECT_EQ(6,r1[1]); // cgt -> acg -> 00110
    EXPECT_EQ(17,r1[2]); // gtc -> gac -> 10001
    EXPECT_EQ(28,r1[3]); // tca -> 11100
}

TEST(kmerencoder_test, encode_ksize_3_ss) {
    encoder32 c(3, true);
    char seq[] = "acgtca";
    vector32 r1 = c.encode(seq);
    ASSERT_EQ(r1.size(),4);
    EXPECT_EQ(6,r1[0]);  // acg -> 000110
    EXPECT_EQ(27,r1[1]); // cgt -> 011011
    EXPECT_EQ(45,r1[2]); // gtc -> 101101
    EXPECT_EQ(52,r1[3]); // tca -> 110100
}

TEST(kmerencoder_test, encode_ksize_15) {
    char seq[] = "gaatctgcccagcac"; // 10 0000 1101 1110 (0)101 0100 1001 0001
    const std::uint32_t r_canonical = 0x106F5491, r_sstrand = 0x20DE5491;
    
    vector32 r32s = encoder32(15,true).encode(seq);
    ASSERT_EQ(r32s.size(),1);
    EXPECT_EQ(r_sstrand,r32s[0]);

    vector32 r32c = encoder32(15).encode(seq);
    ASSERT_EQ(r32c.size(),1);
    EXPECT_EQ(r_canonical,r32c[0]);

    vector64 r64s = encoder64(15,true).encode(seq);
    ASSERT_EQ(r64s.size(),1);
    EXPECT_EQ(r_sstrand,r64s[0]);

    vector64 r64c = encoder64(15).encode(seq);
    ASSERT_EQ(r64c.size(),1);
    EXPECT_EQ(r_canonical,r64c[0]);
}

TEST(kmerencoder_test, encode_ksize_31) {
    char seq[] = "TAAGCGTTTGCTATGCCATCCCATCGGGCCA"; 
    // 11 0000 1001 1011 1111 1001 1100 1110 (0)101 0011 0101 0100 1101 1010 1001 0100
    const std::uint64_t r_canonical = 0x184DFCE75354DA94, r_sstrand = 0x309BF9CE5354DA94;
    
    vector64 r64s = encoder64(31,true).encode(seq);
    ASSERT_EQ(r64s.size(),1);
    EXPECT_EQ(r_sstrand,r64s[0]);

    vector64 r64c = encoder64(31).encode(seq);
    ASSERT_EQ(r64c.size(),1);
    EXPECT_EQ(r_canonical,r64c[0]);
}

// invalid input --------------------------------------------------------

TEST(kmerencoder_test, invalid32) {
    encoder32 c(3);
    char seq[] = "cgn";
    vector32 r1 = c.encode(seq);
    ASSERT_EQ(r1.size(),1);
    EXPECT_TRUE(c.is_invalid(r1[0]));
}

TEST(kmerencoder_test, invalid32_pos1) {
    vector32 r1 = encoder32(5).encode("aaaacn");
    ASSERT_EQ(r1.size(),2);
    EXPECT_EQ(r1[0], 1);
    EXPECT_TRUE(encoder32::is_invalid(r1[1]));
}

TEST(kmerencoder_test, invalid32_mid) {
    vector32 r1 = encoder32(5).encode("aaaacngtttt");
    ASSERT_EQ(r1.size(),7);
    EXPECT_EQ(r1[0], 1);
    EXPECT_TRUE(encoder32::is_invalid(r1[1]));
    EXPECT_TRUE(encoder32::is_invalid(r1[2]));
    EXPECT_TRUE(encoder32::is_invalid(r1[3]));
    EXPECT_TRUE(encoder32::is_invalid(r1[4]));
    EXPECT_TRUE(encoder32::is_invalid(r1[5]));
    EXPECT_EQ(r1[6], 1);
}

TEST(kmerencoder_test, invalid64) {
    encoder64 c(3);
    char seq[] = "cgn";
    std::vector<std::uint64_t> r1 = c.encode(seq);
    EXPECT_EQ(r1.size(),1);
    EXPECT_TRUE(c.is_invalid(r1[0]));
}

TEST(kmerencoder_test, invalid64_mid) {
    vector64 r1 = encoder64(5).encode("aaaacngtttt");
    ASSERT_EQ(r1.size(),7);
    EXPECT_EQ(r1[0], 1);
    EXPECT_TRUE(encoder64::is_invalid(r1[1]));
    EXPECT_TRUE(encoder64::is_invalid(r1[2]));
    EXPECT_TRUE(encoder64::is_invalid(r1[3]));
    EXPECT_TRUE(encoder64::is_invalid(r1[4]));
    EXPECT_TRUE(encoder64::is_invalid(r1[5]));
    EXPECT_EQ(r1[6], 1);
}

// encoding and decoding ------------------------------------------------

const int KBASE = 1000;
static const char dna[KBASE+1] = 
    "gtcagcataagattttgttacgctttcaatcctggagaacggctggccaggtagtgtgccaaaagagtgctgagtaagcgagagtgcgctggtgtgccta"
    "ttcgaaataagttttgcgagaatgctgtcaccatgggggactcatattgattagacgtatgaaggtttcgagttatttacgataaaagcaacttctgtag"
    "tcagtctcgagttgcgtcgaaccacttgtgcccacgtcattcctacgcgcttttcgcattgcttttccctgggtgtacatgtcgagacaatgtctagagg"
    "agtgcagcataaaccgaggatgtatgtgacccctgacatctattggacaaccctcagccgcggctgcgtttggttgtaccaaactcacctaaaattgcta"
    "tgaatgagcgctaacgttaacgtacttatgctttaggcgtaccgactcaagccattcacctaatcccgaataatgcctccgcacaagtggttttaacagg"
    "tataccgcccggtaacgtagggaactctcccttgtgaggaattgtgagttttaccactttgttccgaaaaatcgtgcagccgcccgccatgcctcgctga"
    "cctcacaatgtgaccggtgcgaagagacgtccggcatcttccccctgactacaaacttcgcattcgctcggtcatctattgcggcttttactacaactaa"
    "taaactatcctaagattctcgagaaatgtagcctttgagaagccgatcggcagtgaatgatatcctcagcatcaaaggctctacgtgctttaaaacgtac"
    "tatccgttatgcgcaaccataacgtcggcgcgagggcccttctgaccactctacccatccggtacctccccatggtgttggtattgggcgacaatgtaat"
    "aaactaccacgcggagcacctgatacggaacatttttgatgtttcctggaccgctcccacgcatggcgcaacagagaaagcagggtattccctgaggctt";

TEST(kmerencoder_test, encoder32_1kbase) {

    for (unsigned ksize = 1; ksize <= encoder32::max_ksize; ksize += 2) {   // loop over all possible k-sizes (canonical)
        encoder32 c(ksize);
        vector32 r = c.encode(dna);
        ASSERT_EQ(r.size(), KBASE-ksize+1);

        for (const char *p = dna; p != dna + KBASE - ksize + 1; ++p) {
            std::string s(p, p+ksize);
            char m = p[ksize/2];
            ASSERT_EQ(s, c.decode(r[p-dna], m == 'g' || m == 't'));
        }
    }
}

TEST(kmerencoder_test, encoder32ss_1kbase) {

    for (unsigned ksize = 1; ksize <= encoder32::max_ksize; ++ksize) {  // loop over all k-sizes
        encoder32 c(ksize, true);
        vector32 r = c.encode(dna);
        ASSERT_EQ(r.size(), KBASE-ksize+1);

        for (const char *p = dna; p != dna + KBASE - ksize + 1; ++p) {
            std::string s(p, p+ksize);
            ASSERT_EQ(s, c.decode(r[p-dna]));
        }
    }
}

TEST(kmerencoder_test, encoder64_1kbase) {

    for (unsigned ksize = 1; ksize <= encoder64::max_ksize; ksize += 2) {   // loop over all possible k-sizes (canonical)
        encoder64 c(ksize);
        vector64 r = c.encode(dna);
        ASSERT_EQ(r.size(), KBASE-ksize+1);

        for (const char *p = dna; p != dna + KBASE - ksize + 1; ++p) {
            std::string s(p, p+ksize);
            char m = p[ksize/2];
            ASSERT_EQ(s, c.decode(r[p-dna], m == 'g' || m == 't'));
        }
    }
}

TEST(kmerencoder_test, encoder64ss_1kbase) {

    for (unsigned ksize = 1; ksize <= encoder64::max_ksize; ++ksize) {  // loop over all k-sizes
        encoder64 c(ksize, true);
        vector64 r = c.encode(dna);
        ASSERT_EQ(r.size(), KBASE-ksize+1);

        for (const char *p = dna; p != dna + KBASE - ksize + 1; ++p) {
            std::string s(p, p+ksize);
            ASSERT_EQ(s, c.decode(r[p-dna]));
        }
    }
}

static const char rc_dna[KBASE+1] =
    "aagcctcagggaataccctgctttctctgttgcgccatgcgtgggagcggtccaggaaacatcaaaaatgttccgtatcaggtgctccgcgtggtagttt"
    "attacattgtcgcccaataccaacaccatggggaggtaccggatgggtagagtggtcagaagggccctcgcgccgacgttatggttgcgcataacggata"
    "gtacgttttaaagcacgtagagcctttgatgctgaggatatcattcactgccgatcggcttctcaaaggctacatttctcgagaatcttaggatagttta"
    "ttagttgtagtaaaagccgcaatagatgaccgagcgaatgcgaagtttgtagtcagggggaagatgccggacgtctcttcgcaccggtcacattgtgagg"
    "tcagcgaggcatggcgggcggctgcacgatttttcggaacaaagtggtaaaactcacaattcctcacaagggagagttccctacgttaccgggcggtata"
    "cctgttaaaaccacttgtgcggaggcattattcgggattaggtgaatggcttgagtcggtacgcctaaagcataagtacgttaacgttagcgctcattca"
    "tagcaattttaggtgagtttggtacaaccaaacgcagccgcggctgagggttgtccaatagatgtcaggggtcacatacatcctcggtttatgctgcact"
    "cctctagacattgtctcgacatgtacacccagggaaaagcaatgcgaaaagcgcgtaggaatgacgtgggcacaagtggttcgacgcaactcgagactga"
    "ctacagaagttgcttttatcgtaaataactcgaaaccttcatacgtctaatcaatatgagtcccccatggtgacagcattctcgcaaaacttatttcgaa"
    "taggcacaccagcgcactctcgcttactcagcactcttttggcacactacctggccagccgttctccaggattgaaagcgtaacaaaatcttatgctgac";

TEST(kmerencoder_test, encoder32_1kbase_rc) {

    for (unsigned ksize = 1; ksize <= 16; ksize += 2) {   // loop over all possible k-sizes (canonical)

        encoder32 c(ksize);

        // encode the forward
        vector32 f = c.encode(dna);
        ASSERT_EQ(f.size(), KBASE-ksize+1);

        // encode the reverse
        vector32 r = c.encode(rc_dna);
        ASSERT_EQ(r.size(), KBASE-ksize+1);

        // check that they're equal
        v32iter fi = f.begin();
        r32iter ri = r.rbegin();
        for (; fi != f.end() && ri != r.rend(); ++fi, ++ri)
            ASSERT_EQ(*fi, *ri);
    }
}

TEST(kmerencoder_test, encoder64_1kbase_rc) {

    for (unsigned ksize = 1; ksize <= encoder64::max_ksize; ksize += 2) {   // loop over all possible k-sizes (canonical)

        encoder64 c(ksize);

        // encode the forward
        vector64 f = c.encode(dna);
        ASSERT_EQ(f.size(), KBASE-ksize+1);

        // encode the reverse
        vector64 r = c.encode(rc_dna);
        ASSERT_EQ(r.size(), KBASE-ksize+1);

        // check that they're equal
        v64iter fi = f.begin();
        r64iter ri = r.rbegin();
        for (; fi != f.end() && ri != r.rend(); ++fi, ++ri)
            ASSERT_EQ(*fi, *ri);
    }
}

} // namespace
  // vim: sts=4:sw=4:ai:si:et
