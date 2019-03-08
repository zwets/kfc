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

#include <gtest/gtest.h>
#include "kmercodec.h"
#include "utils.h"

using namespace kfc;

namespace {

const unsigned A_VAL = 0, C_VAL = 1, G_VAL = 2, T_VAL = 3;

typedef kmer_codec<std::uint32_t> codec32;
typedef std::vector<std::uint32_t> vector32;
typedef std::vector<std::uint32_t>::const_iterator v32iter;
typedef std::vector<std::uint32_t>::const_reverse_iterator r32iter;
static std::uint32_t invalid32 = kmer_codec<std::uint32_t>::invalid_kmer;

typedef kmer_codec<std::uint64_t> codec64;
typedef std::vector<std::uint64_t> vector64;
typedef std::vector<std::uint64_t>::const_iterator v64iter;
typedef std::vector<std::uint64_t>::const_reverse_iterator r64iter;
static std::uint64_t invalid64 = kmer_codec<std::uint64_t>::invalid_kmer;

// sizes and limits -----------------------------------------------------

static codec32* make_codec32(int ksize, bool sstrand = false)
{
    codec32 *p;
//    std::cerr << "make_codec32(" << ksize << "," << sstrand << ")" << std::endl;
//    std::cerr.flush();
    p = new codec32(ksize, sstrand);
//    std::cerr << "c.max_kmer() = " << p->max_kmer() << std::endl;
    return p;
}

TEST(kmercodec_test, no_ksize_zero) {
    EXPECT_DEATH({ codec32 c(0); }, ".*");
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
    EXPECT_EQ(c.decode(A_VAL,true), "t");
}

TEST(kmercodec_test, ksize_dec1_c) {
    codec32 c(1);
    EXPECT_EQ(c.decode(C_VAL), "c");
    EXPECT_EQ(c.decode(C_VAL,true), "g");
}

TEST(kmercodec_test, ksize_dec1_g_is_invalid) {
    codec32 c(1);
    EXPECT_EQ(c.decode(G_VAL), "invalid");
    EXPECT_EQ(c.decode(G_VAL,true), "invalid");
}

TEST(kmercodec_test, ksize_dec1_t_is_invalid) {
    codec32 c(1);
    EXPECT_EQ(c.decode(T_VAL), "invalid");
    EXPECT_EQ(c.decode(T_VAL,true), "invalid");
}

TEST(kmercodec_test, ksize_dec1_a_ss) {
    codec32 c(1,true);
    EXPECT_EQ(c.decode(A_VAL), "a");
    EXPECT_EQ(c.decode(A_VAL,true), "t");
}

TEST(kmercodec_test, ksize_dec1_c_ss) {
    codec32 c(1,true);
    EXPECT_EQ(c.decode(C_VAL), "c");
    EXPECT_EQ(c.decode(C_VAL,true), "g");
}

TEST(kmercodec_test, ksize_dec1_g_ss) {
    codec32 c(1,true);
    EXPECT_EQ(c.decode(G_VAL), "g");
    EXPECT_EQ(c.decode(G_VAL,true), "c");
}

TEST(kmercodec_test, ksize_dec1_t_ss) {
    codec32 c(1,true);
    EXPECT_EQ(c.decode(T_VAL), "t");
    EXPECT_EQ(c.decode(T_VAL,true), "a");
}

TEST(kmercodec_test, ksize_dec3_acg) {
    codec32 c(3);
    EXPECT_EQ(c.decode(A_VAL<<3|C_VAL<<2|G_VAL), "acg");
    EXPECT_EQ(c.decode(A_VAL<<3|C_VAL<<2|G_VAL,true), "cgt");
}

TEST(kmercodec_test, ksize_dec3_acg_ss) {
    codec32 c(3,true);
    EXPECT_EQ(c.decode(A_VAL<<4|C_VAL<<2|G_VAL), "acg");
    EXPECT_EQ(c.decode(A_VAL<<4|C_VAL<<2|G_VAL,true), "cgt");
}

TEST(kmercodec_test, ksize_dec3_cgt_ss) {
    codec32 c(3,true);
    EXPECT_EQ(c.decode(C_VAL<<4|G_VAL<<2|T_VAL), "cgt");
    EXPECT_EQ(c.decode(C_VAL<<4|G_VAL<<2|T_VAL,true), "acg");
}

// encoding -------------------------------------------------------------

TEST(kmercodec_test, ksize_is_string) {
    codec32 c(7);
    vector32 r1 = c.encode("acgtaaa");
    EXPECT_EQ(r1.size(), 1);
}

TEST(kmercodec_test, ksize_longer_than_string) {
    codec32 c(7);
    vector32 r1 = c.encode("acgtaa");
    EXPECT_EQ(r1.size(), 0);
}

TEST(kmercodec_test, reverse_long) {
    const int k = 7;
    codec32 c(k);

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

TEST(kmercodec_test, encode_ksize_3) {
    codec32 c(3);
    char seq[] = "acgtca";
    vector32 r1 = c.encode(seq);
    ASSERT_EQ(r1.size(),4); // acg -> 00110
    EXPECT_EQ(6,r1[0]); // acg -> 00110
    EXPECT_EQ(6,r1[1]); // cgt -> acg -> 00110
    EXPECT_EQ(17,r1[2]); // gtc -> gac -> 10001
    EXPECT_EQ(28,r1[3]); // tca -> 11100
}

TEST(kmercodec_test, encode_ksize_3_ss) {
    codec32 c(3, true);
    char seq[] = "acgtca";
    vector32 r1 = c.encode(seq);
    ASSERT_EQ(r1.size(),4);
    EXPECT_EQ(6,r1[0]);  // acg -> 000110
    EXPECT_EQ(27,r1[1]); // cgt -> 011011
    EXPECT_EQ(45,r1[2]); // gtc -> 101101
    EXPECT_EQ(52,r1[3]); // tca -> 110100
}

TEST(kmercodec_test, encode_ksize_15) {
    char seq[] = "gaatctgcccagcac"; // 10 0000 1101 1110 (0)101 0100 1001 0001
    const std::uint32_t r_canonical = 0x106F5491, r_sstrand = 0x20DE5491;
    
    vector32 r32s = codec32(15,true).encode(seq);
    ASSERT_EQ(r32s.size(),1);
    EXPECT_EQ(r_sstrand,r32s[0]);

    vector32 r32c = codec32(15).encode(seq);
    ASSERT_EQ(r32c.size(),1);
    EXPECT_EQ(r_canonical,r32c[0]);

    vector64 r64s = codec64(15,true).encode(seq);
    ASSERT_EQ(r64s.size(),1);
    EXPECT_EQ(r_sstrand,r64s[0]);

    vector64 r64c = codec64(15).encode(seq);
    ASSERT_EQ(r64c.size(),1);
    EXPECT_EQ(r_canonical,r64c[0]);
}

TEST(kmercodec_test, encode_ksize_31) {
    char seq[] = "TAAGCGTTTGCTATGCCATCCCATCGGGCCA"; 
    // 11 0000 1001 1011 1111 1001 1100 1110 (0)101 0011 0101 0100 1101 1010 1001 0100
    const std::uint64_t r_canonical = 0x184DFCE75354DA94, r_sstrand = 0x309BF9CE5354DA94;
    
    vector64 r64s = codec64(31,true).encode(seq);
    ASSERT_EQ(r64s.size(),1);
    EXPECT_EQ(r_sstrand,r64s[0]);

    vector64 r64c = codec64(31).encode(seq);
    ASSERT_EQ(r64c.size(),1);
    EXPECT_EQ(r_canonical,r64c[0]);
}

// invalid input --------------------------------------------------------

TEST(kmercodec_test, invalid32) {
    codec32 c(3);
    char seq[] = "cgn";
    vector32 r1 = c.encode(seq);
    ASSERT_EQ(r1.size(),1);
    EXPECT_EQ(r1[0], invalid32);
}

TEST(kmercodec_test, invalid32_pos1) {
    vector32 r1 = codec32(5).encode("aaaacn");
    ASSERT_EQ(r1.size(),2);
    EXPECT_EQ(r1[0], 1);
    EXPECT_EQ(r1[1], invalid32);
}

TEST(kmercodec_test, invalid32_mid) {
    vector32 r1 = codec32(5).encode("aaaacngtttt");
    ASSERT_EQ(r1.size(),7);
    EXPECT_EQ(r1[0], 1);
    EXPECT_EQ(r1[1], invalid32);
    EXPECT_EQ(r1[2], invalid32);
    EXPECT_EQ(r1[3], invalid32);
    EXPECT_EQ(r1[4], invalid32);
    EXPECT_EQ(r1[5], invalid32);
    EXPECT_EQ(r1[6], 1);
}

TEST(kmercodec_test, invalid64) {
    codec64 c(3);
    char seq[] = "cgn";
    std::vector<std::uint64_t> r1 = c.encode(seq);
    EXPECT_EQ(r1.size(),1);
    EXPECT_EQ(r1[0], invalid64);
}

TEST(kmercodec_test, invalid64_mid) {
    vector64 r1 = codec64(5).encode("aaaacngtttt");
    ASSERT_EQ(r1.size(),7);
    EXPECT_EQ(r1[0], 1);
    EXPECT_EQ(r1[1], invalid64);
    EXPECT_EQ(r1[2], invalid64);
    EXPECT_EQ(r1[3], invalid64);
    EXPECT_EQ(r1[4], invalid64);
    EXPECT_EQ(r1[5], invalid64);
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

TEST(kmercodec_test, encodec32_1kbase) {

    for (int ksize = 1; ksize <= codec32::max_ksize; ksize += 2) {   // loop over all possible k-sizes (canonical)
        codec32 c(ksize);
        vector32 r = c.encode(dna);
        ASSERT_EQ(r.size(), KBASE-ksize+1);

        for (const char *p = dna; p != dna + KBASE - ksize + 1; ++p) {
            std::string s(p, p+ksize);
            char m = p[ksize/2];
            ASSERT_EQ(s, c.decode(r[p-dna], m == 'g' || m == 't'));
        }
    }
}

TEST(kmercodec_test, encodec32ss_1kbase) {

    for (int ksize = 1; ksize <= codec32::max_ksize; ++ksize) {  // loop over all k-sizes
        codec32 c(ksize, true);
        vector32 r = c.encode(dna);
        ASSERT_EQ(r.size(), KBASE-ksize+1);

        for (const char *p = dna; p != dna + KBASE - ksize + 1; ++p) {
            std::string s(p, p+ksize);
            ASSERT_EQ(s, c.decode(r[p-dna]));
        }
    }
}

TEST(kmercodec_test, encodec64_1kbase) {

    for (int ksize = 1; ksize <= codec64::max_ksize; ksize += 2) {   // loop over all possible k-sizes (canonical)
        codec64 c(ksize);
        vector64 r = c.encode(dna);
        ASSERT_EQ(r.size(), KBASE-ksize+1);

        for (const char *p = dna; p != dna + KBASE - ksize + 1; ++p) {
            std::string s(p, p+ksize);
            char m = p[ksize/2];
            ASSERT_EQ(s, c.decode(r[p-dna], m == 'g' || m == 't'));
        }
    }
}

TEST(kmercodec_test, encodec64ss_1kbase) {

    for (int ksize = 1; ksize <= codec64::max_ksize; ++ksize) {  // loop over all k-sizes
        codec64 c(ksize, true);
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

TEST(kmercodec_test, encodec32_1kbase_rc) {

    for (int ksize = 1; ksize <= 16; ksize += 2) {   // loop over all possible k-sizes (canonical)

        codec32 c(ksize);

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

TEST(kmercodec_test, encodec64_1kbase_rc) {

    for (int ksize = 1; ksize <= codec64::max_ksize; ksize += 2) {   // loop over all possible k-sizes (canonical)

        codec64 c(ksize);

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
