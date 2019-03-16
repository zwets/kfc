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
const std::uint32_t INV32 = high_bit<std::uint32_t>;
const std::uint64_t INV64 = high_bit<std::uint64_t>;


// ss_encode_one
// ss_encode
// ss_revcomp
// ss_decode
//
// ds_encode_one
// ds_encode
// 


// decoding 32 -------------------------------------------------------------

std::string ds_decode32_1(std::uint32_t k, bool rc=false) { return ds_decode<std::uint32_t,1>(k,rc); }
std::string ss_decode32_1(std::uint32_t k, bool rc=false) { return ss_decode<std::uint32_t,1>(k,rc); }

std::string ds_decode32_3(std::uint32_t k, bool rc=false) { return ds_decode<std::uint32_t,3>(k,rc); }
std::string ss_decode32_3(std::uint32_t k, bool rc=false) { return ss_decode<std::uint32_t,3>(k,rc); }

std::string ds_decode32_15(std::uint32_t k, bool rc=false) { return ds_decode<std::uint32_t,15>(k,rc); }
std::string ss_decode32_15(std::uint32_t k, bool rc=false) { return ss_decode<std::uint32_t,15>(k,rc); }

TEST(kmercodec_test, dec32_1_a) {
    EXPECT_EQ(ds_decode32_1(A_VAL), "a");
    EXPECT_EQ(ds_decode32_1(A_VAL,true), "t");
}

TEST(kmercodec_test, dec32_1_c) {
    EXPECT_EQ(ds_decode32_1(C_VAL), "c");
    EXPECT_EQ(ds_decode32_1(C_VAL,true), "g");
}

TEST(kmercodec_test, dec32_1_g1_is_a) {
    EXPECT_EQ(ds_decode32_1(G_VAL), "a");
    EXPECT_EQ(ds_decode32_1(G_VAL,true), "t");
}

TEST(kmercodec_test, dec32_1_t1_is_c) {
    EXPECT_EQ(ds_decode32_1(T_VAL), "c");
    EXPECT_EQ(ds_decode32_1(T_VAL,true), "g");
}

TEST(kmercodec_test, dec32_1_a_invalid) {
    EXPECT_EQ(ds_decode32_1(INV32|A_VAL), "X");
    EXPECT_EQ(ds_decode32_1(INV32|A_VAL,true), "X");
}

TEST(kmercodec_test, dec32_1_g1_invalid) {
    EXPECT_EQ(ds_decode32_1(INV32|G_VAL), "X");
    EXPECT_EQ(ds_decode32_1(INV32|G_VAL,true), "X");
}

TEST(kmercodec_test, dec32_1_a_ss) {
    EXPECT_EQ(ss_decode32_1(A_VAL), "a");
    EXPECT_EQ(ss_decode32_1(A_VAL,true), "t");
}

TEST(kmercodec_test, dec32_1_c_ss) {
    EXPECT_EQ(ss_decode32_1(C_VAL), "c");
    EXPECT_EQ(ss_decode32_1(C_VAL,true), "g");
}

TEST(kmercodec_test, dec32_1_g_ss) {
    EXPECT_EQ(ss_decode32_1(G_VAL), "g");
    EXPECT_EQ(ss_decode32_1(G_VAL,true), "c");
}

TEST(kmercodec_test, dec32_1_t_ss) {
    EXPECT_EQ(ss_decode32_1(T_VAL), "t");
    EXPECT_EQ(ss_decode32_1(T_VAL,true), "a");
}

TEST(kmercodec_test, dec32_1_a_ss_invalid) {
    EXPECT_EQ(ss_decode32_1(INV32|A_VAL), "X");
    EXPECT_EQ(ss_decode32_1(INV32|A_VAL,true), "X");
}

TEST(kmercodec_test, dec32_1_t_ss_invalid) {
    EXPECT_EQ(ss_decode32_1(INV32|T_VAL), "X");
    EXPECT_EQ(ss_decode32_1(INV32|T_VAL,true), "X");
}

TEST(kmercodec_test, dec32_3_acg) {
    EXPECT_EQ(ds_decode32_3(A_VAL<<3|C_VAL<<2|G_VAL), "acg");
    EXPECT_EQ(ds_decode32_3(A_VAL<<3|C_VAL<<2|G_VAL,true), "cgt");
}

TEST(kmercodec_test, dec32_3_acg_invalid) {
    EXPECT_EQ(ds_decode32_3(INV32|A_VAL<<3|C_VAL<<2|G_VAL), "XXX");
    EXPECT_EQ(ds_decode32_3(INV32|A_VAL<<3|C_VAL<<2|G_VAL,true), "XXX");
}

TEST(kmercodec_test, dec32_3_acg_ss) {
    EXPECT_EQ(ss_decode32_3(A_VAL<<3|C_VAL<<2|G_VAL), "acg");
    EXPECT_EQ(ss_decode32_3(A_VAL<<3|C_VAL<<2|G_VAL,true), "cgt");
}

TEST(kmercodec_test, dec32_3_cgt_ss) {
    EXPECT_EQ(ss_decode32_3(C_VAL<<4|G_VAL<<2|T_VAL), "cgt");
    EXPECT_EQ(ss_decode32_3(C_VAL<<4|G_VAL<<2|T_VAL,true), "acg");
}

TEST(kmercodec_test, dec32_3_cgt_ss_invalid) {
    EXPECT_EQ(ss_decode32_3(INV32|C_VAL<<4|G_VAL<<2|T_VAL), "XXX");
    EXPECT_EQ(ss_decode32_3(INV32|C_VAL<<4|G_VAL<<2|T_VAL,true), "XXX");
}

TEST(kmercodec_test, dec32_15_ds) {
    std::uint32_t kmer = 
        T_VAL<<27|A_VAL<<25|C_VAL<<23|C_VAL<<21|C_VAL<<19|T_VAL<<17|G_VAL<<15|
        C_VAL<<14|A_VAL<<12|C_VAL<<10|C_VAL<<8 |C_VAL<<6 |A_VAL<<4 |C_VAL<<2 |G_VAL;
    EXPECT_EQ(ds_decode32_15(kmer),"taccctgcacccacg");
    EXPECT_EQ(ds_decode32_15(kmer,true),"cgtgggtgcagggta");
}

TEST(kmercodec_test, dec32_15_ss) {
    std::uint32_t kmer = 
        G_VAL<<28|A_VAL<<26|T_VAL<<24|G_VAL<<22|G_VAL<<20|T_VAL<<18|C_VAL<<16|
        T_VAL<<14|T_VAL<<12|G_VAL<<10|C_VAL<<8 |C_VAL<<6 |C_VAL<<4 |C_VAL<<2 |G_VAL;
    EXPECT_EQ(ss_decode32_15(kmer), "gatggtcttgccccg");
    EXPECT_EQ(ss_decode32_15(kmer,true), "cggggcaagaccatc");
}

TEST(kmercodec_test, dec32_15_ds_inv) {
    EXPECT_EQ(ds_decode32_15(INV32), "XXXXXXXXXXXXXXX");
    EXPECT_EQ(ds_decode32_15(-1), "XXXXXXXXXXXXXXX");
}

TEST(kmercodec_test, dec32_15_ss_inv) {
    EXPECT_EQ(ss_decode32_15(-1), "XXXXXXXXXXXXXXX");
}

TEST(kmercodec_test, dec32_15_ds_invalid) {
    std::uint32_t kmer = INV32|
        T_VAL<<27|A_VAL<<25|C_VAL<<23|C_VAL<<21|C_VAL<<19|T_VAL<<17|G_VAL<<15|
        C_VAL<<14|A_VAL<<12|C_VAL<<10|C_VAL<<8 |C_VAL<<6 |A_VAL<<4 |C_VAL<<2 |G_VAL;
    EXPECT_EQ(ds_decode32_15(kmer),"XXXXXXXXXXXXXXX");
    EXPECT_EQ(ds_decode32_15(kmer,true),"XXXXXXXXXXXXXXX");
}

TEST(kmercodec_test, dec32_15_ss_invalid) {
    std::uint32_t kmer = INV32|
        G_VAL<<28|A_VAL<<26|T_VAL<<24|G_VAL<<22|G_VAL<<20|T_VAL<<18|C_VAL<<16|
        T_VAL<<14|T_VAL<<12|G_VAL<<10|C_VAL<<8 |C_VAL<<6 |C_VAL<<4 |C_VAL<<2 |G_VAL;
    EXPECT_EQ(ss_decode32_15(kmer), "XXXXXXXXXXXXXXX");
    EXPECT_EQ(ss_decode32_15(kmer,true), "XXXXXXXXXXXXXXX");
}

// decoding 64 -------------------------------------------------------------

std::string ds_decode64_1(std::uint64_t k, bool rc=false) { return ds_decode<std::uint64_t,1>(k,rc); }
std::string ss_decode64_1(std::uint64_t k, bool rc=false) { return ss_decode<std::uint64_t,1>(k,rc); }

std::string ds_decode64_3(std::uint64_t k, bool rc=false) { return ds_decode<std::uint64_t,3>(k,rc); }
std::string ss_decode64_3(std::uint64_t k, bool rc=false) { return ss_decode<std::uint64_t,3>(k,rc); }

std::string ds_decode64_31(std::uint64_t k, bool rc=false) { return ds_decode<std::uint64_t,31>(k,rc); }
std::string ss_decode64_31(std::uint64_t k, bool rc=false) { return ss_decode<std::uint64_t,31>(k,rc); }

TEST(kmercodec_test, dec64_1_a) {
    EXPECT_EQ(ds_decode64_1(A_VAL), "a");
    EXPECT_EQ(ds_decode64_1(A_VAL,true), "t");
}

TEST(kmercodec_test, dec64_1_c) {
    EXPECT_EQ(ds_decode64_1(C_VAL), "c");
    EXPECT_EQ(ds_decode64_1(C_VAL,true), "g");
}

TEST(kmercodec_test, dec64_1_g1_is_a) {
    EXPECT_EQ(ds_decode64_1(G_VAL), "a");
    EXPECT_EQ(ds_decode64_1(G_VAL,true), "t");
}

TEST(kmercodec_test, dec64_1_t1_is_c) {
    EXPECT_EQ(ds_decode64_1(T_VAL), "c");
    EXPECT_EQ(ds_decode64_1(T_VAL,true), "g");
}

TEST(kmercodec_test, dec64_1_a_ss) {
    EXPECT_EQ(ss_decode64_1(A_VAL), "a");
    EXPECT_EQ(ss_decode64_1(A_VAL,true), "t");
}

TEST(kmercodec_test, dec64_1_c_ss) {
    EXPECT_EQ(ss_decode64_1(C_VAL), "c");
    EXPECT_EQ(ss_decode64_1(C_VAL,true), "g");
}

TEST(kmercodec_test, dec64_1_g_ss) {
    EXPECT_EQ(ss_decode64_1(G_VAL), "g");
    EXPECT_EQ(ss_decode64_1(G_VAL,true), "c");
}

TEST(kmercodec_test, dec64_1_t_ss) {
    EXPECT_EQ(ss_decode64_1(T_VAL), "t");
    EXPECT_EQ(ss_decode64_1(T_VAL,true), "a");
}

TEST(kmercodec_test, dec64_3_acg) {
    EXPECT_EQ(ds_decode64_3(A_VAL<<3|C_VAL<<2|G_VAL), "acg");
    EXPECT_EQ(ds_decode64_3(A_VAL<<3|C_VAL<<2|G_VAL,true), "cgt");
}

TEST(kmercodec_test, dec64_3_acg_ss) {
    EXPECT_EQ(ss_decode64_3(A_VAL<<3|C_VAL<<2|G_VAL), "acg");
    EXPECT_EQ(ss_decode64_3(A_VAL<<3|C_VAL<<2|G_VAL,true), "cgt");
}

TEST(kmercodec_test, dec64_3_cgt_ss) {
    EXPECT_EQ(ss_decode64_3(C_VAL<<4|G_VAL<<2|T_VAL), "cgt");
    EXPECT_EQ(ss_decode64_3(C_VAL<<4|G_VAL<<2|T_VAL,true), "acg");
}

TEST(kmercodec_test, dec64_31_ds) {
    std::uint64_t kmer = 
        T_VAL<<28|A_VAL<<26|C_VAL<<24|C_VAL<<22|C_VAL<<20|T_VAL<<18|G_VAL<<16|
        C_VAL<<14|A_VAL<<12|C_VAL<<10|C_VAL<<8 |C_VAL<<6 |A_VAL<<4 |C_VAL<<2 |G_VAL;
    EXPECT_EQ(ds_decode64_31(kmer<<31|A_VAL<<30|kmer),"taccctgcacccacg" "a" "taccctgcacccacg");
    EXPECT_EQ(ds_decode64_31(kmer<<31|A_VAL<<30|kmer,true),"cgtgggtgcagggta" "t" "cgtgggtgcagggta");
}

TEST(kmercodec_test, dec64_31_ss) {
    std::uint64_t kmer = 
        G_VAL<<28|A_VAL<<26|T_VAL<<24|G_VAL<<22|G_VAL<<20|T_VAL<<18|C_VAL<<16|
        T_VAL<<14|T_VAL<<12|G_VAL<<10|C_VAL<<8 |C_VAL<<6 |C_VAL<<4 |C_VAL<<2 |G_VAL;
    EXPECT_EQ(ss_decode64_31(kmer<<32|G_VAL<<30|kmer), "gatggtcttgccccg" "g" "gatggtcttgccccg");
    EXPECT_EQ(ss_decode64_31(kmer<<32|G_VAL<<30|kmer,true), "cggggcaagaccatc" "c" "cggggcaagaccatc");
}

TEST(kmercodec_test, dec64_31_ds_inv) {
    EXPECT_EQ(ds_decode64_31(INV64), "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
    EXPECT_EQ(ds_decode64_31(-1), "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
}

TEST(kmercodec_test, dec64_31_ss_inv) {
    EXPECT_EQ(ss_decode64_31(-1), "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
}

TEST(kmercodec_test, dec64_15_ds_invalid) {
    std::uint64_t kmer = INV64|
        T_VAL<<27|A_VAL<<25|C_VAL<<23|C_VAL<<21|C_VAL<<19|T_VAL<<17|G_VAL<<15|
        C_VAL<<14|A_VAL<<12|C_VAL<<10|C_VAL<<8 |C_VAL<<6 |A_VAL<<4 |C_VAL<<2 |G_VAL;
    EXPECT_EQ(ds_decode64_31(kmer),"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
    EXPECT_EQ(ds_decode64_31(kmer,true),"XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
}

TEST(kmercodec_test, dec64_31_ss_invalid) {
    std::uint64_t kmer = INV64|
        G_VAL<<28|A_VAL<<26|T_VAL<<24|G_VAL<<22|G_VAL<<20|T_VAL<<18|C_VAL<<16|
        T_VAL<<14|T_VAL<<12|G_VAL<<10|C_VAL<<8 |C_VAL<<6 |C_VAL<<4 |C_VAL<<2 |G_VAL;
    EXPECT_EQ(ss_decode64_31(kmer), "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
    EXPECT_EQ(ss_decode64_31(kmer,true), "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX");
}

// ss_to_ds 32 ----------------------------------------------------------

TEST(kmercodec_test, ss_to_ds32) {
    std::uint32_t kmer_ds;
    kmer_ds = ss_to_ds<std::uint32_t,1>(A_VAL);
    EXPECT_EQ(kmer_ds, A_VAL);
    kmer_ds = ss_to_ds<std::uint32_t,1>(C_VAL);
    EXPECT_EQ(kmer_ds, C_VAL);
    kmer_ds = ss_to_ds<std::uint32_t,1>(G_VAL);
    EXPECT_EQ(kmer_ds, C_VAL);
    kmer_ds = ss_to_ds<std::uint32_t,1>(T_VAL);
    EXPECT_EQ(kmer_ds, A_VAL);
    std::uint32_t kmer = 
        G_VAL<<28|A_VAL<<26|T_VAL<<24|G_VAL<<22|G_VAL<<20|T_VAL<<18|C_VAL<<16|
        T_VAL<<14|T_VAL<<12|G_VAL<<10|C_VAL<<8 |C_VAL<<6 |C_VAL<<4 |C_VAL<<2 |G_VAL;
    EXPECT_EQ(ss_decode32_15(kmer), "gatggtcttgccccg");
    kmer_ds = ss_to_ds<std::uint32_t,15>(kmer);
    EXPECT_EQ(ds_decode32_15(kmer_ds), "cggggcaagaccatc");
}

TEST(kmercodec_test, ss_to_ds32_3) {
    std::uint32_t kmer_ds;
    kmer_ds = ss_to_ds<std::uint32_t,3>(A_VAL<<4|C_VAL<<2|G_VAL);
    EXPECT_EQ(kmer_ds, A_VAL<<3|C_VAL<<2|G_VAL);
}

TEST(kmercodec_test, ss_to_ds32_3_rc) {
    std::uint32_t kmer_ds;
    kmer_ds = ss_to_ds<std::uint32_t,3>(C_VAL<<4|G_VAL<<2|T_VAL);
    EXPECT_EQ(kmer_ds, A_VAL<<3|C_VAL<<2|G_VAL);
}

TEST(kmercodec_test, ss_to_ds32_3_xchk) {
    std::uint32_t kmer_ds, kmer_ds_rc;
    kmer_ds    = ss_to_ds<std::uint32_t,3>(C_VAL<<4|T_VAL<<2|A_VAL);
    kmer_ds_rc = ss_to_ds<std::uint32_t,3>(T_VAL<<4|A_VAL<<2|G_VAL);
    EXPECT_EQ(kmer_ds, kmer_ds_rc);
}

TEST(kmercodec_test, ss_to_ds32_inv) {
    std::uint32_t kmer_ds;
    kmer_ds = ss_to_ds<std::uint32_t,3>(INV32|A_VAL<<4|C_VAL<<2|G_VAL);
    EXPECT_EQ(kmer_ds & INV32, INV32);
}

// ss_to_ds 64 -------------------------------------------------------------

TEST(kmercodec_test, ss_to_ds64) {
    std::uint64_t kmer_ds;
    kmer_ds = ss_to_ds<std::uint64_t,1>(A_VAL);
    EXPECT_EQ(kmer_ds, A_VAL);
    kmer_ds = ss_to_ds<std::uint64_t,1>(C_VAL);
    EXPECT_EQ(kmer_ds, C_VAL);
    kmer_ds = ss_to_ds<std::uint64_t,1>(G_VAL);
    EXPECT_EQ(kmer_ds, C_VAL);
    kmer_ds = ss_to_ds<std::uint64_t,1>(T_VAL);
    EXPECT_EQ(kmer_ds, A_VAL);

    std::uint64_t kmer = 
        (std::uint64_t)A_VAL<<60|
        (std::uint64_t)G_VAL<<58|(std::uint64_t)A_VAL<<56|(std::uint64_t)T_VAL<<54|
        (std::uint64_t)G_VAL<<52|(std::uint64_t)G_VAL<<50|(std::uint64_t)T_VAL<<48|
        (std::uint64_t)C_VAL<<46|(std::uint64_t)T_VAL<<44|(std::uint64_t)T_VAL<<42|
        (std::uint64_t)G_VAL<<40|(std::uint64_t)C_VAL<<38|(std::uint64_t)C_VAL<<36|
        (std::uint64_t)C_VAL<<34|(std::uint64_t)C_VAL<<32|
        T_VAL<<30|
        G_VAL<<28|A_VAL<<26|T_VAL<<24|G_VAL<<22|G_VAL<<20|T_VAL<<18|C_VAL<<16|
        T_VAL<<14|T_VAL<<12|G_VAL<<10|C_VAL<<8 |C_VAL<<6 |C_VAL<<4 |C_VAL<<2 |G_VAL;
    EXPECT_EQ(ss_decode64_31(kmer), "a""gatggtcttgcccc""t""gatggtcttgccccg");
    kmer_ds = ss_to_ds<std::uint64_t,31>(kmer);
    EXPECT_EQ(ds_decode64_31(kmer_ds), "cggggcaagaccatc""a""ggggcaagaccatc""t");
}

TEST(kmercodec_test, ss_to_ds64_3) {
    std::uint64_t kmer_ds;
    kmer_ds = ss_to_ds<std::uint64_t,3>(A_VAL<<4|C_VAL<<2|G_VAL);
    EXPECT_EQ(kmer_ds, A_VAL<<3|C_VAL<<2|G_VAL);
}

TEST(kmercodec_test, ss_to_ds64_3_rc) {
    std::uint64_t kmer_ds;
    kmer_ds = ss_to_ds<std::uint64_t,3>(C_VAL<<4|G_VAL<<2|T_VAL);
    EXPECT_EQ(kmer_ds, A_VAL<<3|C_VAL<<2|G_VAL);
}

TEST(kmercodec_test, ss_to_ds64_3_xchk) {
    std::uint64_t kmer_ds, kmer_ds_rc;
    kmer_ds    = ss_to_ds<std::uint64_t,3>(C_VAL<<4|T_VAL<<2|A_VAL);
    kmer_ds_rc = ss_to_ds<std::uint64_t,3>(T_VAL<<4|A_VAL<<2|G_VAL);
    EXPECT_EQ(kmer_ds, kmer_ds_rc);
}

TEST(kmercodec_test, ss_to_ds64_inv) {
    std::uint64_t kmer_ds;
    kmer_ds = ss_to_ds<std::uint64_t,3>(INV64|A_VAL<<4|C_VAL<<2|G_VAL);
    EXPECT_EQ(kmer_ds & INV64, INV64);
}

// ss_revcomp -------------------------------------------------------------

TEST(kmercodec_test, ss_revcomp32_1) {
    std::uint32_t kmer_rc;
    kmer_rc = ss_revcomp<std::uint32_t,1>(T_VAL);
    EXPECT_EQ(kmer_rc, A_VAL);
}

TEST(kmercodec_test, ss_revcomp32_1_inv) {
    std::uint32_t kmer_rc;
    kmer_rc = ss_revcomp<std::uint32_t,1>(INV32|T_VAL);
    EXPECT_EQ(kmer_rc & INV32, INV32);
}

TEST(kmercodec_test, ss_revcomp32_3) {
    std::uint32_t kmer_rc;
    kmer_rc = ss_revcomp<std::uint32_t,3>(T_VAL<<4|G_VAL<<2|C_VAL);
    EXPECT_EQ(kmer_rc, G_VAL<<4|C_VAL<<2|A_VAL);
}

TEST(kmercodec_test, ss_revcomp32_3_inv) {
    std::uint32_t kmer_rc;
    kmer_rc = ss_revcomp<std::uint32_t,3>(INV32|T_VAL<<4|G_VAL<<2|C_VAL);
    EXPECT_EQ(kmer_rc & INV32, INV32);
}

// ss_revcomp 64 -------------------------------------------------------------

TEST(kmercodec_test, ss_revcomp64_1) {
    std::uint64_t kmer_rc;
    kmer_rc = ss_revcomp<std::uint32_t,1>(T_VAL);
    EXPECT_EQ(kmer_rc, A_VAL);
}

TEST(kmercodec_test, ss_revcomp64_1_inv) {
    std::uint64_t kmer_rc;
    kmer_rc = ss_revcomp<std::uint64_t,1>(INV64|T_VAL);
    EXPECT_EQ(kmer_rc & INV64, INV64);
}

TEST(kmercodec_test, ss_revcomp64_3) {
    std::uint64_t kmer_rc;
    kmer_rc = ss_revcomp<std::uint64_t,3>(T_VAL<<4|G_VAL<<2|C_VAL);
    EXPECT_EQ(kmer_rc, G_VAL<<4|C_VAL<<2|A_VAL);
}

TEST(kmercodec_test, ss_revcomp64_3_inv) {
    std::uint64_t kmer_rc;
    kmer_rc = ss_revcomp<std::uint64_t,3>(INV64|T_VAL<<4|G_VAL<<2|C_VAL);
    EXPECT_EQ(kmer_rc & INV64, INV64);
}

// encode_one -------------------------------------------------------------

TEST(kmercodec_test, encode32_one_ksize_3) {
    const char seq[] = "cgt";

    std::uint32_t res = ss_encode_one<std::uint32_t,3>(seq);
    EXPECT_EQ(27,res); // cgt -> 011011

    res = ds_encode_one<std::uint32_t,3>(seq);
    EXPECT_EQ(6,res);  // cgt -> acg -> 000110
}

TEST(kmercodec_test, encode32_fwd_and_rev) {
    const char seq[] = "acgattagcgatagggt";
    const char rev[] = "accctatcgctaatcgt";

    for (unsigned i = 0; i < sizeof(seq)-7; ++i) {
        std::uint32_t kmer1 = ds_encode_one<std::uint32_t,7>(seq+i);
        std::uint32_t kmer2 = ds_encode_one<std::uint32_t,7>(rev+10-i);
        EXPECT_EQ(kmer1,kmer2);
    }
}

TEST(kmercodec_test, encode32_ksize_3) {
    const unsigned ksize = 3;
    const char seq[] = "acgtca";
    std::uint32_t res[sizeof(seq)-ksize];

    ss_encode<std::uint32_t,3>(seq, seq+sizeof(seq)-1, res);
    EXPECT_EQ(6,res[0]);  // acg -> 000110
    EXPECT_EQ(27,res[1]); // cgt -> 011011
    EXPECT_EQ(45,res[2]); // gtc -> 101101
    EXPECT_EQ(52,res[3]); // tca -> 110100

    ds_encode<std::uint32_t,3>(seq, seq+sizeof(seq)-1, res);
    EXPECT_EQ(6,res[0]);  // acg -> 00110
    EXPECT_EQ(6,res[1]);  // cgt -> acg -> 00110
    EXPECT_EQ(17,res[2]); // gtc -> gac -> 10001
    EXPECT_EQ(28,res[3]); // tca -> 11100
}

TEST(kmercodec_test, encode32_ksize_15) {
    char seq[] = "gaatctgcccagcac"; // 10 0000 1101 1110 (0)101 0100 1001 0001
    std::uint32_t r32s[sizeof(seq)-15];
    std::uint32_t r32d[sizeof(seq)-15];
    std::uint64_t r64s[sizeof(seq)-15];
    std::uint64_t r64d[sizeof(seq)-15];
    const std::uint32_t r_ds = 0x106F5491, r_ss = 0x20DE5491;
    
    ss_encode<std::uint32_t,15>(seq, seq+sizeof(seq)-1, r32s);
    EXPECT_EQ(r_ss,r32s[0]);

    ds_encode<std::uint32_t,15>(seq, seq+sizeof(seq)-1, r32d);
    EXPECT_EQ(r_ds,r32d[0]);

    ss_encode<std::uint64_t,15>(seq, seq+sizeof(seq)-1, r64s);
    EXPECT_EQ(r_ss,r64s[0]);

    ds_encode<std::uint64_t,15>(seq, seq+sizeof(seq)-1, r64d);
    EXPECT_EQ(r_ds,r64d[0]);
}

TEST(kmercodec_test, encode_ksize_31) {
    char seq[] = "TAAGCGTTTGCTATGCCATCCCATCGGGCCA"; 
    // 11 0000 1001 1011 1111 1001 1100 1110 (0)101 0011 0101 0100 1101 1010 1001 0100
    std::uint64_t r64s[sizeof(seq)-15];
    std::uint64_t r64d[sizeof(seq)-15];
    const std::uint64_t r_ds = 0x184DFCE75354DA94, r_ss = 0x309BF9CE5354DA94;
    
    ss_encode<std::uint64_t,31>(seq, seq+sizeof(seq)-1, r64s);
    EXPECT_EQ(r_ss,r64s[0]);

    ds_encode<std::uint64_t,31>(seq, seq+sizeof(seq)-1, r64d);
    EXPECT_EQ(r_ds,r64d[0]);
}

// invalid input --------------------------------------------------------

TEST(kmercodec_test, invalid32) {
    char seq[] = "cgn";
    std::uint32_t res[sizeof(seq)-3];
    ds_encode<std::uint32_t,3>(seq, seq+sizeof(seq)-1, res);
    EXPECT_TRUE(res[0] & high_bit<std::uint32_t>);
}

TEST(kmercodec_test, invalid32_pos1) {
    char seq[] = "aaaacn";
    std::uint32_t res[5];
    ds_encode<std::uint32_t,5>(seq, seq+sizeof(seq)-1, res);
    EXPECT_EQ(res[0], 1);
    EXPECT_TRUE(res[1] & high_bit<std::uint32_t>);
}

TEST(kmercodec_test, invalid32_mid) {
    char seq[] = "aaaacngtttt";
    std::uint32_t res[5];
    ds_encode<std::uint32_t,5>(seq, seq+sizeof(seq)-1, res);
    EXPECT_EQ(res[0], 1);
    EXPECT_TRUE(res[1] & high_bit<std::uint32_t>);
    EXPECT_TRUE(res[2] & high_bit<std::uint32_t>);
    EXPECT_TRUE(res[3] & high_bit<std::uint32_t>);
    EXPECT_TRUE(res[4] & high_bit<std::uint32_t>);
    EXPECT_TRUE(res[5] & high_bit<std::uint32_t>);
    EXPECT_EQ(res[6], 1);
}

TEST(kmercodec_test, invalid64) {
    char seq[] = "cgn";
    std::uint64_t res[sizeof(seq)-3];
    ds_encode<std::uint64_t,3>(seq, seq+sizeof(seq)-1, res);
    EXPECT_TRUE(res[0] & high_bit<std::uint64_t>);
}

TEST(kmercodec_test, invalid64_mid) {
    char seq[] = "aaaacngtttt";
    std::uint64_t res[5];
    ds_encode<std::uint64_t,5>(seq, seq+sizeof(seq)-1, res);
    EXPECT_EQ(res[0], 1);
    EXPECT_TRUE(res[1] & high_bit<std::uint64_t>);
    EXPECT_TRUE(res[2] & high_bit<std::uint64_t>);
    EXPECT_TRUE(res[3] & high_bit<std::uint64_t>);
    EXPECT_TRUE(res[4] & high_bit<std::uint64_t>);
    EXPECT_TRUE(res[5] & high_bit<std::uint64_t>);
    EXPECT_EQ(res[6], 1);
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
