/* kmercounter-test.cpp
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
#include "kmercounter.h"
#include "utils.h"

using namespace kfc;

namespace {

typedef kmer_counter_list<std::uint32_t> counter32list; // kmer_t 32, count_t 32
typedef kmer_counter_list<std::uint64_t> counter64list; // kmer_t 32, count_t 32
typedef kmer_counter_tally<std::uint32_t,std::uint32_t> counter32tally; // kmer_t 32, count_t 32
typedef kmer_counter_tally<std::uint64_t,std::uint32_t> counter64tally; // kmer_t 32, count_t 32


// sizes and limits -----------------------------------------------------

TEST(kmercounter_test, no_ksize_zero) {
    EXPECT_DEATH(kmer_counter::create(0), ".*");
}

TEST(kmercounter_test, ksize_1) {
    delete(kmer_counter::create(1));
    delete(kmer_counter::create(1, true));
}

TEST(kmercounter_test, no_ksize_16) {
    EXPECT_DEATH(counter32tally::create(16,false,0,0), ".*");
    EXPECT_DEATH(counter32list::create(16,false,0,0,0), ".*");
}

TEST(kmercounter_test, no_ksize_16_ss) {
    EXPECT_DEATH(counter32tally(16, true,0,0), ".*");
    EXPECT_DEATH(counter32list(16, true,0,0,0), ".*");
}

TEST(kmercounter_test, ksize_15) {
    delete(kmer_counter::create(15));
    delete(kmer_counter::create(15,true));
}

TEST(kmercounter_test, no_ksize_32) {
    EXPECT_DEATH(counter32tally(32,false,0,0), ".*");
    EXPECT_DEATH(counter32list(32,false,0,0,0), ".*");
}

TEST(kmercounter_test, no_ksize_32_ss) {
    EXPECT_DEATH(counter32tally(32,true,0,0), ".*");
    EXPECT_DEATH(counter32list(32,true,0,0,0), ".*");
}

TEST(kmercounter_test, ksize_31) {
    delete(kmer_counter::create(31));
    delete(kmer_counter::create(31, true));
}

TEST(kmercounter_test, no_ksize_even) {
    EXPECT_DEATH(kmer_counter::create(6), ".*");
}

TEST(kmercounter_test, ksize_even_ss) {
    delete(kmer_counter::create(6, true));
}

// decoding -------------------------------------------------------------

static unsigned std_opts = output_opts::no_headers;
static unsigned with_zeros = std_opts | output_opts::zeros;
static unsigned with_invalids = std_opts | output_opts::invalids;

static int result_line_count(kmer_counter& c, std::string s, unsigned o = std_opts)
{
    int count = 0; char buf[10240]; std::stringstream ss;

    c.process(s);
    c.write_results(ss, o);
    while (ss.getline(buf, sizeof(buf)))
        ++count;

    return count;
}

TEST(kmercounter_test, process_1_a) {
    counter32tally c(1, false, 0, 0);
    EXPECT_EQ(result_line_count(c, "a"), 1);
    EXPECT_EQ(result_line_count(c, "a", with_zeros), 2);
}

TEST(kmercounter_test, process_1_a_ss) {
    counter32tally c(1, true, 0, 0);
    EXPECT_EQ(result_line_count(c, "a"), 1);
    EXPECT_EQ(result_line_count(c, "a", with_zeros), 4);
}

TEST(kmercounter_test, process_1_aa) {
    counter32tally c(1, false, 0, 0);
    EXPECT_EQ(result_line_count(c, "aa"), 1);
}

TEST(kmercounter_test, process_1_ag) {
    counter32tally c(1, false, 0, 0);
    EXPECT_EQ(result_line_count(c, "ag"), 2);
    EXPECT_EQ(result_line_count(c, "ag", with_zeros), 2);
}

TEST(kmercounter_test, process_1_at) {
    counter32tally c(1, false, 0, 0);
    EXPECT_EQ(result_line_count(c, "at"), 1);
}

TEST(kmercounter_test, process_1_at_ss) {
    counter32tally c(1, true, 0, 0);
    EXPECT_EQ(result_line_count(c, "at"), 2);
}

TEST(kmercounter_test, process_1_invalid) {
    counter32tally c(1, false, 0, 0);
    EXPECT_EQ(result_line_count(c, "x"), 0);
    EXPECT_EQ(result_line_count(c, "x", with_invalids), 1);
    EXPECT_EQ(result_line_count(c, "x", with_zeros), 2);
    EXPECT_EQ(result_line_count(c, "x", with_zeros|with_invalids), 3);
}

TEST(kmercounter_test, process_1_invalid_ss) {
    counter32tally c(1, true, 0, 0);
    EXPECT_EQ(result_line_count(c, "x"), 0);
    EXPECT_EQ(result_line_count(c, "x", with_invalids), 1);
    EXPECT_EQ(result_line_count(c, "x", with_zeros), 4);
    EXPECT_EQ(result_line_count(c, "x", with_zeros|with_invalids), 5);
}

TEST(kmercounter_test, process_3_acgt) {
    counter32tally c(3, false, 0, 0);
    EXPECT_EQ(result_line_count(c, "acgt"), 1);
    EXPECT_EQ(result_line_count(c, "acgt", with_invalids), 1);
    EXPECT_EQ(result_line_count(c, "acgt", with_zeros), 32);
    EXPECT_EQ(result_line_count(c, "acgt", with_zeros|with_invalids), 33);
}

TEST(kmercounter_test, process_3_acgt_ss) {
    counter32tally c(3, true, 0, 0);
    EXPECT_EQ(result_line_count(c, "acgt"), 2);
    EXPECT_EQ(result_line_count(c, "acgt", with_invalids), 2);
    EXPECT_EQ(result_line_count(c, "acgt", with_zeros), 64);
    EXPECT_EQ(result_line_count(c, "acgt", with_zeros|with_invalids), 65);
}

// encoding -------------------------------------------------------------

TEST(kmercounter_test, ksize_is_string) {
    counter32tally c(7, false, 0, 0);
    EXPECT_EQ(result_line_count(c, "acgtaaa", with_invalids), 1);
    EXPECT_EQ(result_line_count(c, "acgtaaa", with_zeros), 8192);
}

TEST(kmercounter_test, ksize_longer_than_string) {
    counter32tally c(7, false, 0, 0);
    EXPECT_EQ(result_line_count(c, "xxxxxx", with_invalids), 0);
    EXPECT_EQ(result_line_count(c, "xxxxxx", with_zeros), 8192);
}

TEST(kmercounter_test, reverse_long) {
    const int k = 7;
    counter32tally c1(k, false, 0, 0);
    counter32tally c2(k, false, 0, 0);

    char seq[] = "acgattagcgatagggt";
    char rev[] = "accctatcgctaatcgt";

    c1.process(seq);
    c2.process(rev);

    std::stringstream r1, r2;
    c1.write_results(r1);
    c2.write_results(r2);

    ASSERT_EQ(r1.str(), r2.str());
}

TEST(kmercounter_test, encode_ksize_15) {
    char seq[] = "gaatctgcccagcac";
    kmer_counter *c = kmer_counter::create(15, false);
    c->process(seq);

    std::stringstream ss;
    c->write_results(ss, output_opts::no_headers);
    delete c;

    std::string dna;
    std::uint64_t ccode;
    std::uint64_t count;
    ss >> dna >> ccode >> count;

    ASSERT_EQ(dna, seq);
    ASSERT_EQ(ccode, 0x106F5491);
    ASSERT_EQ(count, 1);
}

TEST(kmercounter_test, encode_ksize_15_ss) {
    char seq[] = "gaatctgcccagcac";
    kmer_counter *c = kmer_counter::create(15, true);
    c->process(seq);

    std::stringstream ss;
    c->write_results(ss, output_opts::no_headers);
    delete c;

    std::string dna;
    std::uint64_t ccode;
    std::uint64_t count;
    ss >> dna >> ccode >> count;

    ASSERT_EQ(dna, seq);
    ASSERT_EQ(ccode, 0x20DE5491);
    ASSERT_EQ(count, 1);
}

TEST(kmercounter_test, encode_ksize_15_1gb) {
    char seq[] = "gaatctgcccagcac";
    kmer_counter *c = kmer_counter::create(15, false, 1);
    c->process(seq);

    std::stringstream ss;
    c->write_results(ss, output_opts::no_headers);
    delete c;

    std::string dna;
    std::uint64_t ccode;
    std::uint64_t count;
    ss >> dna >> ccode >> count;

    ASSERT_EQ(dna, seq);
    ASSERT_EQ(ccode, 0x106F5491);
    ASSERT_EQ(count, 1);
}

TEST(kmercounter_test, encode_ksize_15_ss_1gb) {
    char seq[] = "gaatctgcccagcac";
    kmer_counter *c = kmer_counter::create(15, true, 1);
    c->process(seq);

    std::stringstream ss;
    c->write_results(ss, output_opts::no_headers);
    delete c;

    std::string dna;
    std::uint64_t ccode;
    std::uint64_t count;
    ss >> dna >> ccode >> count;

    ASSERT_EQ(dna, seq);
    ASSERT_EQ(ccode, 0x20DE5491);
    ASSERT_EQ(count, 1);
}

TEST(kmercounter_test, encode_ksize_31) {
    char seq[] = "taagcgtttgctatgccatcccatcgggcca"; 
    kmer_counter *c = kmer_counter::create(31, false);
    c->process(seq);

    std::stringstream ss;
    c->write_results(ss, output_opts::no_headers);
    delete c;

    std::string dna;
    std::uint64_t ccode;
    std::uint64_t count;
    ss >> dna >> ccode >> count;

    ASSERT_EQ(dna, seq);
    ASSERT_EQ(ccode, 0x184DFCE75354DA94);
    ASSERT_EQ(count, 1);
}

TEST(kmercounter_test, encode_ksize_31_ss) {
    char seq[] = "taagcgtttgctatgccatcccatcgggcca"; 
    kmer_counter *c = kmer_counter::create(31, true);
    c->process(seq);

    std::stringstream ss;
    c->write_results(ss, output_opts::no_headers);
    delete c;

    std::string dna;
    std::uint64_t ccode;
    std::uint64_t count;
    ss >> dna >> ccode >> count;

    ASSERT_EQ(dna, seq);
    ASSERT_EQ(ccode, 0x309BF9CE5354DA94);
    ASSERT_EQ(count, 1);
}

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

TEST(kmercounter_test, counter32_1kbase_rc) {

    kmer_counter *c1 = kmer_counter::create(15, false, 1);
    kmer_counter *c2 = kmer_counter::create(15, false, 1);

    c1->process(dna);
    c2->process(rc_dna);

    std::stringstream ss1;
    std::stringstream ss2;

    c1->write_results(ss1, output_opts::no_headers);
    c2->write_results(ss2, output_opts::no_headers|output_opts::invalids);
    
    ASSERT_EQ(ss1.str(), ss2.str());

    delete c1;
    delete c2;
}

TEST(kmercounter_test, counter64_1kbase_rc) {

    kmer_counter *c1 = kmer_counter::create(31, false, 1);
    kmer_counter *c2 = kmer_counter::create(31, false, 1);

    c1->process(dna);
    c2->process(rc_dna);

    std::stringstream ss1;
    std::stringstream ss2;

    c1->write_results(ss1, output_opts::no_headers);
    c2->write_results(ss2, output_opts::no_headers|output_opts::invalids);
    
    ASSERT_EQ(ss1.str(), ss2.str());

    delete c1;
    delete c2;
}


} // namespace
  // vim: sts=4:sw=4:ai:si:et
