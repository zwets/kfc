/* implpicker.h
 *
 * Copyright (C) 2019  Marco van Zwetselaar <io@zwets.it>
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
#ifndef implpicker_h_INCLUDED
#define implpicker_h_INCLUDED

#include "kmercounter.h"

namespace kfc {


// pick_implementation() - factory method producing "optimal" kmer_counter
//
// This method return a kmer_counter instance which is optimal given the set
// of parameters passed to it.  See implpicker.cpp for the decision logic.
//
// The caller gets ownership of the kmer_counter that is returned.
//
extern kmer_counter* 
pick_implementation(
        int ksize,              // k-mer size
        bool s_strand,          // single strand encoding
        unsigned max_mbp,       // maximum number of millions of bases
        unsigned max_gb,        // maximum memory use in GB 
        char force_impl,        // force vector, list, map implementation 
        unsigned n_threads);    // number of threads


// pick_implementation() - factory method producing optimal kmer_counter
//
// This unit defines pick_implementation(), a factory method that produces
// an instance of kmer_counter which is optimal given a set of parameters
// and constraints.
//
// Implementations: Vector, Map, List
//
// We currently have three kmer_counter implementations: two based on tallying
// the encoded k-mers as they are being processed, of which one uses a vector
// of tallies indexed by the k-mer number and one uses a map of kmer->tally
// pairs, and one which does not tally but collects the list of k-mer numbers
// as-is, then sorts this list when results are requested.[1]
//
// Parameters: K and C (and S)
//
// The primary parameters for determining implementation are the k-mer size
// and k-mer count.  K-mer size (K) determines B, the number of bits needed to
// store a k-mer, and Q, the size of the k-mer space, that is the number of
// possible k-mers of size K.  K-mer count (C) determines the size of the
// type needed to store a single count, and, for the list implementation, the
// size of the list needed to store all k-mers.
//
// A minor parameter is S, the boolean for "single-stranded" encoding.  If
// it is NOT set, then this saves half the k-mer space Q (and one bit from B),
// since in double-stranded encoding a k-mer and its reverse complement have
// a single k-mer number representation.
//
// Constraints: M, T
//
// The most significant constraint is available memory, but it must be noted
// that maximising memory use does not imply speed gains.  For instance, the
// intuition that the vector implementation is fastest because of its O(1)
// lookup time (whereas for a map this is O(log C), and the list incurs an
// O(C log C) sorting step), is wrong.  The main cause is that random access
// to a large vector gives many cache misses, whereas memory use in the map
// and list implementations is highly local.
//
// Another constraint is T, the number of threads available.  Multi-threading
// mostly benefits the list implementation.  Though in all cases the encoding
// (from string to k-mer numbers) can take place on separate threads[2], in
// the tallying implementations there is contention for the map or vector, as
// these need mutex access whenever a tally takes place.
//
// Types: kmer_t, count_t
//
// Implementations are parameterised on types kmer_t and count_t, both unsigned
// integral types of 32 or 64 bits (though they could be any bitsize, and
// count_t could conceivably be floating point).
//
// The count_t type determines the size of the counts that can be kept.  It can
// be 32-bit for counts up to 4Gi, and must be 64-bit beyond that (up to 16Ei).
// Memory impact is factor 2 on Q for the vector implementation, factor 2 on C
// for the map, and none for the list (it stores no counts).
//
// The kmer_t type is determined by B which is 2K when S, and 2K-1 when !S (or
// taking S ∊ {0,1}, it is 2K-!S).  We use one bit for reporting invalid k-mers,
// so for 32 bit kmer_t, K <= 15, and for 64-bit kmer_t, K <= 31.  This is true
// regardless of S, as K must be odd when !S.  The memory impact of 64-bit is
// factor 2 on C for the map and list implementations, none for the vector
// implementation.
//
// Magnitudes: Q, C
//
// The highest memory effect by far is that of K in the vector implementation:
// K determines B = 2K-!S, which determines Q = 2^B = 2^(2K-!S).  Q is the
// length of the vector, so storage is (for !S and 32-bit count_t): 2KB at K=5,
// 2MB at K=10, 2GB at K=15, 4TB at K=20 (remember the kmer_t jump to 64-bit at
// K>15).  When S is true or count_t is 64-bit, these sizes double.
//
// In the map and list implementations, K has no independent effect other than
// the memory doubling at K>15.  For both, memory consumption is linear with C.
// For the list implementation is it exactly 4C resp. 8C (at K<=15 resp. K>15).
// At human genome scale this means ~12-24GB.  The map implementation has a
// likely overhead per C of 24, making total size 32C to 40C, depending on the
// sizes of count_t and kmer_t.
//
// Limits: M, L, given K and S
//
// The user-settable memory limit M (specified in GB, below we convert to MB)
// and count limit L (specified in M kmers, below we use log2) impose the 
// following constraints.  Given K (range 1..31) and S (range 0..1):
//
// - if L > 4096 then count_t must be 64-bits
// - if K > 15 then kmer_t must be 64-bits
//
// * vector: very sensitive to K, not at all to L
//   - L < 32: mem(Q) = 4*2^(2K-!S) = 2^(2K+2-!S)
//   - L >=32: mem(Q) = 8*2^(2K-!S) = 2^(2K+3-!S)
// * list: factor 4 on L, possible 2 on K, but cost of sorting
//   - K <=15: mem(C) = 4*C = 2^(2+L)
//   - K > 15: mem(C) = 8*C = 2^(3+L)
// * map: factor 32-40 on L, marginal to K
//   - K <=15, L < 32: mem(C) = 2^(5+L)      = 2^(5+L)
//   - K <=15, L >=32: mem(C) = 2^(5+L) + 4L = 2^(5+L) + 2^(2+L)
//   - K > 15, L < 32: mem(C) = 2^(5+L) + 4L = 2^(5+L) + 2^(2+L)
//   - K > 15, L >=32: mem(C) = 2^(5+L) + 8L = 2^(5+L) + 2^(3+L)
//
// --- FOOTNOTES 
//
// [1] Another conceivable implementation would be one which does not encode
//     the k-mers, but keeps the full input as-is as one or more strings, then
//     sorts their ksize substrings into a separate index.  This would be more
//     generic as the input could be in any alphabet.
// [2] Or even on the GPU.
// [3] A special case would be binary counting: yes/no presence of each k-mer.
//     This could happen in a bit-vector.
//

static kmer_counter*
make_instance(char impl, bool big_kmer, bool big_count, int ks, bool ss, size_t nk, unsigned nt)
{
    typedef std::uint32_t u32;
    typedef std::uint64_t u64;

    size_t kb = 2*ks-(ss?0:1);

    verbose_emit("kmer_counter instance: impl %c, ksize %d, kbits %lu, max_count %lu, kmer_t %d, count_t %d, nt %u", 
            impl, ks, kb, nk, big_kmer?64:32, big_count?64:32, nt) ;

    switch (impl) {
        case 'v':
            return big_kmer 
                ? big_count 
                    ? (kmer_counter*) new kmer_counter_tally<u64,u64>(new tallyman_vec<u64,u64>(kb), ks, ss, nt)
                    : (kmer_counter*) new kmer_counter_tally<u64,u32>(new tallyman_vec<u64,u32>(kb), ks, ss, nt)
                : big_count
                    ? (kmer_counter*) new kmer_counter_tally<u32,u64>(new tallyman_vec<u32,u64>(kb), ks, ss, nt)
                    : (kmer_counter*) new kmer_counter_tally<u32,u32>(new tallyman_vec<u32,u32>(kb), ks, ss, nt);
        case 'm':
            return big_kmer 
                ? big_count 
                    ? (kmer_counter*) new kmer_counter_tally<u64,u64>(new tallyman_map<u64,u64>(kb), ks, ss, nt)
                    : (kmer_counter*) new kmer_counter_tally<u64,u32>(new tallyman_map<u64,u32>(kb), ks, ss, nt)
                : big_count
                    ? (kmer_counter*) new kmer_counter_tally<u32,u64>(new tallyman_map<u32,u64>(kb), ks, ss, nt)
                    : (kmer_counter*) new kmer_counter_tally<u32,u32>(new tallyman_map<u32,u32>(kb), ks, ss, nt);
        case 'l':
            return big_kmer 
                    ? (kmer_counter*) new kmer_counter_list<u64>(ks, ss, nk, nt)
                    : (kmer_counter*) new kmer_counter_list<u32>(ks, ss, nk, nt);
        default: 
            raise_error("invalid implementation option: %c", impl);
            return 0;
    }
}

static size_t
map_entry_size(bool big_kmer, bool big_count)
{
    typedef std::uint32_t u32;
    typedef std::uint64_t u64;

    return big_kmer
        ? big_count 
            ? sizeof(std::map<u64,u64>::value_type) : sizeof(std::map<u64,u32>::value_type)
        : big_count 
            ? sizeof(std::map<u32,u64>::value_type) : sizeof(std::map<u32,u32>::value_type);
}

kmer_counter* 
pick_implementation(int ksize, bool s_strand, unsigned max_mbp, unsigned max_gb, char force_impl, unsigned n_threads)
{
    bool big_kmer = false;
    bool big_count = false;
    size_t max_mb = 0;
    size_t max_count = 0;
    size_t sz_vec = 0, sz_lst = 0, sz_map = 0;
    size_t cap_count_lst = 0, cap_count_map = 0;
    unsigned k_bits = 2 * ksize - (s_strand ? 0 : 1);

        // determine big_kmer (if true, then 64-bit)

    if (ksize < 1) {
        raise_error("invalid k-mer size: %d", ksize);
    }
    else if (ksize <= 15) {
        verbose_emit("k-mer size is %d, k-bits is %u, storing in 32-bit kmer_t", ksize, k_bits);
        big_kmer = false;
    }
    else if (ksize <= 31) {
        verbose_emit("k-mer size is %d, k-bits is %u, storing in 64-bit kmer_t", ksize, k_bits);
        big_kmer = true;
    }
    else {
        raise_error("k-mer size %d is too large (maximum is 31)", ksize);
    }

        // determine big_count (default is no)

    if (max_mbp > (((size_t)1)<<61) / 125000UL) {
        raise_error("requested count capacity (%uM) too large for this implementation");
    }
    else if (max_mbp > (((size_t)1)<<29) / 125000UL) {
        verbose_emit("user-specified max count %uM requires 64-bit count_t", max_mbp);
        max_count = 1000000UL * max_mbp;
        big_count = true;
    }
    else if (max_mbp) {
        verbose_emit("user-specified max count %uM fits in 32-bit count_t", max_mbp);
        max_count = 1000000UL * max_mbp;
        big_count = false;
    }
    else {
        verbose_emit("no user-specified max count; defaulting to 32-bit count_t");
        big_count = false;
    }

        // determine max_mb from max_gb or physical memory

    if (max_gb) {
        verbose_emit("user-specified maximum memory: %uGB", max_gb);
        max_mb = static_cast<size_t>(max_gb) << 10;
    }
    else {
        size_t phy_mb = get_system_memory() >> 20;
        verbose_emit("defaulting max memory to all%s physical memory", phy_mb > 2048 ? " but 2G" : "");
        max_mb = phy_mb > 2048 ? phy_mb - 2048 : phy_mb;
    }

    verbose_emit("available memory: %luMB", max_mb);

        // determine memory consumption of the vec impl (is independent of count)

    sz_vec = (big_count ? 8 : 4) * (1UL << (k_bits > 20 ? k_bits - 20 : 0));
    verbose_emit("vector implementation requires %luMB", sz_vec);

        // determine memory consumption of the other impls

    if (max_mbp) {
        sz_lst = (big_kmer ? 8 : 4) * (max_count >> 20);
        if (sz_lst == 0) sz_lst = 1;
        verbose_emit("list implementation requires %luMB", sz_lst);

        sz_map = map_entry_size(big_kmer,big_count) * (max_count >> 20);
        if (sz_map == 0) sz_map = 1;
        verbose_emit("map implementation requires %luMB", sz_map);

            // if user specified max_mbp AND max_gb, then bail out if nothing fits

        if (max_gb && sz_vec > max_mb && sz_map > max_mb && sz_lst > max_mb) 
            raise_error("no implementation can count %uM k-mers in %uGB memory", max_mbp, max_gb);
    }
    else {
        cap_count_lst = (max_mb / (big_kmer ? 8 : 4)) << 20;
        verbose_emit("cap of count in list implementation: %luM", cap_count_lst >> 20);

        cap_count_map = (max_mb << 20) / map_entry_size(big_kmer, big_count);
        verbose_emit("cap of count in map implementation: %luM", cap_count_map >> 20);
    }

        // we set max_count to the list cap count if it was not set by user

    if (!max_count)
        max_count = cap_count_lst;

        // if user forced implementation, check that it fits and return it

    if (force_impl) {

        if (force_impl == 'v' && max_gb && sz_vec > max_mb)
            raise_error("requested vector implementation does not fit in %uGB memory", max_gb);
        else if (force_impl == 'l' && max_gb && max_mbp && sz_lst > max_mb)
            raise_error("requested list implementation cannot count %luM k-mers in %UGB memory", max_mbp, max_gb);
        else if (force_impl == 'm' && max_gb && max_mbp && sz_map > max_mb)
            raise_error("requested map implementation cannot count %luM k-mers in %UGB memory", max_mbp, max_gb);

        verbose_emit("user-specified kmer_counter implementation: %c", force_impl);
        return make_instance(force_impl, big_kmer, big_count, ksize, s_strand, max_count, n_threads);
    }

        // now we can pick the implementation

    if (sz_vec <= 512) { // if within half a GB, just go for the vector
        verbose_emit("vector implementation small (%luMB), picking it", sz_vec);
        return make_instance('v', big_kmer, big_count, ksize, s_strand, max_count, n_threads);
    }
    else if (sz_lst != 0) { // we know the size the list would have
        if (sz_lst < 512) {
            verbose_emit("list implementation small (%luMB), picking it", sz_lst);
            return make_instance('l', big_kmer, big_count, ksize, s_strand, max_count, n_threads);
        }
        else if (sz_vec < sz_lst) {
            verbose_emit("vector implementation (%luMB) smaller than list (%luMB)", sz_vec, sz_lst);
            if (sz_vec > max_mb)
                emit("expect trashing: insufficient physical memory (%luMB)", max_mb);
            return make_instance('v', big_kmer, big_count, ksize, s_strand, max_count, n_threads);
        }
        else {
            verbose_emit("list implementation (%luMB) smaller than vector (%luMB)", sz_lst, sz_vec);
            if (sz_lst > max_mb)
                emit("expect trashing: insufficient physical memory (%luMB)", max_mb);
            return make_instance('l', big_kmer, big_count, ksize, s_strand, max_count, n_threads);
        }
    }
    else { // we don't know the count size
        emit("info: unknown input size; use option -l to optimise processing speed");

        if (sz_vec < max_mb) { // vec fits but list may be faster, notify user
            verbose_emit("picking vector implementation (%luMB) as it fits memory (%u), and count size is unknown", sz_vec, max_gb);
            return make_instance('v', big_kmer, big_count, ksize, s_strand, max_count, n_threads);
        }
        else { // vec impossible, need to choose between map or list, lets take list and hope the best
            verbose_emit("picking list implementation as vector would exceed memory, and count size is unknown");
            return make_instance('l', big_kmer, big_count, ksize, s_strand, max_count, n_threads);
        }
    }
}


} // namespace kfc

#endif // implpicker_h_INCLUDED
       // vim: sts=4:sw=4:ai:si:et
