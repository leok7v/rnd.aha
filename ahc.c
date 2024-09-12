// Adaptive Huffman Algorithm R&D
// https://en.wikipedia.org/wiki/Adaptive_Huffman_coding

#ifdef _WIN32
// _WIN32_WINNT_WIN10_* undefined
// since forever...
// https://developercommunity.visualstudio.com/t/several-warnings-in-windows-sdk-100177630-in-windo/435362
#pragma warning(disable: 4668) // replacing with '0' for '#if/#elif'
#include <Windows.h>
#endif

#include "rt.h"

#undef NULL
#undef assert

#define assert(...)  rt_assert(__VA_ARGS__)
#define printf(...)  rt_printf(__VA_ARGS__)

static bool ahc_stats;
static bool ahc_verify;

typedef struct {
    uint64_t freq;
    int32_t  pix;  // parent
    int32_t  lix;  // left
    int32_t  rix;  // right
    int32_t  bits; // 0 for root
    uint64_t path;
} ahc_node_t;

typedef struct {
    int32_t n;
    int32_t depth; // max tree depth seen
    ahc_node_t* node;
    int32_t next;  // next non-terminal nodes in the tree >= n
    int32_t complete; // no more inc freq updates
    // stats:
    struct {
        size_t updates;
        size_t swaps;
        size_t moves;
        size_t each; // report stats on each bytes (e.g. 16 * 1024)
        size_t bits; // encoded bits
    } stats;
} ahc_tree_t;

static int32_t ahc_tree_depth(ahc_tree_t* t) {
    int32_t depth = 0;
    for (int32_t i = 0; i < t->n; i++) {
        if (depth < t->node[i].bits) { depth = t->node[i].bits; }
    }
    assert(depth <= t->depth); // actual depth of the tree
    return depth;
}

static const char* ahc_path2str(uint64_t path, int32_t bits) {
    static char str[64];
    str[bits] = '\0';
    for (int32_t b = 0; b < bits; b++) {
        str[bits - b - 1] = (path & (1ULL << b)) ? '1' : '0';
    }
    return str;
}

static void ahc_print_level(ahc_tree_t* t, int32_t tab, int32_t bits) {
    const int32_t m = t->n * 2 - 1;
    printf("%*c", tab * 6, 0x20);
    if (bits == 0) {
        printf("(R:%03d) ", m - 1); // root
    }
    for (int32_t i = m - 2; i >= 0; i--) {
        if (t->node[i].bits == bits) {
            bool left = t->node[i].pix == -1 || t->node[t->node[i].pix].lix == i;
            printf("(%c:%03d) ", left ? 'L' : 'R', i);
        }
    }
    printf("\n");
    if (tab > 0) { ahc_print_level(t, tab - 1, bits + 1); }
}

static void ahc_print_tree(ahc_tree_t* t) {
    const int32_t m = t->n * 2 - 1;
    const int32_t depth = ahc_tree_depth(t);
    printf("depth: %d\n", depth);
    ahc_print_level(t, depth, 0);
    for (int i = m - 1; i >= 0; i--) {
        printf("[%3d] bits: %3d freq: %6lld \"%s\"\n",
                i, t->node[i].bits, t->node[i].freq,
                ahc_path2str(t->node[i].path, t->node[i].bits));
    }
}

static void ahc_verify_node(ahc_tree_t* t, int32_t ix) {
    const int32_t m = t->n * 2 - 1; (void)m;
    assert(0 <= ix && ix < m);
    const int32_t pix = t->node[ix].pix;
    if (pix == -1) {
        assert(ix == m - 1);
    } else {
        assert(0 <= pix && pix < m);
        assert(t->node[pix].lix == ix || t->node[pix].rix == ix);
    }
    const uint64_t path = t->node[ix].path; (void)path;
    const int32_t  bits = t->node[ix].bits;
    const uint64_t mask = (1ULL << bits) - 1; (void)mask;
    const int32_t  lix = t->node[ix].lix;
    const int32_t  rix = t->node[ix].rix;
    if (0 <= ix && ix < t->n) {
        assert(lix == -1 && rix == -1);
    } else {
        assert(lix >= 0 || rix >= 0);
        if (lix != -1 && rix != -1) {
            assert(0 <= lix && lix < m);
            assert(0 <= rix && rix < m);
            assert(t->node[lix].pix == ix);
            assert(t->node[rix].pix == ix);
            assert(t->node[lix].bits == t->node[rix].bits);
            assert(t->node[lix].bits == bits + 1);
            assert(path == (t->node[lix].path & mask));
            assert(path == (t->node[rix].path & mask));
            assert(t->node[lix].path ==  path);
            assert(t->node[rix].path == (path | (1ULL << bits)));
            assert(t->node[ix].freq == t->node[lix].freq + t->node[rix].freq);
            assert(t->node[lix].freq <= t->node[rix].freq);
            ahc_verify_node(t, lix);
            ahc_verify_node(t, rix);
        } else {
            assert(lix == -1 && rix >= 0);
            assert(0 <= rix && rix < m);
            assert(t->node[rix].pix == ix);
            assert(path == (t->node[rix].path & mask));
            assert(t->node[rix].path == (path | (1ULL << bits)));
            assert(t->node[ix].freq == t->node[rix].freq);
            ahc_verify_node(t, rix);
        }
    }
}

static void ahc_verify_tree(ahc_tree_t* t) {
    if (ahc_verify) {
        const int32_t m = t->n * 2 - 1;
        ahc_verify_node(t, m - 1);
    }
}

static void ahc_update_paths(ahc_tree_t* t, int32_t i) {
    t->stats.updates++;
    const int32_t m = t->n * 2 - 1;
    if (i == m - 1) { t->depth = 0; } // root
    const int32_t  bits = t->node[i].bits;
    const uint64_t path = t->node[i].path;
    assert(bits < (int32_t)sizeof(uint64_t) * 8 - 1);
    assert((path & (~((1ULL << (bits + 1)) - 1))) == 0);
    const int32_t lix = t->node[i].lix;
    const int32_t rix = t->node[i].rix;
    if (lix != -1) {
        t->node[lix].bits = bits + 1;
        t->node[lix].path = path;
        ahc_update_paths(t, lix);
    }
    if (rix != -1) {
        t->node[rix].bits = bits + 1;
        t->node[rix].path = path | (1ULL << bits);
        ahc_update_paths(t, rix);
    }
    if (bits > t->depth) { t->depth = bits; }
}

static int32_t ahc_swap_siblings_if_necessary(ahc_tree_t* t, const int32_t ix) {
    const int32_t m = t->n * 2 - 1;
    assert(0 <= ix && ix < m);
    if (ix < m - 1) { // not root
        const int32_t pix = t->node[ix].pix;
        assert(pix >= t->n); // parent (cannot be a leaf)
        const int32_t lix = t->node[pix].lix;
        const int32_t rix = t->node[pix].rix;
        if (lix >= 0 && rix >= 0) {
            assert(0 <= lix && lix < m - 1 && 0 <= rix && rix < m - 1);
            if (t->node[lix].freq > t->node[rix].freq) { // swap
                t->stats.swaps++;
                t->node[pix].lix = rix;
                t->node[pix].rix = lix;
                ahc_update_paths(t, pix); // because swap changed all path below
                return ix == lix ? rix : lix;
            }
        }
    }
    return ix;
}

static void ahc_frequency_changed(ahc_tree_t* t, int32_t i);

static void ahc_update_freq(ahc_tree_t* t, int32_t i) {
    const int32_t lix = t->node[i].lix;
    const int32_t rix = t->node[i].rix;
    assert(lix != -1 || rix != -1); // at least one leaf present
    t->node[i].freq = (lix >= 0 ? t->node[lix].freq : 0) +
                      (rix >= 0 ? t->node[rix].freq : 0);
}

static void ahc_move_up(ahc_tree_t* t, int32_t i) {
    const int32_t pix = t->node[i].pix; // parent
    assert(pix != -1);
    const int32_t gix = t->node[pix].pix; // grandparent
    assert(gix != -1);
    assert(t->node[pix].rix == i);
    // Is parent grandparent`s left or right child?
    const bool parent_is_left_child = pix == t->node[gix].lix;
    const int32_t psx = parent_is_left_child ? // parent sibling index
        t->node[gix].rix : t->node[gix].lix;   // aka auntie/uncle
    if (t->node[i].freq > t->node[psx].freq) {
        // Move grandparents left or right subtree to be
        // parents right child instead of 'i'.
        t->stats.moves++;
        t->node[i].pix = gix;
        if (parent_is_left_child) {
            t->node[gix].rix = i;
        } else {
            t->node[gix].lix = i;
        }
        t->node[pix].rix = psx;
        t->node[psx].pix = pix;
        ahc_update_freq(t, pix);
        ahc_update_freq(t, gix);
        ahc_swap_siblings_if_necessary(t, i);
        ahc_swap_siblings_if_necessary(t, psx);
        ahc_swap_siblings_if_necessary(t, pix);
        ahc_update_paths(t, gix);
        ahc_frequency_changed(t, gix);
    }
}

static void ahc_frequency_changed(ahc_tree_t* t, int32_t i) {
    const int32_t m = t->n * 2 - 1; (void)m;
    const int32_t pix = t->node[i].pix;
    if (pix == -1) { // `i` is root
        assert(i == m - 1);
        ahc_update_freq(t, i);
        i = ahc_swap_siblings_if_necessary(t, i);
    } else {
        assert(0 <= pix && pix < m);
        ahc_update_freq(t, pix);
        i = ahc_swap_siblings_if_necessary(t, i);
        ahc_frequency_changed(t, pix);
    }
    if (pix != -1 && t->node[pix].pix != -1 && i == t->node[pix].rix) {
        assert(t->node[i].freq >= t->node[t->node[pix].lix].freq);
        ahc_move_up(t, i);
    }
}

static void ahc_insert(ahc_tree_t* t, int32_t i) {
    const int32_t root = t->n * 2 - 1 - 1;
    int32_t ipx = root;
    assert(t->node[i].pix == -1 && t->node[i].lix == -1 && t->node[i].rix == -1);
    assert(t->node[i].freq == 0 && t->node[i].bits == 0 && t->node[i].path == 0);
    t->node[i].freq = 1;
    while (ipx >= t->n) {
        if (t->node[ipx].rix == -1) {
            t->node[ipx].rix = i;
            t->node[i].pix = ipx;
            break;
        } else if (t->node[ipx].lix == -1) {
            t->node[ipx].lix = i;
            t->node[i].pix = ipx;
            break;
        } else {
            assert(t->node[ipx].lix >= 0);
            assert(t->node[i].freq <= t->node[t->node[ipx].lix].freq);
            ipx = t->node[ipx].lix;
        }
    }
    if (ipx >= t->n) { // not a leaf, inserted
        t->node[ipx].freq++;
        ahc_swap_siblings_if_necessary(t, i);
        assert(t->node[ipx].lix == i || t->node[ipx].rix);
        assert(t->node[ipx].freq ==
                (t->node[ipx].rix >= 0 ? t->node[t->node[ipx].rix].freq : 0) +
                (t->node[ipx].lix >= 0 ? t->node[t->node[ipx].lix].freq : 0));
    } else { // leaf
        assert(t->next > t->n);
        if (t->next == t->n) {
            t->complete = true;
        } else {
            t->next--;
            int32_t nix = t->next;
            t->node[nix] = (ahc_node_t){
                .freq = t->node[ipx].freq,
                .lix = ipx,
                .rix = -1,
                .pix = t->node[ipx].pix,
                .bits = t->node[ipx].bits,
                .path = t->node[ipx].path
            };
            if (t->node[ipx].pix != -1) {
                if (t->node[t->node[ipx].pix].lix == ipx) {
                    t->node[t->node[ipx].pix].lix = nix;
                } else {
                    t->node[t->node[ipx].pix].rix = nix;
                }
            }
            t->node[ipx].pix = nix;
            t->node[ipx].bits++;
            t->node[ipx].path = t->node[nix].path;
            t->node[nix].rix = i;
            t->node[i].pix = nix;
            t->node[i].bits = t->node[nix].bits + 1;
            t->node[i].path = t->node[nix].path | (1ULL << t->node[nix].bits);
            ahc_update_freq(t, nix);
            ipx = nix;
        }
    }
    ahc_frequency_changed(t, i);
    ahc_update_paths(t, ipx);
    assert(t->node[i].freq != 0 && t->node[i].bits != 0);
    ahc_verify_tree(t);
}

static void ahc_inc_node_frequency(ahc_tree_t* t, int32_t i) {
    assert(0 <= i && i < t->n); // terminal
    // If input sequence frequencies are severely skewed (e.g. Lucas numbers
    // similar to Fibonacci numbers) and input sequence is long enough.
    // The depth of the tree will grow past 64 bits.
    // The first Lucas number that exceeds 2^64 is
    // L(81) = 18,446,744,073,709,551,616 not actually realistic but
    // better be safe than sorry:
    if (t->node[i].pix == -1) {
        ahc_insert(t, i); // Unseen terminal node.
    } else if (!t->complete && t->depth < 63 && t->node[i].freq < UINT64_MAX - 1) {
        t->node[i].freq++;
        ahc_frequency_changed(t, i);
    } else {
        // ignore future frequency updates
        t->complete = 1;
    }
}

static void ahc_init(ahc_tree_t* t, ahc_node_t nodes[], const int32_t m) {
    assert(m > 3); // must pow(2, bits_per_symbol) * 2 - 1
    const int32_t n = (m + 1) / 2;
    assert(n > 4 && (n & (n - 1)) == 0); // must be power of 2
    memset(&t->stats, 0x00, sizeof(t->stats));
    t->node = nodes;
    t->n = n;
    t->next = m - 1 - 1; // next (after root) non-terminal node
    t->depth = 0;
    t->complete = 0;
    for (int32_t i = 0; i < m; i++) {
        t->node[i] = (ahc_node_t){
            .freq = 0, .pix = -1, .lix = -1, .rix = -1, .bits = 0, .path = 0
        };
    }
}

static uint32_t ahc_random32(uint32_t* state) {
    // https://gist.github.com/tommyettinger/46a874533244883189143505d203312c
    static thread_local bool started; // first seed must be odd
    if (!started) { started = true; *state |= 1; }
    uint32_t z = (*state += 0x6D2B79F5UL);
    z = (z ^ (z >> 15)) * (z | 1UL);
    z ^= z + (z ^ (z >> 7)) * (z | 61UL);
    return z ^ (z >> 14);
}

static double ahc_random(void) { // [0..1[
    static uint32_t seed = 1;
    return (double)ahc_random32(&seed) / ((double)UINT32_MAX + 1);
}

static int32_t ahc_next_random_symbol(uint64_t cf[], uint64_t sum, int32_t n) {
    uint64_t r = (uint64_t)(ahc_random() * sum);
    for (int32_t i = 0; i < n; i++) {
        if (r <= cf[i]) { return i; }
    }
    // fallback, should not occur if freq[] is properly populated
    assert(false);
    return n;
}

static void ahc_shuffle(uint64_t a[], int32_t n) {
    for (int32_t i = n - 1; i > 0; i--) {
        uint64_t j = (uint64_t)(ahc_random() * (i + 1));
        uint64_t swap = a[i]; a[i] = a[j]; a[j] = swap;
//      rt_swap(a[i], a[j]);
    }
}

// Shannon ahc_entropy of frequencies distribution

static double ahc_entropy(uint64_t freq[], int32_t n) {
    double total = 0;
    double ahc_entropy = 0.0;
    for (int32_t i = 0; i < n; i++) { total += (double)freq[i]; }
    for (int32_t i = 0; i < n; i++) {
        if (freq[i] > 0) {
            double p_i = (double)freq[i] / total;
            ahc_entropy += p_i * log2(p_i);
        }
    }
    return -ahc_entropy;
}

static void ahc_generate_geometric(uint64_t freq[], int32_t n,
        double base, int64_t initial_value) {
    freq[0] = initial_value;
    for (int32_t i = 1; i < n; i++) {
        freq[i] = (int64_t)(freq[i-1] * base);
    }
}

static void test(void) {
    enum { bps = 9, n = 1U << bps, m = n * 2 - 1 };
    uint64_t freq[n];
    double base = 1.05;  // Adjust this base to control how quickly values grow
    int64_t initial_value = 100;  // Adjust the initial value to scale the frequencies
    ahc_generate_geometric(freq, n, base, initial_value);
    ahc_shuffle(freq, n);
    uint64_t cf[n] = {0}; // cumulative frequency
    uint64_t total = 0; // total sum of all frequencies
    for (int32_t i = 0; i < n; i++) {
        total += freq[i];
        cf[i] = total;
    }
    ahc_tree_t tree = {0};
    ahc_node_t node[m];
    ahc_tree_t* t = &tree;
    ahc_init(t, node, sizeof(node) / sizeof(node[0]));
    // insert NYE (Not Yet Encoded) marker symbol
    ahc_inc_node_frequency(t, n / 2);
    int32_t  count = 0;
    uint64_t depth_sum = 0;
    uint64_t max_freq  = 2;
    uint64_t bits = 0;
    uint64_t actual[n] = {0};
    enum { N = 10 * 1000 * 1000 };
    assert(t->node != null); // for IntelliSense
    for (int32_t i = 0; i < N; i++) {
        int32_t s = ahc_next_random_symbol(cf, total, n);
        assert(0 <= s && s < n);
        actual[s]++;
        bits += t->node[s].bits;
        ahc_inc_node_frequency(t, s);
        max_freq = t->node[s].freq;
        int32_t depth = ahc_tree_depth(t);
        depth_sum += depth;
        count++;
    }
//  ahc_print_tree(t);
    #ifdef HUFFMAN_DUMP_ACTUAL_FREQUENCIES
        uint64_t total = 0;
        for (int32_t i = 0; i < n; i++) {
            total += actual[i];
        }
        for (int32_t i = 0; i < n; i++) {
            printf("[%3d] %.9f\n", i, (double)actual[i] / (double)total);
        }
    #endif
    printf("\xF0\x9F\x93\xA6 depth max: %d bps: %.1f "
           "Shannon Entropy H: %.3f \n",
           t->depth, (double)bits / (double)N, ahc_entropy(freq, n));
    printf("\n");
}

typedef struct bitstream_s bitstream_t;

typedef struct bitstream_s {
    uint64_t b64;  // bit shifting buffer
    int32_t  bits; // bit count inside b64
    int32_t  padding;
    uint8_t* data; // large for testing purposes
    size_t   capacity;
    size_t   bytes; // number of bytes written
    size_t   read;  // number of bytes read
    errno_t (*write_bit)(bitstream_t* bs, int32_t bit);
    errno_t (*read_bit)(bitstream_t* bs, bool *bit);
} bitstream_t;

static errno_t bitstream_write_bit(bitstream_t* bs, int32_t bit) {
    errno_t r = 0;
    bs->b64 <<= 1;
    bs->b64 |= (bit & 1);
    bs->bits++;
    if (bs->bits == 64) {
        for (int i = 0; i < 8 && r == 0; i++) {
            if (bs->bytes == bs->capacity) { r = E2BIG; }
            bs->data[bs->bytes++] = (bs->b64 >> ((7 - i) * 8)) & 0xFF;
        }
        bs->bits = 0;
        bs->b64 = 0;
    }
    return r;
}

static errno_t bitstream_read_bit(bitstream_t* bs, bool *bit) {
    errno_t r = 0;
    if (bs->bits == 0) {
        bs->b64 = 0;
        for (int i = 0; i < 8 && r == 0; i++) {
            if (bs->read == bs->bytes) {
                r = E2BIG;
            } else {
                const uint64_t byte = (bs->data[bs->read] & 0xFF);
                bs->b64 |= byte << ((7 - i) * 8);
                bs->read++;
            }
        }
        bs->bits = 64;
    }
    bool b = ((int64_t)bs->b64) < 0; // same as (bs->b64 >> 63) & 1;
    bs->b64 <<= 1;
    bs->bits--;
    *bit = (bool)b;
    return r;
}

static errno_t bitstream_create(bitstream_t* bs, size_t capacity) {
    errno_t r = 0;
    memset(bs, 0x00, sizeof(*bs));
    bs->data = (uint8_t*)malloc(capacity);
    if (bs->data != null) {
        bs->capacity  = capacity;
        bs->write_bit = bitstream_write_bit;
        bs->read_bit  = bitstream_read_bit;
    } else {
        r = ENOMEM;
    }
    return r;
}

static void bitstream_dispose(bitstream_t* bs) {
    free(bs->data);
    memset(bs, 0x00, sizeof(*bs));
}

static errno_t ahc_encode(ahc_tree_t* t, const uint8_t data[],
                          size_t bytes, bitstream_t* bs) {
    errno_t r = 0;
    size_t sc = t->stats.swaps;
    size_t mc = t->stats.moves;
    size_t uc = t->stats.updates;
    for (size_t i = 0; r == 0 && i < bytes; i++) {
        int32_t sym = data[i];
        ahc_node_t* s = &t->node[sym];
        if (s->bits == 0) { // first seen, escape with NYE (not yet encoded):
            ahc_node_t* nye = &t->node[t->n / 2];
            for (int32_t b = 0; r == 0 && b < nye->bits; b++) {
                r = bs->write_bit(bs, (nye->path >> b) & 1);
                t->stats.bits++;
            }
//          printf("%02X: ", sym);
            const int32_t bps = rt_const_log2_of_pow2(t->n) - 1;
            for (int32_t b = 0; r == 0 && b < bps; b++) {
                r = bs->write_bit(bs, (sym >> b) & 1);
//              printf("%d", (sym >> b) & 1);
                t->stats.bits++;
            }
//          printf("\n");
        } else {
            for (int32_t b = 0; r == 0 && b < s->bits; b++) {
                r = bs->write_bit(bs, (s->path >> b) & 1);
                t->stats.bits++;
            }
        }
        ahc_inc_node_frequency(t, sym);
        #ifdef DEBUG
            ahc_verify_tree(t);
        #endif
        if (t->stats.each > 0 && i % t->stats.each == t->stats.each - 1) {
            const double d = (double)t->stats.each;
            size_t ds = t->stats.swaps - sc;
            size_t dm = t->stats.moves - mc;
            size_t du = t->stats.updates - uc;
            printf("[%07lld] swap: %5d \xF0\x9F\x94\x84 %.3f "
                            "move: %5d \xF0\x9F\x94\x83 %.3f "
                            "path: %5d \xF0\x9F\x8C\xB3 %.3f\n",
                    i,
                    ds, ds / d,
                    dm, dm / d,
                    du, du / d);
            sc = t->stats.swaps;
            mc = t->stats.moves;
            uc = t->stats.updates;
        }
    }
    while (r == 0 && bs->bits > 0) {
        r = bs->write_bit(bs, 0);
        t->stats.bits++;
    }
    return r;
}

static errno_t ahc_decode(ahc_tree_t* t, bitstream_t* bs, int32_t *symbol) {
    const int32_t m = t->n * 2 - 1;
    int32_t i = m - 1; // root
    bool bit = 0;
    errno_t r = bs->read_bit(bs, &bit);
    while (r == 0) {
        i = bit ? t->node[i].rix : t->node[i].lix;
        assert(0 <= i && i < m);
        if (t->node[i].lix < 0 && t->node[i].rix < 0) { break; } // leaf
        r = bs->read_bit(bs, &bit);
    }
    int32_t sym = i;
    if (sym == t->n / 2) { // NYE (Not Yet Encoded) escape:
//      printf("NYE: ");
        sym = 0;
        const int32_t bps = rt_const_log2_of_pow2(t->n) - 1;
        for (int32_t b = 0; b < bps && r == 0; b++) {
            r = bs->read_bit(bs, &bit);
            sym |= (bit << b);
//          printf("%d", bit);
        }
//      printf(" %02X\n", sym);
    }
    ahc_inc_node_frequency(t, sym);
    #ifdef DEBUG
        ahc_verify_tree(t);
    #endif
    *symbol = sym;
    return r;
}

static bool ahc_file_exist(const char* filename) {
    struct stat st = {0};
    return stat(filename, &st) == 0;
}

static errno_t ahc_file_size(FILE* f, size_t* size) {
    // on error returns (fpos_t)-1 and sets errno
    fpos_t pos = 0;
    if (fgetpos(f, &pos) != 0)      { return errno; }
    if (fseek(f, 0, SEEK_END) != 0) { return errno; }
    fpos_t eof = 0;
    if (fgetpos(f, &eof) != 0)      { return errno; }
    if (fseek(f, 0, SEEK_SET) != 0) { return errno; }
    if ((uint64_t)eof > SIZE_MAX)   { return E2BIG; }
    *size = (size_t)eof;
    return 0;
}

static errno_t ahc_read_fully(FILE* f, const uint8_t* *data, size_t *bytes) {
    size_t size = 0;
    errno_t r = ahc_file_size(f, &size);
    if (r != 0) { return r; }
    if (size > SIZE_MAX) { return E2BIG; }
    uint8_t* p = (uint8_t*)malloc(size); // does set errno on failure
    if (p == null) { return errno; }
    if (fread(p, 1, size, f) != (size_t)size) { free(p); return errno; }
    *data = p;
    *bytes = (size_t)size;
    return 0;
}

static errno_t ahc_read_whole_file(const char* fn,
                                   const uint8_t* *data, size_t *bytes) {
    FILE* f = null;
    errno_t r = fopen_s(&f, fn, "rb");
    if (r != 0) {
        printf("Failed to open file \"%s\": %s\n", fn, strerror(r));
        return r;
    }
    r = ahc_read_fully(f, data, bytes); // to the bh
    if (r != 0) {
        printf("Failed to read file \"%s\": %s\n", fn, strerror(r));
        fclose(f);
        return r;
    }
    // file was open for reading fclose() should not fail
    return fclose(f) == 0 ? 0 : errno;
}

static double ahc_shannon_entropy(ahc_tree_t* t) {
    enum { bps = 9, n = 1U << bps }; // bits per symbol
    assert(t->n == n);
    uint64_t freq[n];
    for (int32_t i = 0; i < n; i++) { freq[i] = t->node[i].freq; }
    return ahc_entropy(freq, n);
}

static errno_t ahc_test_encode_decode(const char* pathname,
                                      const uint8_t* data, size_t bytes) {
    enum { bps = 8 + 1, n = 1U << bps, m = n * 2 - 1 }; // `bps` bits per symbol
    bitstream_t bitstream;
    errno_t r = bitstream_create(&bitstream, bytes * 2);
    ahc_tree_t* t = null;
    if (r == 0) {
        t = (ahc_tree_t*)calloc(1, sizeof(ahc_tree_t));
        if (t == null) { r = ENOMEM; }
    }
    ahc_node_t* nodes = null;
    if (r == 0) {
        nodes = (ahc_node_t*)calloc(m, sizeof(ahc_node_t));
        if (nodes == null) { r = ENOMEM; }
    }
    if (r == 0) {
        ahc_init(t, nodes, m);
        // insert NYE (Not Yet Encoded) marker symbol
        ahc_inc_node_frequency(t, n / 2);
        if (ahc_stats) {
            t->stats.each = bytes / 16; // print stats 15 times
        }
        r = ahc_encode(t, (const uint8_t*)data, bytes, &bitstream);
    }
    if (r == 0) {
        assert(t->stats.bits % 64 == 0);
        int32_t depth = ahc_tree_depth(t);
        double percent = 100.0 * bitstream.bytes / (double)bytes;
        const char* fn = strrchr(pathname, '\\'); // shorten
        if (fn == null) { fn = strrchr(pathname, '/'); }
        if (fn == null) { fn = pathname; } else { fn++; }
        percent = 100.0 * t->stats.bits / ((double)bytes * 8.0);
        const double H = ahc_shannon_entropy(t);
        const double bps = (double)t->stats.bits / (double)bytes;
        printf("\xF0\x9F\x93\xA6 %7lld bytes into %7lld (%.1f%%) "
               "depth: %d/%d bps: %.2f H: %.2f \"%s\"\n",
                (uint64_t)bytes, (uint64_t)t->stats.bits / 8, percent,
                depth, t->depth, bps, H, fn);
    }
    free(nodes); nodes = null;
    free(t); t = null;
    if (r == 0) {
        t = (ahc_tree_t*)calloc(1, sizeof(ahc_tree_t));
        if (t == null) { r = ENOMEM; }
    }
    if (r == 0) {
        nodes = (ahc_node_t*)calloc(m, sizeof(ahc_node_t));
        if (nodes == null) { r = ENOMEM; }
    }
    if (r == 0) {
        ahc_init(t, nodes, m);
        // insert NYE (Not Yet Encoded) marker symbol
        ahc_inc_node_frequency(t, n / 2);
        uint8_t* decoded = (uint8_t*)malloc(bytes);
        if (decoded == null) {
            r = ENOMEM;
        } else {
            size_t decoded_bytes = 0;
            for (size_t i = 0; i < bytes; i++) {
                int32_t symbol = 0;
                r = ahc_decode(t, &bitstream, &symbol);
                decoded[i] = (uint8_t)symbol;
                decoded_bytes++;
            }
            if (memcmp(data, decoded, bytes) == 0) {
                printf("\xF0\x9F\x91\x8D %7lld bytes decoded to %7lld\n",
                       (uint64_t)bytes, (uint64_t)decoded_bytes);
            } else {
                size_t i = 0;
                while (i < bytes && data[i] == decoded[i]) { i++; }
                printf("Error @%lld of %lld\n", i, bytes);
                r = EINVAL;
            }
            free(decoded); decoded = null;
        }
    }
    free(nodes); nodes = null; // free(null) is OK
    free(t); t = null;
    if (bitstream.data != null) { bitstream_dispose(&bitstream); }
    return r;
}

static errno_t ahc_test_file(const char* filename) {
    errno_t r = ENOENT;
    if (ahc_file_exist(filename)) {
        const uint8_t* data = null;
        size_t bytes = 0;
        r = ahc_read_whole_file(filename, &data, &bytes);
        assert(r == 0);
        if (r == 0) {
            r = ahc_test_encode_decode(filename, data, bytes);
            free((void*)data);
        }
    }
    return r;
}

static errno_t locate_test_folder(void) {
    // on Unix systems with "make" executable usually resided
    // and is run from root of repository... On Windows with
    // MSVC it is buried inside bin/... folder depths
    // on X Code in MacOS it can be completely out of tree.
    // So we need to find the test files.
    for (;;) {
        if (ahc_file_exist("test/bible.txt")) { return 0; }
        if (chdir("..") != 0) { return errno; }
    }
}

int main(int argc, const char* argv[]) {
    printf("Hello\xF0\x9F\x91\x8B "
           "Adaptive Huffman Coding\xF0\x9F\x8C\x8D!\n");
    bool t = false; // --test
    bool s = false; // --stat
    bool c = false; // --check (verify tree on every byte)
    for (int i = 1; i < argc; i++) {
        t |= strcmp(argv[i], "--test") == 0;
        s |= strcmp(argv[i], "--stat") == 0;
        c |= strcmp(argv[i], "--check") == 0;
    }
    ahc_stats = s;
    ahc_verify = s;
    if (t) { test(); }
    errno_t r = locate_test_folder();
    if (r == 0) { r = ahc_test_file("test/hhgttg.txt"); }
    if (r == 0) { r = ahc_test_file("test/bible.txt"); }
    if (r == 0) { r = ahc_test_file("test/arm64.elf"); }
    if (r == 0) { r = ahc_test_file("test/snow-crash.bmp"); }
    if (r == 0) { r = ahc_test_file("test/sqlite3.c.txt"); }
    if (r == 0) { r = ahc_test_file(__FILE__); }
    if (r == 0) { r = ahc_test_file(argv[0]); }
    if (r != 0) {
        printf("Goodbye \xF0\x9F\x98\x88 cruel \xF0\x9F\x98\xB1 "
               "Universe \xF0\x9F\x8C\xA0\xF0\x9F\x8C\x8C..."
               "\xF0\x9F\x92\xA4\n");
        printf("%s\n", strerror(r));
    }
    return r;
}
