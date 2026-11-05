// The haplotype assembler (src/hap_asm.c) on hand-built fragment sets with a
// known answer. Build and run with `make test`.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hap_asm.h"

// Globals main.c defines for the program; the test links everything else.
const char PROG[20] = "test_hap_asm";
int LONGCALLD_VERBOSE = 0;
char *CMD = "test_hap_asm";

static int n_pass = 0, n_fail = 0;
static void chk(const char *name, int ok) {
    printf("  %s  %s\n", ok ? "PASS" : "FAIL", name);
    if (ok) ++n_pass; else ++n_fail;
}

// A fragment set under construction.
typedef struct { ha_t h; int m_frag, m_call; } bld_t;

static void bld_init(bld_t *b, int n_var) {
    memset(b, 0, sizeof(*b));
    b->h.n_var = n_var;
    b->h.var = (ha_var_t*)calloc(n_var, sizeof(ha_var_t));
    for (int v = 0; v < n_var; ++v) {
        b->h.var[v].pos = 1000 + 100 * v; b->h.var[v].phasable = 1; b->h.var[v].geno = 2;
        b->h.var[v].g[0] = 0; b->h.var[v].g[1] = 1; // seed every site het
    }
}
static void bld_frag(bld_t *b) {
    if (b->h.n_frag == b->m_frag) { b->m_frag = b->m_frag ? 2 * b->m_frag : 64; b->h.frag = (ha_frag_t*)realloc(b->h.frag, b->m_frag * sizeof(ha_frag_t)); }
    b->h.frag[b->h.n_frag].off = b->h.n_call; b->h.frag[b->h.n_frag].n = 0; b->h.frag[b->h.n_frag].smp = 0;
    ++b->h.n_frag;
}
static void bld_call(bld_t *b, int v, int a, double e, int tag_only) {
    if (b->h.n_call == b->m_call) { b->m_call = b->m_call ? 2 * b->m_call : 256; b->h.call = (ha_call_t*)realloc(b->h.call, b->m_call * sizeof(ha_call_t)); }
    ha_call_t *c = b->h.call + b->h.n_call++;
    c->v = v; c->a = (uint8_t)a; c->lm = (float)log(1 - e); c->lx = (float)log(e);
    c->ok = !tag_only && e < pow(10, -0.7); c->tag_only = (uint8_t)tag_only;
    b->h.frag[b->h.n_frag - 1].n++;
}

// Deterministic noise.
static unsigned rng = 12345;
static double urand(void) { rng = rng * 1103515245u + 12345u; return (rng >> 8) / 16777216.0; }

int main(void) {
    ha_opt_t o; ha_opt_init(&o);

    // 1. A known haplotype over 40 hets, 30 reads per haplotype each spanning
    //    8 consecutive sites, 5% read error, seeded with every orientation
    //    scrambled. The assembly has to recover it up to a global flip.
    {
        int nv = 40; int truth[40];
        for (int v = 0; v < nv; ++v) truth[v] = urand() < 0.5;
        bld_t b; bld_init(&b, nv);
        for (int v = 0; v < nv; ++v) if (urand() < 0.5) { b.h.var[v].g[0] = 1; b.h.var[v].g[1] = 0; }
        for (int r = 0; r < 330; ++r) { // even depth: every start 10 times, both haplotypes
            int s = r % (nv - 7), hap = (r / (nv - 7)) & 1;
            bld_frag(&b);
            for (int v = s; v < s + 8; ++v) {
                int a = hap ? 1 - truth[v] : truth[v];
                if (urand() < 0.05) a = 1 - a;
                bld_call(&b, v, a, 0.05, 0);
            }
        }
        ha_run(&o, &b.h);
        int same = 0, flip = 0, het = 0;
        for (int v = 0; v < nv; ++v) {
            if (b.h.var[v].g[0] == b.h.var[v].g[1]) continue;
            ++het; if (b.h.var[v].g[0] == truth[v]) ++same; else ++flip;
        }
        chk("truth haplotype recovered from scrambled seeds", het == nv && (same == nv || flip == nv));
        if (het != nv || (same != nv && flip != nv)) {
            for (int v = 0; v < nv; ++v)
                printf("    v%d truth=%d g=%d%d ps=%ld pq=%.1f sq=%.1f incons=%d\n", v, truth[v], b.h.var[v].g[0], b.h.var[v].g[1],
                       (long)b.h.var[v].ps, b.h.var[v].pq, b.h.var[v].sq, b.h.var[v].inconsistent);
        }
        int one_block = 1;
        for (int v = 1; v < nv; ++v) if (b.h.var[v].ps != b.h.var[0].ps) one_block = 0;
        chk("read-linked hets share one block named after the first", one_block && b.h.var[0].ps == 1000);
        int tagged = 0;
        for (int f = 0; f < b.h.n_frag; ++f) tagged += b.h.frag_hap[f] != 0;
        chk("every read is tagged", tagged == b.h.n_frag);
        ha_free(&b.h);
    }

    // 2. Site 2 has its alt on both haplotypes (an artifact), sites 0,1,3,4 are
    //    a clean het chain. Site 2 must not be phased; the rest stay one block.
    {
        bld_t b; bld_init(&b, 5);
        for (int r = 0; r < 40; ++r) {
            int hap = r & 1;
            bld_frag(&b);
            for (int v = 0; v < 5; ++v) {
                int a = v == 2 ? (r % 4 < 2) : hap; // site 2: alt on reads 0,1 of every 4, whatever their haplotype
                bld_call(&b, v, a, 0.02, 0);
            }
        }
        ha_run(&o, &b.h);
        chk("alt spread over both haplotypes -> HapInconsistent", b.h.var[2].inconsistent && b.h.var[2].ps == 0);
        chk("its neighbours stay one block", b.h.var[0].ps == b.h.var[4].ps && b.h.var[0].ps > 0 &&
                                             !b.h.var[0].inconsistent && !b.h.var[4].inconsistent);
        ha_free(&b.h);
    }

    // 3. Two hets no read links to each other: each is a block of its own.
    {
        bld_t b; bld_init(&b, 2);
        for (int r = 0; r < 20; ++r) { bld_frag(&b); bld_call(&b, r < 10 ? 0 : 1, r & 1, 0.02, 0); }
        ha_run(&o, &b.h);
        chk("an isolated het is its own block", b.h.var[0].ps == 1000 && b.h.var[1].ps == 1100);
        ha_free(&b.h);
    }

    // 4. A site every read shows alt at, seeded het, is genotyped 1/1.
    {
        bld_t b; bld_init(&b, 3);
        for (int r = 0; r < 30; ++r) {
            bld_frag(&b);
            bld_call(&b, 0, r & 1, 0.02, 0); bld_call(&b, 1, 1, 0.02, 0); bld_call(&b, 2, r & 1, 0.02, 0);
        }
        ha_run(&o, &b.h);
        chk("hom-alt site seeded het -> 1/1", b.h.var[1].g[0] == 1 && b.h.var[1].g[1] == 1);
        chk("and the hets either side are still one block", b.h.var[0].ps > 0 && b.h.var[0].ps == b.h.var[2].ps);
        ha_free(&b.h);
    }

    // 5. Tagging: two agreeing calls at 10% error give 0.988 and are tagged at
    //    the 0.9 default; one call at 15% error (0.85) is not.
    {
        bld_t b; bld_init(&b, 2);
        for (int r = 0; r < 20; ++r) { bld_frag(&b); bld_call(&b, 0, r & 1, 0.02, 0); bld_call(&b, 1, r & 1, 0.02, 0); }
        bld_frag(&b); bld_call(&b, 0, 1, 0.1, 0); bld_call(&b, 1, 1, 0.1, 0);
        bld_frag(&b); bld_call(&b, 0, 1, 0.15, 0);
        ha_run(&o, &b.h);
        chk("two agreeing calls at 10% error are tagged", b.h.frag_hap[20] != 0);
        chk("one call at 15% error is not", b.h.frag_hap[21] == 0);
        ha_free(&b.h);
    }

    // 6. Tag-only calls tag a read but never link: two blocks stay two even
    //    though a skipped read covers both.
    {
        bld_t b; bld_init(&b, 4);
        for (int r = 0; r < 20; ++r) { bld_frag(&b); bld_call(&b, 0, r & 1, 0.02, 0); bld_call(&b, 1, r & 1, 0.02, 0); }
        for (int r = 0; r < 20; ++r) { bld_frag(&b); bld_call(&b, 2, r & 1, 0.02, 0); bld_call(&b, 3, r & 1, 0.02, 0); }
        bld_frag(&b); bld_call(&b, 0, 1, 0.02, 1); bld_call(&b, 1, 1, 0.02, 1); bld_call(&b, 2, 1, 0.02, 1);
        ha_run(&o, &b.h);
        chk("tag-only calls do not link blocks", b.h.var[0].ps != b.h.var[2].ps);
        chk("but they do tag the read", b.h.frag_hap[40] != 0);
        ha_free(&b.h);
    }

    // 6. Five hets read by 20 reads per haplotype at 2% error, except that at
    //    site 2 a third of haplotype 1's reads show haplotype 2's allele. With
    //    the concordance test on, site 2 is HapDiscordant and leaves the
    //    phasing; the others still link (PHASING). With it off, nothing is.
    for (int on = 1; on >= 0; --on) {
        ha_opt_t o2; ha_opt_init(&o2);
        o2.min_hap_conc = on ? 0.9 : 0; o2.conc_err[0] = o2.conc_err[1] = 0.05;
        int nv = 5; bld_t b; bld_init(&b, nv);
        for (int r = 0; r < 40; ++r) {
            int hap = r & 1; // 0: haplotype carrying 0,1,0,1,0; 1: the other
            bld_frag(&b);
            for (int v = 0; v < nv; ++v) {
                int a = (v & 1) ^ hap;
                if (urand() < 0.02) a ^= 1;
                if (v == 2 && hap == 0 && r % 6 == 0) a ^= 1;
                bld_call(&b, v, a, 0.02, 0);
            }
        }
        ha_run(&o2, &b.h);
        int others_linked = 1;
        for (int v = 0; v < nv; ++v) if (v != 2 && !b.h.var[v].linked) others_linked = 0;
        if (on) {
            chk("a third of a haplotype's reads on the other allele -> HapDiscordant", b.h.var[2].discordant);
            chk("and it is not used for phasing", !b.h.var[2].linked);
            chk("its neighbours are (PHASING)", others_linked);
            chk("and are not discordant", !b.h.var[0].discordant && !b.h.var[4].discordant);
        } else {
            chk("concordance off: nothing is discordant", !b.h.var[2].discordant);
        }
        ha_free(&b.h);
    }

    printf("%d passed, %d failed\n", n_pass, n_fail);
    return n_fail ? 1 : 0;
}
