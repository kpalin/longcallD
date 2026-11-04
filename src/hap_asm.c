#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "hap_asm.h"
#include "utils.h"
#include "bam_utils.h"
#include "collect_var.h"
#include "call_var_main.h"
#include "var_qual.h"

#define HA_LN10 2.302585092994046
#define HA_EPS 1e-9
#define HA_MAX_SEARCH 200   // local-search iterations per round
#define HA_MAX_GENO_PASS 10 // genotyping passes per round
#define HA_MAX_PRUNE 3      // prune / re-assemble cycles

void ha_opt_init(ha_opt_t *o) {
    memset(o, 0, sizeof(*o));
    o->max_miscall = pow(10.0, -7.0 / 10.0); // Longshot --min_allele_qual 7
    o->min_pq = 10.0;
    o->hap_assign_p = 0.9;
    // Longshot's per-site priors: het/hom SNV 1e-3 / 5e-4, het/hom indel 1e-5 / 5e-5.
    o->prior[0][1] = 1e-3; o->prior[0][2] = 5e-4;
    o->prior[1][1] = 1e-5; o->prior[1][2] = 5e-5;
    for (int t = 0; t < 2; ++t) o->prior[t][0] = 1.0 - o->prior[t][1] - o->prior[t][2];
    o->geno_smp = -1;
    o->max_rounds = 10;
}

void ha_free(ha_t *h) {
    free(h->var); free(h->frag); free(h->call); free(h->frag_ps); free(h->frag_hap);
    memset(h, 0, sizeof(*h));
}

static inline double lse(double a, double b) { // ln(e^a + e^b)
    if (a < b) { double t = a; a = b; b = t; }
    if (b == -INFINITY) return a;
    return a + log1p(exp(b - a));
}

static inline double cterm(const ha_call_t *c, int allele) { return allele == c->a ? c->lm : c->lx; }

typedef struct {
    const ha_opt_t *o; ha_t *h;
    int *vc_off, *vc_idx;   // variant -> its calls (indices into h->call)
    int *call_frag;         // call -> fragment
    uint8_t *inP;           // in the phasing graph this round: phasable, unpruned, het
    int *uf, *comp;         // union-find; variant -> block root (its first variant), -1 if not in P
    int *fblk;              // fragment -> block root, -1
    double *fp0, *fp1;      // per fragment: ln P(f|H1), ln P(f|H2) over its block's trusted calls
    int *rank, *sw_off;     // variant -> rank in its block; block root -> offset of its switch array
    int *bsz;               // block root -> number of variants
    double *sw;             // switch gains, per block, one per cut
} ws_t;

static int uf_find(int *uf, int x) {
    while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; }
    return x;
}

static void ws_init(ws_t *w, const ha_opt_t *o, ha_t *h) {
    memset(w, 0, sizeof(*w));
    w->o = o; w->h = h;
    int nv = h->n_var, nf = h->n_frag, nc = h->n_call;
    w->vc_off = (int*)_err_calloc(nv + 1, sizeof(int));
    w->vc_idx = (int*)_err_malloc((nc > 0 ? nc : 1) * sizeof(int));
    w->call_frag = (int*)_err_malloc((nc > 0 ? nc : 1) * sizeof(int));
    for (int f = 0; f < nf; ++f)
        for (int k = h->frag[f].off; k < h->frag[f].off + h->frag[f].n; ++k) w->call_frag[k] = f;
    for (int k = 0; k < nc; ++k) w->vc_off[h->call[k].v + 1]++;
    for (int v = 0; v < nv; ++v) w->vc_off[v + 1] += w->vc_off[v];
    int *fill = (int*)_err_calloc(nv > 0 ? nv : 1, sizeof(int));
    for (int k = 0; k < nc; ++k) { int v = h->call[k].v; w->vc_idx[w->vc_off[v] + fill[v]++] = k; }
    free(fill);
    w->inP = (uint8_t*)_err_calloc(nv > 0 ? nv : 1, 1);
    w->uf = (int*)_err_malloc((nv > 0 ? nv : 1) * sizeof(int));
    w->comp = (int*)_err_malloc((nv > 0 ? nv : 1) * sizeof(int));
    w->rank = (int*)_err_malloc((nv > 0 ? nv : 1) * sizeof(int));
    w->sw_off = (int*)_err_malloc((nv > 0 ? nv : 1) * sizeof(int));
    w->bsz = (int*)_err_calloc(nv > 0 ? nv : 1, sizeof(int));
    w->sw = (double*)_err_malloc((nv > 0 ? nv : 1) * sizeof(double));
    w->fblk = (int*)_err_malloc((nf > 0 ? nf : 1) * sizeof(int));
    w->fp0 = (double*)_err_calloc(nf > 0 ? nf : 1, sizeof(double));
    w->fp1 = (double*)_err_calloc(nf > 0 ? nf : 1, sizeof(double));
}

static void ws_free(ws_t *w) {
    free(w->vc_off); free(w->vc_idx); free(w->call_frag); free(w->inP); free(w->uf); free(w->comp);
    free(w->rank); free(w->sw_off); free(w->bsz); free(w->sw); free(w->fblk); free(w->fp0); free(w->fp1);
}

static inline int is_het(const ha_var_t *v) { return v->g[0] != v->g[1]; }

// The phasing graph: trusted calls of one read link the hets they fall on,
// consecutive pairs being enough since linkage is transitive. Blocks are its
// connected components, named after their first variant.
static void build_blocks(ws_t *w) {
    ha_t *h = w->h;
    for (int v = 0; v < h->n_var; ++v) {
        w->inP[v] = h->var[v].phasable && !h->var[v].pruned && is_het(h->var + v);
        w->uf[v] = v;
    }
    for (int f = 0; f < h->n_frag; ++f) {
        int last = -1;
        for (int k = h->frag[f].off; k < h->frag[f].off + h->frag[f].n; ++k) {
            const ha_call_t *c = h->call + k;
            if (!c->ok || !w->inP[c->v]) continue;
            if (last >= 0) {
                int a = uf_find(w->uf, last), b = uf_find(w->uf, c->v);
                if (a != b) { if (a < b) w->uf[b] = a; else w->uf[a] = b; }
            }
            last = c->v;
        }
    }
    // Union by smaller index keeps every root at its block's first variant.
    for (int v = 0; v < h->n_var; ++v) w->comp[v] = w->inP[v] ? uf_find(w->uf, v) : -1;
    for (int f = 0; f < h->n_frag; ++f) {
        w->fblk[f] = -1;
        for (int k = h->frag[f].off; k < h->frag[f].off + h->frag[f].n; ++k) {
            const ha_call_t *c = h->call + k;
            if (c->ok && w->inP[c->v]) { w->fblk[f] = w->comp[c->v]; break; }
        }
    }
    // Rank of each variant within its block, and where the block's cuts live.
    int off = 0;
    memset(w->bsz, 0, (h->n_var > 0 ? h->n_var : 1) * sizeof(int));
    for (int v = 0; v < h->n_var; ++v) if (w->comp[v] >= 0) w->rank[v] = w->bsz[w->comp[v]]++;
    for (int v = 0; v < h->n_var; ++v) if (w->comp[v] == v) { w->sw_off[v] = off; off += w->bsz[v]; }
}

static void frag_sums(ws_t *w) {
    ha_t *h = w->h;
    for (int f = 0; f < h->n_frag; ++f) {
        double s0 = 0, s1 = 0;
        if (w->fblk[f] >= 0)
            for (int k = h->frag[f].off; k < h->frag[f].off + h->frag[f].n; ++k) {
                const ha_call_t *c = h->call + k;
                if (!c->ok || !w->inP[c->v]) continue;
                s0 += cterm(c, h->var[c->v].g[0]);
                s1 += cterm(c, h->var[c->v].g[1]);
            }
        w->fp0[f] = s0; w->fp1[f] = s1;
    }
}

// Log-likelihood gain of swapping variant v's two haplotype alleles.
static double flip_gain(ws_t *w, int v) {
    ha_t *h = w->h; const ha_var_t *x = h->var + v;
    double d = 0;
    for (int i = w->vc_off[v]; i < w->vc_off[v + 1]; ++i) {
        int k = w->vc_idx[i]; const ha_call_t *c = h->call + k;
        if (!c->ok) continue;
        int f = w->call_frag[k];
        double t1 = cterm(c, x->g[0]), t2 = cterm(c, x->g[1]);
        double a0 = w->fp0[f], a1 = w->fp1[f];
        d += lse(a0 - t1 + t2, a1 - t2 + t1) - lse(a0, a1);
    }
    return d;
}

static void apply_flip(ws_t *w, int v) {
    ha_t *h = w->h; ha_var_t *x = h->var + v;
    for (int i = w->vc_off[v]; i < w->vc_off[v + 1]; ++i) {
        int k = w->vc_idx[i]; const ha_call_t *c = h->call + k;
        if (!c->ok) continue;
        int f = w->call_frag[k];
        double t1 = cterm(c, x->g[0]), t2 = cterm(c, x->g[1]);
        w->fp0[f] += t2 - t1; w->fp1[f] += t1 - t2;
    }
    int8_t t = x->g[0]; x->g[0] = x->g[1]; x->g[1] = t;
}

// Gain of flipping every variant of a block after each cut. A fragment changes
// only if it has calls on both sides of the cut, so the whole pass is linear
// in the calls: each fragment adds one value over each run of cuts it spans.
static void switch_gains(ws_t *w) {
    ha_t *h = w->h;
    int n = 0;
    for (int v = 0; v < h->n_var; ++v) if (w->comp[v] >= 0) ++n;
    double *diff = (double*)_err_calloc(n + 1, sizeof(double));
    double *L0 = NULL, *L1 = NULL; int *r = NULL, m_tmp = 0;
    for (int f = 0; f < h->n_frag; ++f) {
        int b = w->fblk[f];
        if (b < 0) continue;
        int k0 = h->frag[f].off, nk = 0;
        if (h->frag[f].n + 1 > m_tmp) {
            m_tmp = h->frag[f].n + 1;
            L0 = (double*)_err_realloc(L0, m_tmp * sizeof(double));
            L1 = (double*)_err_realloc(L1, m_tmp * sizeof(double));
            r = (int*)_err_realloc(r, m_tmp * sizeof(int));
        }
        double s0 = 0, s1 = 0;
        for (int k = k0; k < k0 + h->frag[f].n; ++k) {
            const ha_call_t *c = h->call + k;
            if (!c->ok || !w->inP[c->v]) continue;
            s0 += cterm(c, h->var[c->v].g[0]); s1 += cterm(c, h->var[c->v].g[1]);
            L0[nk] = s0; L1[nk] = s1; r[nk] = w->rank[c->v]; ++nk;
        }
        if (nk < 2) continue;
        double base = lse(s0, s1);
        for (int i = 0; i < nk - 1; ++i) {
            double n0 = L0[i] + (s1 - L1[i]), n1 = L1[i] + (s0 - L0[i]);
            double d = lse(n0, n1) - base;
            diff[w->sw_off[b] + r[i]] += d;
            diff[w->sw_off[b] + r[i + 1]] -= d;
        }
    }
    // Prefix sums per block: sw[off + i] is the gain of the cut after rank i.
    for (int v = 0; v < h->n_var; ++v) {
        if (w->comp[v] != v) continue;
        double acc = 0;
        for (int i = w->sw_off[v]; i < w->sw_off[v] + w->bsz[v]; ++i) { acc += diff[i]; w->sw[i] = acc; }
    }
    free(diff); free(L0); free(L1); free(r);
}

// One local-search round: single-variant flips, then the best switch of each
// block, until neither improves the fragment likelihood.
static void local_search(ws_t *w) {
    ha_t *h = w->h;
    const int *bsize = w->bsz;
    int *best_cut = (int*)_err_malloc((h->n_var > 0 ? h->n_var : 1) * sizeof(int));
    for (int it = 0; it < HA_MAX_SEARCH; ++it) {
        int changed = 0;
        for (int v = 0; v < h->n_var; ++v)
            if (w->comp[v] >= 0 && flip_gain(w, v) > HA_EPS) { apply_flip(w, v); changed = 1; }
        switch_gains(w);
        for (int v = 0; v < h->n_var; ++v) {
            best_cut[v] = -1;
            if (w->comp[v] != v || bsize[v] < 2) continue;
            double best = 1e-6;
            for (int i = 0; i < bsize[v] - 1; ++i)
                if (w->sw[w->sw_off[v] + i] > best) { best = w->sw[w->sw_off[v] + i]; best_cut[v] = i; }
        }
        for (int v = 0; v < h->n_var; ++v) {
            int b = w->comp[v];
            if (b < 0 || best_cut[b] < 0 || w->rank[v] <= best_cut[b]) continue;
            int8_t t = h->var[v].g[0]; h->var[v].g[0] = h->var[v].g[1]; h->var[v].g[1] = t;
            changed = 1;
        }
        frag_sums(w);
        if (!changed) break;
    }
    free(best_cut);
}

static double lprior(const ha_opt_t *o, const ha_var_t *x, int a0, int a1) {
    const double *p = o->prior[x->is_indel ? 1 : 0];
    if (a0 != a1) return log(p[1] * 0.5);
    return log(a0 ? p[2] : p[0]);
}

// The block a variant is genotyped against: its own if it is in the graph,
// otherwise the one most of its reads belong to.
static int ref_block(ws_t *w, int v) {
    if (w->comp[v] >= 0) return w->comp[v];
    int best = -1, nb = 0, cur = -1, nc = 0;
    // reads of one block are not contiguous here, so a Boyer-Moore majority
    for (int i = w->vc_off[v]; i < w->vc_off[v + 1]; ++i) {
        int b = w->fblk[w->call_frag[w->vc_idx[i]]];
        if (b < 0) continue;
        if (nc == 0) { cur = b; nc = 1; } else if (b == cur) ++nc; else --nc;
    }
    if (cur < 0) return -1;
    for (int i = w->vc_off[v]; i < w->vc_off[v + 1]; ++i)
        if (w->fblk[w->call_frag[w->vc_idx[i]]] == cur) ++nb;
    best = nb > 0 ? cur : -1;
    return best;
}

// P(reads at v | ordered genotype) on the reads' haplotype likelihoods, the
// variant's own term swapped for the candidate's. Returns the four ordered
// log-posteriors (00, 01, 10, 11), unnormalised.
static void geno_ll(ws_t *w, int v, double ll[4]) {
    ha_t *h = w->h; const ha_opt_t *o = w->o; const ha_var_t *x = h->var + v;
    static const int G[4][2] = {{0,0},{0,1},{1,0},{1,1}};
    int B = ref_block(w, v);
    for (int g = 0; g < 4; ++g) ll[g] = lprior(o, x, G[g][0], G[g][1]);
    for (int i = w->vc_off[v]; i < w->vc_off[v + 1]; ++i) {
        int k = w->vc_idx[i]; const ha_call_t *c = h->call + k;
        int f = w->call_frag[k];
        if (c->tag_only || (o->geno_smp >= 0 && h->frag[f].smp != o->geno_smp)) continue;
        if (B >= 0 && w->fblk[f] == B) {
            double b0 = w->fp0[f], b1 = w->fp1[f];
            if (w->inP[v] && c->ok) { b0 -= cterm(c, x->g[0]); b1 -= cterm(c, x->g[1]); }
            double base = lse(b0, b1);
            for (int g = 0; g < 4; ++g)
                ll[g] += lse(b0 + cterm(c, G[g][0]), b1 + cterm(c, G[g][1])) - base;
        } else {
            for (int g = 0; g < 4; ++g)
                ll[g] += lse(cterm(c, G[g][0]), cterm(c, G[g][1])) - log(2.0);
        }
    }
}

// One greedy pass of haplotype-informed genotyping. Returns how many changed.
static int genotype_pass(ws_t *w) {
    ha_t *h = w->h;
    static const int G[4][2] = {{0,0},{0,1},{1,0},{1,1}};
    int n_changed = 0;
    for (int v = 0; v < h->n_var; ++v) {
        ha_var_t *x = h->var + v;
        if (x->geno == 0) continue;
        double ll[4];
        geno_ll(w, v, ll);
        int allow[4] = {1, 1, 1, 1};
        if (x->geno == 1) { // zygosity fixed: only the orientation is decided here
            int het = x->ug == 1;
            allow[0] = !het && x->ug == 0; allow[3] = !het && x->ug == 2; allow[1] = allow[2] = het;
        }
        int best = -1;
        for (int g = 0; g < 4; ++g) if (allow[g] && (best < 0 || ll[g] > ll[best])) best = g;
        double tot = -INFINITY;
        for (int g = 0; g < 4; ++g) tot = lse(tot, ll[g]);
        double p_best = exp(ll[best] - tot);
        if (best == 1 || best == 2) p_best += exp(ll[3 - best] - tot); // phase is not part of the genotype
        double e = 1.0 - p_best; if (e < 1e-30) e = 1e-30;
        x->gq = (float)(-10.0 * log10(e));
        if (G[best][0] == x->g[0] && G[best][1] == x->g[1]) continue;
        // Keep the reads' haplotype sums in step with a variant already in them.
        if (w->inP[v]) {
            for (int i = w->vc_off[v]; i < w->vc_off[v + 1]; ++i) {
                int k = w->vc_idx[i]; const ha_call_t *c = h->call + k;
                if (!c->ok) continue;
                int f = w->call_frag[k];
                w->fp0[f] += cterm(c, G[best][0]) - cterm(c, x->g[0]);
                w->fp1[f] += cterm(c, G[best][1]) - cterm(c, x->g[1]);
            }
        }
        x->g[0] = (int8_t)G[best][0]; x->g[1] = (int8_t)G[best][1];
        ++n_changed;
    }
    return n_changed;
}

static double total_ll(ws_t *w) {
    ha_t *h = w->h; double s = 0;
    for (int v = 0; v < h->n_var; ++v)
        if (h->var[v].geno) s += lprior(w->o, h->var + v, h->var[v].g[0], h->var[v].g[1]);
    for (int f = 0; f < h->n_frag; ++f) if (w->fblk[f] >= 0) s += lse(w->fp0[f], w->fp1[f]);
    return s;
}

// Number of reads that tie v to another het of its block.
static int n_linking(ws_t *w, int v) {
    ha_t *h = w->h; int n = 0;
    for (int i = w->vc_off[v]; i < w->vc_off[v + 1]; ++i) {
        int k = w->vc_idx[i];
        if (!h->call[k].ok) continue;
        int f = w->call_frag[k], m = 0;
        for (int j = h->frag[f].off; j < h->frag[f].off + h->frag[f].n && m < 2; ++j)
            if (h->call[j].ok && w->inP[h->call[j].v]) ++m;
        if (m >= 2) ++n;
    }
    return n;
}

static void assemble(ws_t *w) {
    double prev = -INFINITY;
    for (int round = 0; round < w->o->max_rounds; ++round) {
        build_blocks(w);
        frag_sums(w);
        local_search(w);
        int n_changed = 0;
        for (int p = 0; p < HA_MAX_GENO_PASS; ++p) {
            int c = genotype_pass(w);
            n_changed += c;
            if (c == 0) break;
        }
        double cur = total_ll(w);
        if (n_changed == 0) break;
        if (prev > -INFINITY && fabs(cur - prev) <= 1e-4 * fabs(prev)) break;
        prev = cur;
    }
    build_blocks(w);
    frag_sums(w);
}

void ha_run(const ha_opt_t *o, ha_t *h) {
    ws_t w; ws_init(&w, o, h);
    // Genotype before phasing: every read on equal terms (Longshot's
    // call_genotypes_no_haplotypes). A seeded het keeps its orientation.
    for (int v = 0; v < h->n_var; ++v) {
        ha_var_t *x = h->var + v;
        x->pq = x->sq = x->gq = -1; x->ps = 0; x->pruned = x->inconsistent = 0;
        if (x->geno == 0) { x->ug = (int8_t)(x->g[0] + x->g[1]); continue; }
        if (x->geno == 1) { x->ug = (int8_t)(x->g[0] + x->g[1]); continue; }
        double ll[3];
        for (int z = 0; z < 3; ++z) ll[z] = z == 1 ? log(o->prior[x->is_indel][1]) : log(o->prior[x->is_indel][z]);
        for (int i = w.vc_off[v]; i < w.vc_off[v + 1]; ++i) {
            const ha_call_t *c = h->call + w.vc_idx[i];
            if (c->tag_only || (o->geno_smp >= 0 && h->frag[w.call_frag[w.vc_idx[i]]].smp != o->geno_smp)) continue;
            ll[0] += cterm(c, 0); ll[2] += cterm(c, 1);
            ll[1] += lse(cterm(c, 0), cterm(c, 1)) - log(2.0);
        }
        int z = 0;
        for (int k = 1; k < 3; ++k) if (ll[k] > ll[z]) z = k;
        x->ug = (int8_t)z;
        if (z == 1) { if (!is_het(x)) { x->g[0] = 0; x->g[1] = 1; } }
        else x->g[0] = x->g[1] = (int8_t)(z == 2);
    }
    assemble(&w);
    // Prune the hets the reads cannot agree on, then assemble again without them.
    for (int cyc = 0; cyc < HA_MAX_PRUNE; ++cyc) {
        int n_pruned = 0;
        for (int v = 0; v < h->n_var; ++v) {
            if (w.comp[v] < 0) continue;
            double pq = -flip_gain(&w, v) * 10.0 / HA_LN10;
            h->var[v].pq = (float)pq;
            if (pq < o->min_pq && n_linking(&w, v) >= 2) { h->var[v].pruned = 1; ++n_pruned; }
        }
        if (n_pruned == 0) break;
        assemble(&w);
    }
    // A het in the graph with no trusted link of its own, whose reads are
    // nevertheless on another block, belongs to that block: genotype it there
    // rather than leave it a block of one inside another's span.
    {
        uint8_t *solo = (uint8_t*)_err_calloc(h->n_var > 0 ? h->n_var : 1, 1);
        int n_solo = 0;
        for (int v = 0; v < h->n_var; ++v) {
            if (w.comp[v] != v || w.bsz[v] != 1) continue;
            for (int i = w.vc_off[v]; i < w.vc_off[v + 1]; ++i) {
                int b = w.fblk[w.call_frag[w.vc_idx[i]]];
                if (b >= 0 && b != v) { solo[v] = 1; ++n_solo; break; }
            }
        }
        if (n_solo > 0) {
            for (int v = 0; v < h->n_var; ++v) if (solo[v]) h->var[v].phasable = 0;
            build_blocks(&w); frag_sums(&w);
            for (int p = 0; p < HA_MAX_GENO_PASS && genotype_pass(&w) > 0; ++p) ;
            build_blocks(&w); frag_sums(&w);
            for (int v = 0; v < h->n_var; ++v) if (solo[v]) h->var[v].phasable = 1;
        }
        free(solo);
    }
    // Phase quality and switch quality on the final haplotypes.
    switch_gains(&w);
    for (int v = 0; v < h->n_var; ++v) {
        ha_var_t *x = h->var + v;
        if (w.comp[v] >= 0) {
            x->pq = (float)(-flip_gain(&w, v) * 10.0 / HA_LN10);
            if (w.rank[v] > 0) x->sq = (float)(-w.sw[w.sw_off[w.comp[v]] + w.rank[v] - 1] * 10.0 / HA_LN10);
            x->ps = h->var[w.comp[v]].pos;
        } else if (is_het(x) && !x->pruned) {
            // A het outside the graph is phased against the block its reads are
            // on; with none, it is a block of its own -- unless the phasing must
            // be anchored on particular sites (paired mode: normal hets) and
            // this is not one of them.
            int b = ref_block(&w, v);
            x->ps = b >= 0 ? h->var[b].pos : (o->solo_needs_phasable && !x->phasable ? 0 : x->pos);
        }
        // HapInconsistent: pruned, or called before phasing but carried by no
        // haplotype once the reads are phased -- its alt reads are spread over
        // both. Reported with the unphased genotype and no phase set.
        if (x->pruned || (x->geno == 2 && x->ug > 0 && x->g[0] == 0 && x->g[1] == 0)) {
            x->inconsistent = 1; x->ps = 0;
            if (x->ug == 2) x->g[0] = x->g[1] = 1; else { x->g[0] = 0; x->g[1] = 1; }
        }
    }
    // Read haplotypes: the posterior over the read's calls at the phased hets
    // of the block it has most calls in. Every call counts here, trusted or
    // not -- the trust floor decides what may link two variants, not whether
    // a read that plainly sits on one haplotype is told so.
    h->frag_ps = (hts_pos_t*)_err_calloc(h->n_frag > 0 ? h->n_frag : 1, sizeof(hts_pos_t));
    h->frag_hap = (int8_t*)_err_calloc(h->n_frag > 0 ? h->n_frag : 1, sizeof(int8_t));
    for (int f = 0; f < h->n_frag; ++f) {
        hts_pos_t cur = 0; int nc = 0; // Boyer-Moore majority over the read's phased hets
        for (int k = h->frag[f].off; k < h->frag[f].off + h->frag[f].n; ++k) {
            const ha_var_t *x = h->var + h->call[k].v;
            if (x->ps <= 0 || !is_het(x)) continue;
            if (nc == 0) { cur = x->ps; nc = 1; } else if (x->ps == cur) ++nc; else --nc;
        }
        if (cur <= 0) continue;
        hts_pos_t ps = cur; double s0 = 0, s1 = 0; int n = 0;
        for (int k = h->frag[f].off; k < h->frag[f].off + h->frag[f].n; ++k) {
            const ha_call_t *c = h->call + k; const ha_var_t *x = h->var + c->v;
            if (x->ps != cur || !is_het(x)) continue;
            s0 += cterm(c, x->g[0]); s1 += cterm(c, x->g[1]); ++n;
        }
        if (n == 0) continue;
        double p1 = 1.0 / (1.0 + exp(s1 - s0));
        int hap = p1 > o->hap_assign_p ? 1 : (1.0 - p1 > o->hap_assign_p ? 2 : 0);
        if (hap == 0) continue;
        h->frag_hap[f] = (int8_t)hap;
        h->frag_ps[f] = ps;
    }
    ws_free(&w);
}

// ------------------------------------------------------------------ chunks --

void hap_asm_collect_skipped(const call_var_opt_t *opt, bam_chunk_t *chunk) {
    (void)opt;
    if (chunk->digars == NULL || chunk->n_cand_vars <= 0 || chunk->var_i_to_cate == NULL) return;
    int m = 0;
    for (int r = 0; r < chunk->n_reads; ++r) {
        if (chunk->is_skipped[r] != 0 && chunk->is_skipped[r] != BAM_RECORD_WRONG_MAP) continue;
        const digar_t *d = chunk->digars + r;
        // A read in calling has a profile; only the sites it holds no ref/alt
        // call at -- inside its own noisy regions, mostly -- are read here.
        const read_var_profile_t *p1 = (chunk->is_skipped[r] == 0 && chunk->read_var_profile != NULL)
                                     ? chunk->read_var_profile + r : NULL;
        if (d->n_digar <= 0 || d->bseq == NULL) continue;
        // first candidate at or after the read's start
        int lo = 0, hi = chunk->n_cand_vars;
        while (lo < hi) { int mid = (lo + hi) / 2; if (chunk->cand_vars[mid].pos < d->beg) lo = mid + 1; else hi = mid; }
        int gi = 0;
        for (int v = lo; v < chunk->n_cand_vars && chunk->cand_vars[v].pos <= d->end; ++v) {
            const cand_var_t *cv = chunk->cand_vars + v;
            if (cv->var_type != BAM_CDIFF || (chunk->var_i_to_cate[v] & LONGCALLD_CAND_GERMLINE_VAR_CATE) == 0) continue;
            if (cv->alt_seq == NULL) continue;
            if (p1 != NULL && p1->start_var_idx >= 0 && v >= p1->start_var_idx && v <= p1->end_var_idx &&
                v - p1->start_var_idx < p1->cap) {
                int pa = p1->alleles[v - p1->start_var_idx];
                if (pa == 0 || pa == 1) continue;
            }
            while (gi < d->n_digar) { // the digar reaching pos, or the first after it
                const digar1_t *g = d->digars + gi;
                int span = (g->type == BAM_CEQUAL || g->type == BAM_CDIFF || g->type == BAM_CDEL) ? g->len : 0;
                if (g->pos + span - 1 >= cv->pos) break;
                ++gi;
            }
            int a = 0, qi = -1;
            if (gi < d->n_digar && d->digars[gi].pos <= cv->pos) {
                const digar1_t *g = d->digars + gi;
                int off = (int)(cv->pos - g->pos);
                if (g->type == BAM_CEQUAL) { a = 0; qi = g->qi + off; }
                else if (g->type == BAM_CDIFF) { // alt_seq of both is 0..3 (ACGT); bseq is BAM's packed 4-bit
                    qi = g->qi + off;
                    a = (g->alt_seq != NULL && g->alt_seq[off] == cv->alt_seq[0]) ? 1 : -1;
                }
                else a = -1; // deleted, or an insertion at the base
            }
            if (a < 0) continue;
            if (chunk->n_skip_calls == m) { m = m ? m * 2 : 256; chunk->skip_calls = (skip_call_t*)_err_realloc(chunk->skip_calls, m * sizeof(skip_call_t)); }
            skip_call_t *s = chunk->skip_calls + chunk->n_skip_calls++;
            s->read_i = r; s->var_i = v; s->a = (uint8_t)a;
            s->q = (qi >= 0 && qi < d->qlen && d->qual != NULL) ? d->qual[qi] : LONGCALLD_QUAL_UNKNOWN;
        }
    }
}

typedef struct { int gid, v; uint8_t a, ok, tag_only; float lm, lx; } raw_call_t;

static int cmp_raw_call(const void *pa, const void *pb) {
    const raw_call_t *a = (const raw_call_t*)pa, *b = (const raw_call_t*)pb;
    if (a->gid != b->gid) return a->gid < b->gid ? -1 : 1;
    return (a->v > b->v) - (a->v < b->v);
}

static inline double clamp_p(double p) { return p < 1e-6 ? 1e-6 : (p > 0.5 ? 0.5 : p); }

// The quality model's per-read error (var_qual.c): base error toward the one
// other base plus mismapping, and for an indel the repeat-context error in
// each direction on top.
static void call_err(const call_var_opt_t *opt, int q, int mq, int is_indel, double e_r, double e_a,
                     double *p_alt_if_ref, double *p_ref_if_alt) {
    if (q == LONGCALLD_QUAL_UNKNOWN) q = opt->vq_default_bq;
    if (q > opt->vq_bq_cap) q = opt->vq_bq_cap;
    if (mq < 0) mq = 0;
    if (mq > 60) mq = 60;
    double e_bq = pow(10.0, -q / 10.0), e_mq = pow(10.0, -mq / 10.0);
    double e = (1.0 - e_mq) * e_bq / 3.0 + e_mq * 0.5;
    if (!is_indel) { *p_alt_if_ref = *p_ref_if_alt = clamp_p(e); return; }
    *p_alt_if_ref = clamp_p(1.0 - (1.0 - e_r) * (1.0 - e));
    *p_ref_if_alt = clamp_p(1.0 - (1.0 - e_a) * (1.0 - e));
}

// Position a variant is written at, the same rule make_variants() uses to
// decide which chunk owns it.
static inline hts_pos_t record_pos(const cand_var_t *cv) {
    return (cv->var_type == BAM_CDEL || cv->var_type == BAM_CINS) ? cv->pos - 1 : cv->pos;
}

static void hap_asm_group(const call_var_opt_t *opt, bam_chunk_t *chunks, int beg, int end) {
    int n_ch = end - beg, paired = opt->n_samples > 1;
    int **vmap = (int**)_err_calloc(n_ch, sizeof(int*)), **gid = (int**)_err_calloc(n_ch, sizeof(int*));
    ha_t h; memset(&h, 0, sizeof(h));
    int m_var = 0; double *er = NULL, *ea = NULL;
    // variants: every germline candidate a chunk owns
    for (int k = 0; k < n_ch; ++k) {
        bam_chunk_t *c = chunks + beg + k;
        vmap[k] = (int*)_err_malloc((c->n_cand_vars > 0 ? c->n_cand_vars : 1) * sizeof(int));
        for (int i = 0; i < c->n_cand_vars; ++i) {
            vmap[k][i] = -1;
            const cand_var_t *cv = c->cand_vars + i;
            if (c->var_i_to_cate == NULL || (c->var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE) == 0) continue;
            if (cv->hap_to_cons_alle == NULL) continue;
            hts_pos_t rp = record_pos(cv);
            if (rp < c->reg_beg || rp > c->reg_end) continue;
            if (h.n_var == m_var) {
                m_var = m_var ? m_var * 2 : 1024;
                h.var = (ha_var_t*)_err_realloc(h.var, m_var * sizeof(ha_var_t));
                er = (double*)_err_realloc(er, m_var * sizeof(double));
                ea = (double*)_err_realloc(ea, m_var * sizeof(double));
            }
            ha_var_t *x = h.var + h.n_var; memset(x, 0, sizeof(*x));
            x->pos = rp; x->is_indel = cv->var_type != BAM_CDIFF; x->phasable = 1;
            x->geno = (c->var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_CLEAN_VAR_CATE) ? 2 : 1;
            int h1 = cv->hap_to_cons_alle[1], h2 = cv->hap_to_cons_alle[2];
            if (h1 == -1 && h2 == -1) h1 = h2 = cv->hap_to_cons_alle[0];
            x->g[0] = h1 == 1; x->g[1] = h2 == 1;
            er[h.n_var] = ea[h.n_var] = 0;
            if (x->is_indel) {
                uint8_t per, tr; char unit[8];
                ref_repeat_context(c->ref_seq, c->ref_beg, c->ref_end, cv->pos, 0, &per, &tr, unit);
                indel_ctx_error(opt->is_ont, per, tr, er + h.n_var, ea + h.n_var);
                // Only an indel outside any repeat, called in a clean region, may
                // join the phasing graph (Longshot phases substitutions alone). A
                // homopolymer or tandem-repeat indel is phased against its block
                // but cannot tie two blocks together: on HG002 ONT chr20 runs of
                // such GQ ~1 "hets" carried a 121 kb block across a switch.
                double e_r0, e_a0;
                indel_ctx_error(opt->is_ont, 0, 1, &e_r0, &e_a0);
                x->phasable = x->geno == 2 && ea[h.n_var] <= e_a0;
            }
            vmap[k][i] = h.n_var++;
        }
    }
    // reads: one fragment per read, its copies in adjacent chunks joined
    int n_gid = 0, m_raw = 0, n_raw = 0, m_gid = 0; raw_call_t *raw = NULL; uint8_t *gid_smp = NULL;
    double max_mis = pow(10.0, -opt->min_allele_qual / 10.0);
    for (int k = 0; k < n_ch; ++k) {
        bam_chunk_t *c = chunks + beg + k;
        gid[k] = (int*)_err_malloc((c->n_reads > 0 ? c->n_reads : 1) * sizeof(int));
        for (int r = 0; r < c->n_reads; ++r) gid[k][r] = -1;
        if (k > 0) {
            bam_chunk_t *p = c - 1;
            for (int b = 0; b < opt->n_in_bam_fn; ++b) {
                if (c->n_up_ovlp_reads == NULL || p->n_down_ovlp_reads == NULL) break;
                if (c->n_up_ovlp_reads[b] != p->n_down_ovlp_reads[b]) continue;
                for (int j = 0; j < c->n_up_ovlp_reads[b]; ++j)
                    gid[k][c->up_ovlp_read_i[b][j]] = gid[k-1][p->down_ovlp_read_i[b][j]];
            }
        }
        for (int r = 0; r < c->n_reads; ++r) {
            if (gid[k][r] < 0) {
                if (n_gid == m_gid) { m_gid = m_gid ? m_gid * 2 : 4096; gid_smp = (uint8_t*)_err_realloc(gid_smp, m_gid); }
                gid_smp[n_gid] = c->read_sample_idx != NULL ? c->read_sample_idx[r] : 0;
                gid[k][r] = n_gid++;
            }
            if (c->is_skipped[r] || c->read_var_profile == NULL) continue;
            const read_var_profile_t *p1 = c->read_var_profile + r;
            if (p1->start_var_idx < 0 || p1->end_var_idx < 0) continue;
            int mq = c->read_mapq != NULL ? c->read_mapq[r] : 60;
            for (int vi = p1->start_var_idx; vi <= p1->end_var_idx && vi < c->n_cand_vars; ++vi) {
                int j = vi - p1->start_var_idx;
                if (j >= p1->cap) break;
                int a = p1->alleles[j], gv = vmap[k][vi];
                if ((a != 0 && a != 1) || gv < 0) continue;
                int q = p1->qual != NULL ? p1->qual[j] : LONGCALLD_QUAL_UNKNOWN;
                double par, pra;
                call_err(opt, q, mq, h.var[gv].is_indel, er[gv], ea[gv], &par, &pra);
                double lm = a ? log(1.0 - pra) : log(1.0 - par), lx = a ? log(par) : log(pra);
                if (n_raw == m_raw) { m_raw = m_raw ? m_raw * 2 : 65536; raw = (raw_call_t*)_err_realloc(raw, m_raw * sizeof(raw_call_t)); }
                raw_call_t *rc = raw + n_raw++;
                rc->gid = gid[k][r]; rc->v = gv; rc->a = (uint8_t)a; rc->lm = (float)lm; rc->lx = (float)lx;
                // Trusted: right with probability above the floor, on a well-mapped read.
                double eps = 1.0 / (1.0 + exp(lm - lx));
                rc->ok = eps < max_mis && mq >= opt->min_phase_mapq;
                rc->tag_only = 0;
            }
        }
        for (int i = 0; i < c->n_skip_calls; ++i) {
            const skip_call_t *s = c->skip_calls + i;
            int gv = s->var_i < c->n_cand_vars ? vmap[k][s->var_i] : -1;
            if (gv < 0) continue;
            int mq = c->read_mapq != NULL ? c->read_mapq[s->read_i] : 60;
            double par, pra;
            call_err(opt, s->q, mq, 0, 0, 0, &par, &pra);
            if (n_raw == m_raw) { m_raw = m_raw ? m_raw * 2 : 65536; raw = (raw_call_t*)_err_realloc(raw, m_raw * sizeof(raw_call_t)); }
            raw_call_t *rc = raw + n_raw++;
            rc->gid = gid[k][s->read_i]; rc->v = gv; rc->a = s->a; rc->ok = 0; rc->tag_only = 1;
            rc->lm = (float)(s->a ? log(1.0 - pra) : log(1.0 - par)); rc->lx = (float)(s->a ? log(par) : log(pra));
        }
    }
    if (raw != NULL) qsort(raw, n_raw, sizeof(raw_call_t), cmp_raw_call);
    int *fmap = (int*)_err_malloc((n_gid > 0 ? n_gid : 1) * sizeof(int));
    for (int g = 0; g < n_gid; ++g) fmap[g] = -1;
    h.call = (ha_call_t*)_err_malloc((n_raw > 0 ? n_raw : 1) * sizeof(ha_call_t));
    h.frag = (ha_frag_t*)_err_malloc((n_gid > 0 ? n_gid : 1) * sizeof(ha_frag_t));
    for (int i = 0; i < n_raw; ++i) {
        if (i > 0 && raw[i].gid == raw[i-1].gid && raw[i].v == raw[i-1].v) continue; // a read seen twice
        if (fmap[raw[i].gid] < 0) {
            fmap[raw[i].gid] = h.n_frag;
            h.frag[h.n_frag].off = h.n_call; h.frag[h.n_frag].n = 0; h.frag[h.n_frag].smp = gid_smp[raw[i].gid];
            ++h.n_frag;
        }
        ha_call_t *hc = h.call + h.n_call++;
        hc->v = raw[i].v; hc->a = raw[i].a; hc->ok = raw[i].ok; hc->tag_only = raw[i].tag_only;
        hc->lm = raw[i].lm; hc->lx = raw[i].lx;
        h.frag[fmap[raw[i].gid]].n++;
    }
    free(raw);
    // Paired: only a site that looks like a germline het in the normal may
    // join the phasing graph, so tumour reads link hets but cannot make them.
    if (paired) {
        int *ndp = (int*)_err_calloc(h.n_var > 0 ? h.n_var : 1, sizeof(int));
        int *nalt = (int*)_err_calloc(h.n_var > 0 ? h.n_var : 1, sizeof(int));
        for (int f = 0; f < h.n_frag; ++f) {
            if (h.frag[f].smp != LONGCALLD_SAMPLE_NORMAL) continue;
            for (int i = h.frag[f].off; i < h.frag[f].off + h.frag[f].n; ++i) {
                if (h.call[i].tag_only) continue;
                ndp[h.call[i].v]++; nalt[h.call[i].v] += h.call[i].a;
            }
        }
        for (int v = 0; v < h.n_var; ++v) {
            double af = ndp[v] > 0 ? (double)nalt[v] / ndp[v] : 0;
            h.var[v].phasable = ndp[v] >= opt->min_normal_phase_dp && af >= opt->min_normal_phase_af && af <= opt->max_normal_phase_af;
        }
        free(ndp); free(nalt);
    }
    ha_opt_t ho; ha_opt_init(&ho);
    ho.max_miscall = max_mis;
    ho.min_pq = opt->min_phase_qual;
    ho.hap_assign_p = 1.0 - pow(10.0, -opt->hap_assign_qual / 10.0);
    ho.geno_smp = paired ? LONGCALLD_SAMPLE_NORMAL : -1;
    ho.solo_needs_phasable = paired;
    ha_run(&ho, &h);
    // Write back: genotypes and blocks onto the candidates, haplotypes onto the reads.
    for (int k = 0; k < n_ch; ++k) {
        bam_chunk_t *c = chunks + beg + k;
        if (c->n_cand_vars > 0) {
            free(c->var_hap_flag); free(c->var_sq); free(c->var_pq);
            c->var_hap_flag = (uint8_t*)_err_calloc(c->n_cand_vars, 1);
            c->var_sq = (float*)_err_malloc(c->n_cand_vars * sizeof(float));
            c->var_pq = (float*)_err_malloc(c->n_cand_vars * sizeof(float));
        }
        for (int i = 0; i < c->n_cand_vars; ++i) {
            c->var_sq[i] = c->var_pq[i] = -1;
            int gv = vmap[k][i];
            if (gv < 0) continue;
            const ha_var_t *x = h.var + gv; cand_var_t *cv = c->cand_vars + i;
            cv->hap_to_cons_alle[1] = x->g[0]; cv->hap_to_cons_alle[2] = x->g[1];
            cv->phase_set = (x->g[0] != x->g[1] && x->ps > 0) ? x->ps : 0;
            c->var_hap_flag[i] = x->inconsistent ? LONGCALLD_HAP_INCONSISTENT : 0;
            c->var_sq[i] = x->sq; c->var_pq[i] = x->pq;
        }
        for (int r = 0; r < c->n_reads; ++r) {
            int f = fmap[gid[k][r]];
            // A skipped read keeps a tag here only through its tag-only calls;
            // the calling code tests is_skipped before it looks at haps[].
            if (f < 0 || h.frag_hap[f] == 0) { c->haps[r] = 0; c->phase_sets[r] = -1; continue; }
            c->haps[r] = h.frag_hap[f]; c->phase_sets[r] = h.frag_ps[f];
        }
        free(vmap[k]); free(gid[k]);
    }
    free(vmap); free(gid); free(fmap); free(gid_smp); free(er); free(ea);
    ha_free(&h);
}

void hap_asm_chunks(const call_var_opt_t *opt, bam_chunk_t *chunks, int n_chunks) {
    for (int s = 0; s < n_chunks; ) {
        int e = s + 1;
        while (e < n_chunks && chunks[e].tid == chunks[s].tid) ++e;
        hap_asm_group(opt, chunks, s, e);
        s = e;
    }
}
