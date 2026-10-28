#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <pthread.h>
#include "main.h"
#include "utils.h"
#include "bam_utils.h"
#include "collect_var.h"
#include "call_var_main.h"
#include "var_qual.h"

// Index of each fraction is fixed: the somatic models below address the grid
// by position, and every haplotype fraction phi there has phi/2 on the grid
// for the unassigned reads, which carry half of it.
const double LONGCALLD_QGRID[LONGCALLD_QGRID_N] = {
    0.0, 0.005, 0.01, 0.02, 0.025, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2,
    0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0
};
#define QG_0    0
#define QG_CONT 3   // 0.02: tumour-in-normal contamination
#define QG_P10  8   // 0.1 and 0.9: a haplotype carrying an allele, up to mis-phased reads
#define QG_HALF 16
#define QG_P90  20
#define QG_1    21
// Somatic within-haplotype fractions, and the index of half of each.
static const int QG_PHI[]      = { 2, 3, 5, 6, 8, 9, 10, 12, 14, 16, 17, 18, 19, 20, 21 };
static const int QG_PHI_HALF[] = { 1, 2, 3, 4, 6, 7,  8,  9, 10, 11, 12, 13, 14, 15, 16 };
#define QG_N_PHI ((int)(sizeof(QG_PHI) / sizeof(QG_PHI[0])))

const char *LONGCALLD_VQF_NAMES[VQF_N] = {
    "bias", "sor", "one_strand", "mq_rs", "mq_diff", "bq_rs", "bq_low",
    "hp_indel", "str_indel", "hp_snv", "unph_frac", "unph_alt_excess", "log_unph",
    "alt_minor_hap", "alt_conflict", "unphased_call", "dp_high", "dp_low",
    "log_alt", "vaf", "lik", "lik_low"
};
const char *LONGCALLD_VQ_CLASS_NAMES[VQ_N_CLASS] = { "germ_snv", "germ_indel", "som_snv", "som_indel" };

static inline double lse2(double a, double b) { // log10(10^a + 10^b)
    if (a < b) { double t = a; a = b; b = t; }
    if (b == -INFINITY) return a;
    return a + log10(1.0 + pow(10.0, b - a));
}

// ---------------------------------------------------------------- step 1 --

// Read position of reference base `pos`, or -1 if the read deletes it. Walks
// forward from *hint, so the caller must ask for non-decreasing positions.
static int digar_qi_at(const digar_t *d, hts_pos_t pos, int *hint) {
    int i = *hint;
    for (; i < d->n_digar; ++i) {
        const digar1_t *g = d->digars + i;
        // I, S and H take no reference bases: they end before they start.
        int takes_ref = (g->type == BAM_CEQUAL || g->type == BAM_CDIFF || g->type == BAM_CDEL);
        hts_pos_t g_end = takes_ref ? g->pos + g->len - 1 : g->pos - 1;
        if (g_end >= pos) break;
    }
    *hint = i;
    if (i >= d->n_digar) return -1;
    const digar1_t *g = d->digars + i;
    if (g->pos > pos) return -1;
    if (g->type == BAM_CEQUAL || g->type == BAM_CDIFF) return g->qi + (int)(pos - g->pos);
    return -1;
}

static inline int qual_at(const digar_t *d, int qi) {
    if (qi < 0 || qi >= d->qlen) return -1;
    return d->qual[qi];
}

// Mean of the two bases either side of an indel's event point, the bases an
// aligner can trade for the gap. One side missing: the other alone.
static int flank_qual(const digar_t *d, int qi) {
    int a = qual_at(d, qi - 1), b = qual_at(d, qi);
    if (a < 0) return b;
    if (b < 0) return a;
    return (a + b) / 2;
}

// The quality an allele call rests on. Mirrors get_alt_qual() in assign_hap.c,
// without its allocation: substitution, the base; insertion, the median
// inserted base; deletion, the flanking pair. A reference read at an indel is
// judged on the same flanking pair the alt would have been.
static int alt_cell_qual(const digar_t *d, const cand_var_t *cv, int alt_qi) {
    if (alt_qi < 0) return -1;
    if (cv->var_type == BAM_CDIFF) return qual_at(d, alt_qi);
    if (cv->var_type == BAM_CDEL) return flank_qual(d, alt_qi);
    int n = cv->alt_len > 0 ? cv->alt_len : 1, buf[64];
    if (n > 64) n = 64;
    for (int i = 0; i < n; ++i) { int q = qual_at(d, alt_qi + i); buf[i] = q < 0 ? 0 : q; }
    // insertion sort; n is small
    for (int i = 1; i < n; ++i) { int x = buf[i], j = i - 1; while (j >= 0 && buf[j] > x) { buf[j+1] = buf[j]; --j; } buf[j+1] = x; }
    return buf[(n - 1) / 2];
}

void fill_read_var_profile_qual(const call_var_opt_t *opt, bam_chunk_t *chunk) {
    (void)opt;
    read_var_profile_t *p = chunk->read_var_profile;
    if (p == NULL || chunk->digars == NULL) return;
    for (int k = 0; k < chunk->n_reads; ++k) {
        int read_i = chunk->ordered_read_ids[k];
        if (chunk->is_skipped[read_i]) continue;
        read_var_profile_t *p1 = p + read_i;
        if (p1->start_var_idx < 0 || p1->end_var_idx < 0 || p1->qual == NULL) continue;
        const digar_t *d = chunk->digars + read_i;
        if (d->qual == NULL) continue;
        int hint = 0;
        for (int var_i = p1->start_var_idx; var_i <= p1->end_var_idx; ++var_i) {
            int j = var_i - p1->start_var_idx;
            if (j >= p1->cap || var_i >= chunk->n_cand_vars) break;
            const cand_var_t *cv = chunk->cand_vars + var_i;
            int a = p1->alleles[j], q = -1;
            if (a == 1 || a == -2) q = alt_cell_qual(d, cv, p1->alt_qi[j]);
            else if (a == 0) {
                int qi = digar_qi_at(d, cv->pos, &hint);
                if (qi >= 0) q = cv->var_type == BAM_CDIFF ? qual_at(d, qi) : flank_qual(d, qi);
            }
            p1->qual[j] = (q < 0) ? LONGCALLD_QUAL_UNKNOWN : (uint8_t)(q > 254 ? 254 : q);
        }
    }
}

// --------------------------------------------------------- repeat context --

void ref_repeat_context(const char *ref_seq, hts_pos_t beg, hts_pos_t end, hts_pos_t pos, int is_snv,
                        uint8_t *period, uint8_t *tract_len, char *unit) {
    *period = 0; *tract_len = 1; unit[0] = '\0';
    if (ref_seq == NULL || pos < beg || pos > end) return;
    const char *s = ref_seq - beg; // s[x] is the base at reference position x
    // Homopolymer at pos; for a substitution also either neighbour, since a
    // base beside a run is where run-length errors masquerade as mismatches.
    int best_hp = 0; hts_pos_t hp_at = pos;
    for (int off = is_snv ? -1 : 0; off <= (is_snv ? 1 : 0); ++off) {
        hts_pos_t x = pos + off;
        if (x < beg || x > end) continue;
        hts_pos_t a = x, z = x;
        while (a > beg && toupper(s[a-1]) == toupper(s[x])) --a;
        while (z < end && toupper(s[z+1]) == toupper(s[x])) ++z;
        if (z - a + 1 > best_hp) { best_hp = (int)(z - a + 1); hp_at = x; }
    }
    // Shortest-period tandem repeat through pos with the longest tract.
    int best_tr = 0, best_p = 0;
    for (int p = 2; p <= 6; ++p) {
        if (pos + p - 1 > end) break;
        int copies = 1;
        for (hts_pos_t j = pos + p; j + p - 1 <= end; j += p) {
            int same = 1;
            for (int k = 0; k < p; ++k) if (toupper(s[j+k]) != toupper(s[pos+k])) { same = 0; break; }
            if (!same) break;
            ++copies;
        }
        for (hts_pos_t j = pos - p; j >= beg; j -= p) {
            int same = 1;
            for (int k = 0; k < p; ++k) if (toupper(s[j+k]) != toupper(s[pos+k])) { same = 0; break; }
            if (!same) break;
            ++copies;
        }
        // A period-p tract made of one base is a homopolymer, already counted.
        int uniform = 1;
        for (int k = 1; k < p; ++k) if (toupper(s[pos+k]) != toupper(s[pos])) { uniform = 0; break; }
        if (uniform) continue;
        if (copies >= 2 && copies * p > best_tr) { best_tr = copies * p; best_p = p; }
    }
    if (best_hp >= 3 && best_hp >= best_tr) {
        *period = 1; *tract_len = (uint8_t)(best_hp > 255 ? 255 : best_hp);
        unit[0] = (char)toupper(s[hp_at]); unit[1] = '\0';
    } else if (best_p > 0) {
        *period = (uint8_t)best_p; *tract_len = (uint8_t)(best_tr > 255 ? 255 : best_tr);
        for (int k = 0; k < best_p; ++k) unit[k] = (char)toupper(s[pos+k]);
        unit[best_p] = '\0';
    }
}

// Per-read indel error at this repeat context, in two directions that differ
// by an order of magnitude: e_a, a read from the indel allele showing the
// reference (slippage undoing the event), and e_r, a reference read showing
// this particular indel. Measured on GIAB HG002 chr20 against v4.2.1:
// hom-alt VAF = 1 - e_a, het VAF = (1 - e_a + e_r) / 2. In HiFi homopolymers
// e_a climbs from ~2% at runs of 5 to 9% at 10, 15% at 15 and ~25% at 20,
// while e_r stays at 1-4%; ONT has e_a ~6% to a run of 8, then 25-45%.
// A single symmetric rate, as here before, took a T15 het at VAF 0.31 for
// reference plus slippage and gave it QUAL 0.
static void indel_ctx_error(int is_ont, int period, int tract, double *e_r, double *e_a) {
    int L = tract;
    if (period == 1 && L >= 3) {
        if (!is_ont) { *e_a = 0.01 + 0.0135 * (L > 4 ? L - 4 : 0); *e_r = 0.001 + 0.002 * (L > 4 ? L - 4 : 0); }
        else { *e_a = 0.06 + 0.08 * (L > 8 ? L - 8 : 0); *e_r = 0.005 + 0.02 * (L > 8 ? L - 8 : 0); }
    } else if (period == 2 && L >= 4) {
        if (!is_ont) { *e_a = 0.015 + 0.0025 * (L - 4); *e_r = 0.001 + 0.0005 * (L - 4); }
        else { *e_a = 0.04 + 0.006 * (L - 4); *e_r = 0.003 + 0.002 * (L - 4); }
    } else if (period >= 3 && L >= 6) {
        if (!is_ont) { *e_a = 0.02 + 0.001 * (L - 6); *e_r = 0.001; }
        else { *e_a = 0.06 + 0.002 * (L - 6); *e_r = 0.003; }
    } else {
        *e_a = is_ont ? 0.045 : 0.01; *e_r = is_ont ? 0.002 : 0.0003;
    }
    if (*e_a > 0.5) *e_a = 0.5;
    if (*e_r > 0.15) *e_r = 0.15;
}

// ---------------------------------------------------------------- step 2 --

var_qual_acc_t *collect_var_qual_acc(const call_var_opt_t *opt, bam_chunk_t *chunk, const int *want) {
    if (chunk->n_cand_vars <= 0 || chunk->read_var_profile == NULL) return NULL;
    int n = chunk->n_cand_vars;
    var_qual_acc_t *acc = (var_qual_acc_t*)_err_calloc(n, sizeof(var_qual_acc_t));
    // Per-variant constant: an indel's context error, both directions. Applied
    // when the grid is filled, so indel reads bin on their base error alone.
    for (int i = 0; i < n; ++i) {
        if (want != NULL && !want[i]) continue;
        const cand_var_t *cv = chunk->cand_vars + i;
        if (cv->var_type == BAM_CDIFF) continue;
        uint8_t per, tr; char unit[8]; double e_r, e_a;
        ref_repeat_context(chunk->ref_seq, chunk->ref_beg, chunk->ref_end, cv->pos, 0, &per, &tr, unit);
        indel_ctx_error(opt->is_ont, per, tr, &e_r, &e_a);
        acc[i].is_indel = 1; acc[i].e_r = (float)e_r; acc[i].e_a = (float)e_a;
    }
    int n_samples = chunk->n_samples > 0 ? chunk->n_samples : 1;
    double bq_cap = opt->vq_bq_cap;
    // Substitution error phred for every (base quality, mapping quality) pair,
    // so the per-cell work is a lookup; indels depend on their context too.
    static const int NQ = LONGCALLD_QHIST_N, NM = 61;
    uint8_t snv_qe[LONGCALLD_QHIST_N][61];
    double e_bq_tab[LONGCALLD_QHIST_N], e_mq_tab[61];
    for (int q = 0; q < NQ; ++q) e_bq_tab[q] = pow(10.0, -(q < bq_cap ? q : bq_cap) / 10.0);
    for (int m = 0; m < NM; ++m) e_mq_tab[m] = pow(10.0, -m / 10.0);
    for (int q = 0; q < NQ; ++q)
        for (int m = 0; m < NM; ++m) {
            double e = (1.0 - e_mq_tab[m]) * e_bq_tab[q] / 3.0 + e_mq_tab[m] * 0.5;
            int qe = (int)(-10.0 * log10(e) + 0.5);
            snv_qe[q][m] = (uint8_t)(qe < 0 ? 0 : (qe >= LONGCALLD_QE_N ? LONGCALLD_QE_N - 1 : qe));
        }
    for (int k = 0; k < chunk->n_reads; ++k) {
        int read_i = chunk->ordered_read_ids[k];
        if (chunk->is_skipped[read_i]) continue;
        read_var_profile_t *p1 = chunk->read_var_profile + read_i;
        if (p1->start_var_idx < 0 || p1->end_var_idx < 0) continue;
        int smp = chunk->read_sample_idx != NULL ? chunk->read_sample_idx[read_i] : 0;
        if (smp < 0 || smp >= n_samples || smp >= 2) continue;
        // The same slot rule collect_all_var_slot_counts() applies, so the
        // likelihood is decided on exactly the grid the classifier saw.
        int hap = chunk->haps[read_i];
        if (hap < 0 || hap >= LONGCALLD_N_HAP_SLOTS) hap = 0;
        if (chunk->phase_sets[read_i] <= 0) hap = 0;
        hts_pos_t ps = hap != 0 ? chunk->phase_sets[read_i] : 0;
        int mq = chunk->read_mapq != NULL ? chunk->read_mapq[read_i] : 60;
        if (mq > 60) mq = 60;
        int is_rev = chunk->read_is_rev != NULL ? chunk->read_is_rev[read_i] : 0;
        int conflict = chunk->n_clean_conflict_snps != NULL && chunk->n_clean_conflict_snps[read_i] >= 2;
        for (int var_i = p1->start_var_idx; var_i <= p1->end_var_idx; ++var_i) {
            if (var_i < 0 || var_i >= n) continue;
            if (want != NULL && !want[var_i]) continue;
            int j = var_i - p1->start_var_idx;
            int a = p1->alleles[j];
            if (a != 0 && a != 1) continue; // third allele, or alt failing the base-quality floor
            int h = (hap != 0 && chunk->var_slots != NULL && ps == chunk->var_slots[var_i].ps) ? hap : 0;
            int q = p1->qual != NULL ? p1->qual[j] : LONGCALLD_QUAL_UNKNOWN;
            int known = q != LONGCALLD_QUAL_UNKNOWN;
            if (!known) q = opt->vq_default_bq;
            int qb = q < NQ ? q : NQ - 1;
            // Base error e_bq/3 (the one specific other base) and mismapping; an
            // indel's context error is composed with it in fill_ll().
            int qe = snv_qe[qb][mq < 0 ? 0 : mq];
            if (acc[var_i].ne[smp][h][a][qe] < UINT16_MAX) acc[var_i].ne[smp][h][a][qe]++;
            acc[var_i].n_strand[smp][a][is_rev ? 1 : 0]++;
            if (h == 0) acc[var_i].n_unph[smp][a]++;
            if (a == 1 && conflict) acc[var_i].n_alt_conflict[smp]++;
            if (known) { int b = q < LONGCALLD_QHIST_N ? q : LONGCALLD_QHIST_N - 1; if (acc[var_i].bq[smp][a][b] < UINT16_MAX) acc[var_i].bq[smp][a][b]++; }
            { int b = mq < LONGCALLD_QHIST_N ? mq : LONGCALLD_QHIST_N - 1; if (acc[var_i].mq[smp][a][b] < UINT16_MAX) acc[var_i].mq[smp][a][b]++; }
        }
    }
    return acc;
}

// Mann-Whitney z of the alt histogram against the ref one, with the tie
// correction -- MAPQ is mostly one value, and without it the variance is
// overstated and every z collapses towards 0. Negative: alt ranks lower.
static double rank_sum_z(const int *alt, const int *ref) {
    double n1 = 0, n2 = 0;
    for (int b = 0; b < LONGCALLD_QHIST_N; ++b) { n1 += alt[b]; n2 += ref[b]; }
    if (n1 < 1 || n2 < 1) return 0.0;
    double u = 0, ref_below = 0, ties = 0, n = n1 + n2;
    for (int b = 0; b < LONGCALLD_QHIST_N; ++b) {
        u += alt[b] * (ref_below + 0.5 * ref[b]);
        ref_below += ref[b];
        double t = alt[b] + ref[b];
        ties += t * t * t - t;
    }
    double var = n1 * n2 / 12.0 * ((n + 1) - ties / (n * (n - 1)));
    if (var <= 1e-9) return 0.0;
    return (u - n1 * n2 / 2.0) / sqrt(var);
}

static double hist_mean(const int *h) {
    double s = 0, n = 0;
    for (int b = 0; b < LONGCALLD_QHIST_N; ++b) { s += (double)b * h[b]; n += h[b]; }
    return n > 0 ? s / n : 0.0;
}

static int hist_median(const int *h) {
    int n = 0, c = 0;
    for (int b = 0; b < LONGCALLD_QHIST_N; ++b) n += h[b];
    if (n == 0) return -1;
    for (int b = 0; b < LONGCALLD_QHIST_N; ++b) { c += h[b]; if (2 * c >= n) return b; }
    return LONGCALLD_QHIST_N - 1;
}

// log10 P(read | f) for a read of error phred qe showing ref ([0]) or alt ([1]).
static double qe_ll_tab[LONGCALLD_QE_N][LONGCALLD_QGRID_N][2];
static pthread_once_t qe_ll_once = PTHREAD_ONCE_INIT;
static void init_qe_ll_tab(void) {
    for (int qe = 0; qe < LONGCALLD_QE_N; ++qe) {
        double e = pow(10.0, -qe / 10.0);
        if (e > 0.5) e = 0.5;
        for (int g = 0; g < LONGCALLD_QGRID_N; ++g) {
            double f = LONGCALLD_QGRID[g];
            qe_ll_tab[qe][g][0] = log10(f * e + (1.0 - f) * (1.0 - e));
            qe_ll_tab[qe][g][1] = log10(f * (1.0 - e) + (1.0 - f) * e);
        }
    }
}

static void fill_ll(var_qual_acc_t *acc) {
    pthread_once(&qe_ll_once, init_qe_ll_tab);
    for (int s = 0; s < 2; ++s)
        for (int h = 0; h < 3; ++h) {
            double ll[LONGCALLD_QGRID_N] = {0};
            for (int a = 0; a < 2; ++a)
                for (int qe = 0; qe < LONGCALLD_QE_N; ++qe) {
                    int n = acc->ne[s][h][a][qe];
                    if (n == 0) continue;
                    if (!acc->is_indel) {
                        for (int g = 0; g < LONGCALLD_QGRID_N; ++g) ll[g] += n * qe_ll_tab[qe][g][a];
                        continue;
                    }
                    // Indel: the context errors, one per direction, on top of the
                    // read's own base and mapping error.
                    double eb = pow(10.0, -qe / 10.0);
                    double p_alt_if_ref = 1.0 - (1.0 - acc->e_r) * (1.0 - eb);
                    double p_ref_if_alt = 1.0 - (1.0 - acc->e_a) * (1.0 - eb);
                    for (int g = 0; g < LONGCALLD_QGRID_N; ++g) {
                        double f = LONGCALLD_QGRID[g];
                        double p = (a == 1) ? f * (1.0 - p_ref_if_alt) + (1.0 - f) * p_alt_if_ref
                                            : f * p_ref_if_alt + (1.0 - f) * (1.0 - p_alt_if_ref);
                        ll[g] += n * log10(p);
                    }
                }
            for (int g = 0; g < LONGCALLD_QGRID_N; ++g) acc->ll[s][h][g] = (float)ll[g];
        }
}

static inline double slot_ll(const var_qual_acc_t *acc, int s, int h, int g) { return acc->ll[s][h][g]; }
static inline double smp_ll(const var_qual_acc_t *acc, int s, int g) {
    return acc->ll[s][0][g] + acc->ll[s][1][g] + acc->ll[s][2][g];
}

// Diploid genotype log-likelihoods of one sample. The het term mixes two
// readings of the same reads: phased, where each haplotype is pure ref or pure
// alt and only the unassigned reads sit at 1/2, and unphased, every read at
// 1/2. A clean segregation makes the phased reading win by up to 2^n; alt
// scattered across both haplotypes makes it lose, and the unphased reading
// carries the het. Phasing is read off the same reads, so it is mixed rather
// than trusted outright. `any_frac`: also let the het be any fraction up to 1/2
// -- a tumour, or a mosaic sample, is not diploid at the site.
static void sample_gl(const call_var_opt_t *opt, const var_qual_acc_t *acc, int s, int grid_ok, int any_frac,
                      double gl[3]) {
    gl[0] = smp_ll(acc, s, QG_0);
    gl[2] = smp_ll(acc, s, QG_1);
    double het = smp_ll(acc, s, QG_HALF);
    int n_phased = acc->n_strand[s][0][0] + acc->n_strand[s][0][1] + acc->n_strand[s][1][0] + acc->n_strand[s][1][1]
                 - acc->n_unph[s][0] - acc->n_unph[s][1];
    if (grid_ok && n_phased > 0) {
        double o1 = slot_ll(acc, s, 1, QG_1) + slot_ll(acc, s, 2, QG_0) + slot_ll(acc, s, 0, QG_HALF);
        double o2 = slot_ll(acc, s, 1, QG_0) + slot_ll(acc, s, 2, QG_1) + slot_ll(acc, s, 0, QG_HALF);
        double ph = lse2(o1, o2) - log10(2.0);
        double w = opt->vq_phase_weight;
        het = lse2(ph + log10(w), het + log10(1.0 - w));
    }
    if (any_frac)
        for (int g = 1; g <= QG_HALF; ++g) { double l = smp_ll(acc, s, g); if (l > het) het = l; }
    gl[1] = het;
}

// PL is tempered by `scale` (< 1): reads are not independent -- a shared
// mapping or a mis-phased block moves many at once -- so the raw genotype
// likelihoods are overconfident by a roughly constant factor, measured against
// GIAB genotypes by `eval/qual_eval.py gqfit`. Tempering keeps the order of the
// genotypes, so it never changes which one is most likely.
static void gl_to_pl(const double gl[3], double scale, int16_t pl[3]) {
    double mx = gl[0];
    if (gl[1] > mx) mx = gl[1];
    if (gl[2] > mx) mx = gl[2];
    for (int g = 0; g < 3; ++g) {
        double v = -10.0 * scale * (gl[g] - mx);
        pl[g] = (int16_t)(v > 9999 ? 9999 : (int)(v + 0.5));
    }
}

int var_qual_gq(const int16_t pl[3], int dosage) {
    if (dosage < 0 || dosage > 2) return 0;
    int other = INT32_MAX;
    for (int g = 0; g < 3; ++g) if (g != dosage && pl[g] < other) other = pl[g];
    int gq = other - pl[dosage];
    return gq < 0 ? 0 : (gq > 99 ? 99 : gq);
}

static double germline_prior(const call_var_opt_t *opt, const var1_t *v) {
    if (v->type == BAM_CDIFF) return opt->vq_theta_snv;
    // Tandem repeats mutate orders of magnitude faster than unique sequence.
    int tr = v->vq.rep_len, per = v->vq.rep_period;
    if ((per == 1 && tr >= 6) || (per >= 2 && tr >= 2 * per + 4)) return opt->vq_theta_indel * 10.0;
    return opt->vq_theta_indel;
}

// log10 P(no variant | D) for a germline record, against het and hom-alt.
static double germline_lfalse(const call_var_opt_t *opt, const var_qual_acc_t *acc, const var1_t *v, int n_samples,
                              int grid_ok) {
    double gl[3] = {0, 0, 0}, g1[3];
    for (int s = 0; s < n_samples; ++s) {
        sample_gl(opt, acc, s, grid_ok, n_samples > 1 && s == LONGCALLD_SAMPLE_TUMOR, g1);
        for (int g = 0; g < 3; ++g) gl[g] += g1[g];
    }
    double th = germline_prior(opt, v);
    double a0 = log10(1.0 - 1.5 * th) + gl[0], a1 = log10(th) + gl[1], a2 = log10(th / 2.0) + gl[2];
    return a0 - lse2(lse2(a0, a1), a2);
}

// log10 P(not somatic | D). Hypotheses:
//   absent:    every alt read in either sample is error
//   germline:  normal het or hom; tumour as for somatic, or at any fraction
//   somatic:   normal clean or lightly contaminated by tumour; tumour either on
//              one haplotype (h at within-haplotype fraction phi, the other
//              haplotype clean, unassigned reads at phi/2) or at any pooled
//              fraction up to 1
// The pooled branch is not optional. Without it a somatic call at tumour VAF 1
// -- copy-neutral loss of heterozygosity, common in HCC1395 -- could only be
// germline, and a call whose alt reads the phasing split over both haplotypes
// could only be absent: both came out at QUAL 0 with 40 clean alt reads
// against a clean normal. The tumour cannot separate somatic from germline in
// any case; only the normal can, so both hypotheses read the tumour alike.
// Single-sample (mosaic) calls drop the normal terms; their germline
// alternative is the diploid one.
static double somatic_lfalse(const call_var_opt_t *opt, const var_qual_acc_t *acc, const var1_t *v, int n_samples,
                             int grid_ok) {
    const int T = LONGCALLD_SAMPLE_TUMOR, N = LONGCALLD_SAMPLE_NORMAL;
    int paired = n_samples > 1;
    double lt_abs = smp_ll(acc, T, QG_0);
    double lt_hap = -INFINITY, lt_pool = -INFINITY;
    if (grid_ok)
        for (int j = 0; j < QG_N_PHI; ++j)
            for (int h = 1; h <= 2; ++h) {
                double l = slot_ll(acc, T, h, QG_PHI[j]) + slot_ll(acc, T, 3 - h, QG_0) + slot_ll(acc, T, 0, QG_PHI_HALF[j]);
                lt_hap = lse2(lt_hap, l - log10(2.0 * QG_N_PHI));
            }
    for (int g = 1; g < LONGCALLD_QGRID_N; ++g) lt_pool = lse2(lt_pool, smp_ll(acc, T, g));
    lt_pool -= log10((double)(LONGCALLD_QGRID_N - 1));
    double lt_som = grid_ok ? lse2(lt_hap, lt_pool) - log10(2.0) : lt_pool;
    double h_abs, h_germ, h_som;
    if (paired) {
        double ln_abs = smp_ll(acc, N, QG_0);
        double ln_clean = lse2(ln_abs + log10(0.95), smp_ll(acc, N, QG_CONT) + log10(0.05));
        // The normal's het is read on its haplotypes: a germline allele would
        // fill one of them, and every read of both showing reference is far
        // less likely than the same reads pooled at 1/2. Unlike a germline
        // record's own het (sample_gl), there is no circularity to hedge
        // against -- a somatic candidate never anchors phasing -- so the phased
        // reading leads. It is 0.9/0.1 rather than 1/0, allowing reads the
        // phasing put on the wrong haplotype, and still mixed with the pooled
        // reading for a block the phasing got wrong outright.
        double n_het = smp_ll(acc, N, QG_HALF);
        int n_ph = acc->n_strand[N][0][0] + acc->n_strand[N][0][1] + acc->n_strand[N][1][0] + acc->n_strand[N][1][1]
                 - acc->n_unph[N][0] - acc->n_unph[N][1];
        if (grid_ok && n_ph > 0) {
            double o1 = slot_ll(acc, N, 1, QG_P90) + slot_ll(acc, N, 2, QG_P10) + slot_ll(acc, N, 0, QG_HALF);
            double o2 = slot_ll(acc, N, 1, QG_P10) + slot_ll(acc, N, 2, QG_P90) + slot_ll(acc, N, 0, QG_HALF);
            n_het = lse2(lse2(o1, o2) - log10(2.0) + log10(0.9), n_het + log10(0.1));
        }
        double ln_germ = lse2(n_het + log10(2.0 / 3.0), smp_ll(acc, N, QG_1) + log10(1.0 / 3.0));
        // In the tumour a germline het also sits on one haplotype, so it gets the
        // same tumour model the somatic hypothesis does, plus a lost allele
        // (fraction 0). Scoring the germline tumour on pooled reads alone handed
        // the somatic hypothesis the tumour's haplotype segregation as if it
        // were evidence -- 3.7 of the 5 log10 units of a clonal test call.
        double lt_germ = lse2(lt_som + log10(21.0 / 22.0), lt_abs + log10(1.0 / 22.0));
        h_abs = ln_abs + lt_abs; h_germ = ln_germ + lt_germ; h_som = ln_clean + lt_som;
    } else {
        double gl[3];
        sample_gl(opt, acc, T, grid_ok, 0, gl);
        h_abs = gl[0]; h_germ = lse2(gl[1] + log10(2.0 / 3.0), gl[2] + log10(1.0 / 3.0)); h_som = lt_som;
    }
    double pg = 1.5 * germline_prior(opt, v), pi_s = opt->vq_som_prior;
    double a_abs = log10(1.0 - pg - pi_s) + h_abs, a_g = log10(pg) + h_germ, a_s = log10(pi_s) + h_som;
    double num = lse2(a_abs, a_g);
    return num - lse2(num, a_s);
}

static void add_hist(int *dst, const uint16_t *src) { for (int b = 0; b < LONGCALLD_QHIST_N; ++b) dst[b] += src[b]; }

static void artifact_features(const var_qual_acc_t *acc, const var1_t *v, int n_samples, int is_som, int med_dp,
                              float *x) {
    // Somatic: judged on the tumour, where the alt lives. Germline: on all reads.
    int s0 = 0, s1 = n_samples - 1;
    if (is_som) s0 = s1 = LONGCALLD_SAMPLE_TUMOR;
    int st[2][2] = {{0,0},{0,0}}, unph[2] = {0,0}, conflict = 0, a1 = 0, a2 = 0;
    int bq[2][LONGCALLD_QHIST_N], mq[2][LONGCALLD_QHIST_N];
    memset(bq, 0, sizeof(bq)); memset(mq, 0, sizeof(mq));
    for (int s = s0; s <= s1; ++s) {
        for (int a = 0; a < 2; ++a) {
            st[a][0] += acc->n_strand[s][a][0]; st[a][1] += acc->n_strand[s][a][1];
            unph[a] += acc->n_unph[s][a];
            add_hist(bq[a], acc->bq[s][a]); add_hist(mq[a], acc->mq[s][a]);
        }
        conflict += acc->n_alt_conflict[s];
        a1 += v->slots.alt[s][1]; a2 += v->slots.alt[s][2];
    }
    int n_ref = st[0][0] + st[0][1], n_alt = st[1][0] + st[1][1], dp = n_ref + n_alt;
    for (int f = 0; f < VQF_N; ++f) x[f] = 0;
    x[VQF_BIAS] = 1;
    {
        double rf = st[0][0] + 1, rr = st[0][1] + 1, af = st[1][0] + 1, ar = st[1][1] + 1;
        double r = (rf * ar) / (rr * af);
        double sor = log(r + 1.0 / r) + log((rf < rr ? rf : rr) / (rf > rr ? rf : rr))
                   - log((af < ar ? af : ar) / (af > ar ? af : ar));
        x[VQF_SOR] = (float)sor;
    }
    x[VQF_ONE_STRAND] = (n_alt >= 4 && (st[1][0] == 0 || st[1][1] == 0)) ? 1.0f : 0.0f;
    double mz = rank_sum_z(mq[1], mq[0]), bz = rank_sum_z(bq[1], bq[0]);
    x[VQF_MQ_RS] = (float)(mz < 0 ? -mz : 0);
    if (n_ref > 0 && n_alt > 0) {
        double d = hist_mean(mq[0]) - hist_mean(mq[1]);
        x[VQF_MQ_DIFF] = (float)(d > 0 ? d / 10.0 : 0);
    }
    x[VQF_BQ_RS] = (float)(bz < 0 ? -bz : 0);
    int bmed = hist_median(bq[1]);
    if (bmed >= 0 && bmed < 20) x[VQF_BQ_LOW] = (float)((20 - bmed) / 10.0);
    int is_indel = v->type != BAM_CDIFF, per = v->vq.rep_period, tr = v->vq.rep_len;
    if (is_indel && per == 1) x[VQF_HP_INDEL] = (float)(tr < 20 ? tr : 20);
    if (is_indel && per >= 2) x[VQF_STR_INDEL] = (float)((tr < 40 ? tr : 40) / 2.0);
    if (!is_indel && per == 1 && tr >= 4) x[VQF_HP_SNV] = (float)(tr < 15 ? tr : 15);
    int n_unph = unph[0] + unph[1];
    if (dp > 0) x[VQF_UNPH_FRAC] = (float)n_unph / dp;
    if (n_alt > 0 && n_ref > 0) {
        double ex = (double)unph[1] / n_alt - (double)unph[0] / n_ref;
        x[VQF_UNPH_ALT_EXCESS] = (float)(ex > 0 ? ex : 0);
    } else if (n_alt > 0) x[VQF_UNPH_ALT_EXCESS] = (float)unph[1] / n_alt;
    x[VQF_LOG_UNPH] = (float)log1p(n_unph);
    int is_hom = v->GT[0] != 0 && v->GT[1] != 0;
    if (!is_hom && a1 + a2 > 0) x[VQF_ALT_MINOR_HAP] = (float)(a1 < a2 ? a1 : a2) / (a1 + a2);
    if (n_alt > 0) x[VQF_ALT_CONFLICT] = (float)conflict / n_alt;
    x[VQF_UNPHASED_CALL] = (n_samples > 1 && v->is_unphased_call) ? 1.0f : 0.0f;
    if (med_dp > 0 && v->DP > 0) {
        double r = log2((double)v->DP / med_dp);
        if (r > 0) x[VQF_DP_HIGH] = (float)r; else x[VQF_DP_LOW] = (float)-r;
    }
    x[VQF_LOG_ALT] = (float)log1p(n_alt);
    if (dp > 0) x[VQF_VAF] = (float)n_alt / dp;
}

void score_var_quality(const call_var_opt_t *opt, var_qual_acc_t *acc, int med_dp, var1_t *v) {
    var_qual_t *q = &v->vq;
    fill_ll(acc);
    int n_samples = v->slots.n_samples > 0 ? v->slots.n_samples : 1;
    if (n_samples > 2) n_samples = 2;
    int paired = n_samples > 1;
    int is_som = paired ? (v->paired_class == LONGCALLD_PAIRED_SOMATIC) : v->is_somatic;
    // A grid is usable when the site's reads were counted against one phase
    // block, and -- for a paired verdict -- when the classifier decided on it.
    int grid_ok = v->slots.ps > 0 && !(paired && v->is_unphased_call);
    int is_indel = v->type != BAM_CDIFF;
    q->vq_class = (uint8_t)(is_som ? (is_indel ? VQ_CLASS_SOM_INDEL : VQ_CLASS_SOM_SNV)
                                   : (is_indel ? VQ_CLASS_GERM_INDEL : VQ_CLASS_GERM_SNV));
    for (int s = 0; s < 2; ++s) {
        q->pl[s][0] = q->pl[s][1] = q->pl[s][2] = 0;
        q->n_dp[s] = q->n_unph[s] = q->n_unph_alt[s] = 0;
    }
    for (int s = 0; s < n_samples; ++s) {
        double gl[3];
        // A tumour, or a mosaic call's sample, need not be diploid at the site.
        int any_frac = (paired && s == LONGCALLD_SAMPLE_TUMOR) || (!paired && v->is_somatic);
        sample_gl(opt, acc, s, v->slots.ps > 0, any_frac, gl);
        gl_to_pl(gl, opt->vq_pl_scale, q->pl[s]);
        int dp = acc->n_strand[s][0][0] + acc->n_strand[s][0][1] + acc->n_strand[s][1][0] + acc->n_strand[s][1][1];
        q->n_dp[s] = (int16_t)(dp > INT16_MAX ? INT16_MAX : dp);
        int u = acc->n_unph[s][0] + acc->n_unph[s][1];
        q->n_unph[s] = (int16_t)(u > INT16_MAX ? INT16_MAX : u);
        q->n_unph_alt[s] = (int16_t)(acc->n_unph[s][1] > INT16_MAX ? INT16_MAX : acc->n_unph[s][1]);
    }
    double lfalse = is_som ? somatic_lfalse(opt, acc, v, n_samples, grid_ok)
                           : germline_lfalse(opt, acc, v, n_samples, v->slots.ps > 0);
    double qlik = -10.0 * lfalse;
    if (!(qlik >= 0)) qlik = 0; // NaN and a negative zero both land here
    q->qlik = (float)(qlik > 9999 ? 9999 : qlik);

    artifact_features(acc, v, n_samples, is_som, med_dp, q->feat);
    {
        int st[2][2] = {{0,0},{0,0}};
        int s0 = is_som ? LONGCALLD_SAMPLE_TUMOR : 0, s1 = is_som ? LONGCALLD_SAMPLE_TUMOR : n_samples - 1;
        int mq[2][LONGCALLD_QHIST_N], bq[2][LONGCALLD_QHIST_N];
        memset(mq, 0, sizeof(mq)); memset(bq, 0, sizeof(bq));
        for (int s = s0; s <= s1; ++s)
            for (int a = 0; a < 2; ++a) {
                st[a][0] += acc->n_strand[s][a][0]; st[a][1] += acc->n_strand[s][a][1];
                add_hist(mq[a], acc->mq[s][a]); add_hist(bq[a], acc->bq[s][a]);
            }
        q->sor = q->feat[VQF_SOR];
        q->mq_rs = (float)rank_sum_z(mq[1], mq[0]);
        q->bq_rs = (float)rank_sum_z(bq[1], bq[0]);
    }
    q->feat[VQF_LIK] = (float)((qlik < 60 ? qlik : 60) / 10.0);
    q->feat[VQF_LIK_LOW] = (float)((qlik < 20 ? qlik : 20) / 10.0);
    q->p_art = 0;
    int tech = opt->is_ont ? 1 : 0, cls = q->vq_class;
    const var_qual_weights_t *wt = opt->vq_weights;
    double qual = qlik;
    // The somatic weights were fitted on tumour/normal pairs; a single-sample
    // mosaic call has no normal behind its features, so it keeps QLIK.
    int fitted_here = !(is_som && !paired);
    if (wt != NULL && wt->set[tech][cls] && fitted_here) {
        double z = 0;
        for (int f = 0; f < VQF_N; ++f) z += wt->w[tech][cls][f] * q->feat[f];
        // -10 log10 sigmoid(z), without underflow for very negative z.
        double ln_p = z < 0 ? z - log1p(exp(z)) : -log1p(exp(-z));
        q->p_art = (float)exp(ln_p);
        qual = -10.0 * ln_p / M_LN10;
    }
    if (qual < 0) qual = 0;
    if (qual > opt->max_qual) qual = opt->max_qual;
    v->QUAL = (float)(floor(qual * 10.0 + 0.5) / 10.0); // one decimal, as printed
    v->GQ = var_qual_gq(q->pl[0], (v->GT[0] != 0) + (v->GT[1] != 0));
    q->is_set = 1;
}

// ---------------------------------------------------------------- weights --

static void weights_clear(var_qual_weights_t *wt) { memset(wt, 0, sizeof(*wt)); }

void var_qual_weights_default(var_qual_weights_t *wt) {
    weights_clear(wt);
#include "var_qual_weights.inc"
}

static int lookup_name(const char *s, const char **names, int n) {
    for (int i = 0; i < n; ++i) if (strcmp(s, names[i]) == 0) return i;
    return -1;
}

int var_qual_weights_load(var_qual_weights_t *wt, const char *fn) {
    if (strcmp(fn, "none") == 0) { weights_clear(wt); return 0; }
    FILE *fp = fopen(fn, "r");
    if (fp == NULL) { _err_error("Error: cannot open --qual-weights file %s\n", fn); return -1; }
    char line[512], tech_s[64], cls_s[64], feat_s[64]; double w; int n_line = 0, n_set = 0;
    // A cell named in the file is replaced whole: features it omits become 0,
    // so a file states a model rather than patching one.
    uint8_t cleared[VQ_N_TECH][VQ_N_CLASS]; memset(cleared, 0, sizeof(cleared));
    while (fgets(line, sizeof(line), fp) != NULL) {
        ++n_line;
        if (line[0] == '#' || line[0] == '\n') continue;
        if (sscanf(line, "%63s %63s %63s %lf", tech_s, cls_s, feat_s, &w) != 4) {
            _err_error("Error: %s:%d: expected \"tech class feature weight\"\n", fn, n_line); fclose(fp); return -1;
        }
        int t = strcmp(tech_s, "hifi") == 0 ? 0 : (strcmp(tech_s, "ont") == 0 ? 1 : -1);
        int c = lookup_name(cls_s, LONGCALLD_VQ_CLASS_NAMES, VQ_N_CLASS);
        int f = lookup_name(feat_s, LONGCALLD_VQF_NAMES, VQF_N);
        if (t < 0 || c < 0 || f < 0) {
            _err_error("Error: %s:%d: unknown tech, class or feature in \"%s %s %s\"\n", fn, n_line, tech_s, cls_s, feat_s);
            fclose(fp); return -1;
        }
        if (!cleared[t][c]) { memset(wt->w[t][c], 0, sizeof(wt->w[t][c])); cleared[t][c] = 1; }
        wt->w[t][c][f] = w; wt->set[t][c] = 1; ++n_set;
    }
    fclose(fp);
    return n_set;
}
