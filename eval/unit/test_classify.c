// Table-driven test of the paired verdict (classify_paired_var) on hand-built
// slot grids, at the paired HiFi defaults. Build and run with `make test`.
//
// Each case states the reads in each (sample, slot) as alt/ref counts: slot 0
// is the unphased bucket, slots 1 and 2 the haplotypes of one phase block.
#include <stdio.h>
#include <string.h>
#include "call_var_main.h"
#include "collect_var.h"

// Globals main.c defines for the program; the test links everything else.
const char PROG[20] = "test_classify";
int LONGCALLD_VERBOSE = 0;
char *CMD = "test_classify";

static const char *CLS[] = {"UNSET", "SOMATIC", "GERMLINE", "NORMAL_ONLY", "LOH", "LOWDP", "AMBIG"};
enum { T = LONGCALLD_SAMPLE_TUMOR, N = LONGCALLD_SAMPLE_NORMAL };

typedef struct {
    const char *name;
    int n[3][2], t[3][2];  // [slot][alt, ref]
    int type;              // BAM_CDIFF / BAM_CINS / BAM_CDEL
    int in_repeat, tsd_len, polya_len, te_seq_i;
    int want_cls, want_unphased, want_som_hap;
} case_t;

static int run_case(const call_var_opt_t *opt, const case_t *c) {
    var1_t v; memset(&v, 0, sizeof(v));
    v.type = c->type; v.in_repeat = c->in_repeat;
    v.tsd_len = c->tsd_len; v.polya_len = c->polya_len; v.te_seq_i = c->te_seq_i;
    v.slots.n_samples = 2; v.slots.ps = 100;
    for (int h = 0; h < 3; ++h) {
        v.slots.alt[N][h] = c->n[h][0]; v.slots.ref[N][h] = c->n[h][1];
        v.slots.alt[T][h] = c->t[h][0]; v.slots.ref[T][h] = c->t[h][1];
    }
    classify_paired_var(opt, &v);
    int ok = v.paired_class == c->want_cls && v.is_unphased_call == c->want_unphased &&
             (c->want_som_hap < 0 || v.som_hap == c->want_som_hap);
    printf("  %s  %-48s %s%s som_hap=%d\n", ok ? "PASS" : "FAIL", c->name, CLS[v.paired_class],
           v.is_unphased_call ? " (unphased)" : "", v.som_hap);
    return ok;
}

int main(void) {
    call_var_opt_t *opt = call_var_init_para();
    // The paired HiFi presets that finalise_sample_layout() applies.
    opt->min_af = LONGCALLD_PAIRED_MIN_CAND_AF;
    opt->som_pval = LONGCALLD_PAIRED_SOM_PVAL;
    opt->min_hap_vaf = LONGCALLD_PAIRED_MIN_HAP_VAF;
    //  name                                   normal {s0,h1,h2}            tumour {s0,h1,h2}           type       rep tsd pA  te   want                       unph hap
    const case_t cases[] = {
        {"normal too shallow -> LowDP",          {{0,0},{0,3},{0,3}},     {{0,0},{5,5},{0,10}},    BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_LOWDP,       0, -1},
        {"alt on one hap in both -> germline",   {{0,0},{10,0},{0,10}},   {{0,0},{10,0},{0,10}},   BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_GERMLINE,    0, -1},
        {"tumour hap 1 only -> somatic",         {{0,0},{0,10},{0,10}},   {{0,0},{6,4},{0,10}},    BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_SOMATIC,     0, 1},
        {"stray alt within 20% of other hap",    {{0,0},{0,10},{0,10}},   {{0,0},{0,10},{6,4}},    BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_SOMATIC,     0, 2},
        {"alt scattered on both tumour haps",    {{0,0},{0,10},{0,10}},   {{0,0},{3,7},{3,7}},     BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_AMBIG,       0, -1},
        {"normal-only, tumour depleted 4:1 on hap", {{0,0},{10,0},{0,10}}, {{0,0},{0,4},{0,20}},   BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_LOH,         0, -1},
        // Known gap, pinned rather than endorsed: with no tumour reads at all
        // on the allele's haplotype the grid branch never runs (it requires
        // tumour depth on that haplotype), so its "tumour lacks the haplotype"
        // LOH test cannot fire; the pooled fallback keeps a phase-set anchor
        // GERMLINE.
        {"normal-only, tumour lost the hap entirely", {{0,0},{10,0},{0,10}}, {{0,0},{0,0},{0,15}},  BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_GERMLINE,    1, -1},
        {"normal-only, tumour kept the hap",     {{0,0},{10,0},{0,10}},   {{0,0},{0,10},{0,10}},   BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_NORMAL_ONLY, 0, -1},
        {"SNV below the tumour VAF floor",       {{0,0},{0,10},{0,10}},   {{0,0},{2,8},{0,40}},    BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_AMBIG,       0, -1},
        {"repeat indel below its 0.30 floor",    {{0,0},{0,10},{0,10}},   {{0,0},{8,12},{0,20}},   BAM_CINS,  1, 0, 0, -1, LONGCALLD_PAIRED_AMBIG,       0, -1},
        {"same indel outside a repeat",          {{0,0},{0,10},{0,10}},   {{0,0},{8,12},{0,20}},   BAM_CINS,  0, 0, 0, -1, LONGCALLD_PAIRED_SOMATIC,     0, 1},
        {"one-read MEI with TSD and poly-A",     {{0,0},{0,10},{0,10}},   {{0,0},{1,9},{0,10}},    BAM_CINS,  0, 9, 25, 0, LONGCALLD_PAIRED_SOMATIC,     0, 1},
        {"one-read insertion, no TE signature",  {{0,0},{0,10},{0,10}},   {{0,0},{1,9},{0,10}},    BAM_CINS,  0, 0, 0, -1, LONGCALLD_PAIRED_AMBIG,       0, -1},
        {"pooled: clean normal -> somatic",      {{0,10},{0,0},{0,0}},    {{6,14},{0,0},{0,0}},    BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_SOMATIC,     1, -1},
        {"pooled: both carry it -> germline",    {{5,5},{0,0},{0,0}},     {{5,5},{0,0},{0,0}},     BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_GERMLINE,    1, -1},
        {"pooled: normal het, tumour lacks it",  {{5,5},{0,0},{0,0}},     {{0,20},{0,0},{0,0}},    BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_NORMAL_ONLY, 1, -1},
        {"pooled: normal noise -> ambiguous",    {{2,58},{0,0},{0,0}},    {{0,20},{0,0},{0,0}},    BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_AMBIG,       1, -1},
        {"pooled: phase-set anchor stays germline", {{0,0},{8,0},{0,8}},  {{0,20},{0,0},{0,0}},    BAM_CDIFF, 0, 0, 0, -1, LONGCALLD_PAIRED_GERMLINE,    1, -1},
    };
    int n = sizeof(cases) / sizeof(cases[0]), n_ok = 0;
    for (int i = 0; i < n; ++i) n_ok += run_case(opt, cases + i);
    printf("\n  classify_paired_var: %d/%d cases passed\n", n_ok, n);
    call_var_free_para(opt);
    return n_ok == n ? 0 : 1;
}
