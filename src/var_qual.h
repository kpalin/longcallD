#ifndef LONGCALLD_VAR_QUAL_H
#define LONGCALLD_VAR_QUAL_H

#include <stdint.h>
#include "htslib/sam.h"

// Variant (QUAL) and genotype (GQ/PL) quality.
//
// Two error sources are modelled separately, because they behave differently
// with depth:
//
//  1. Sampling error: the alt reads are sequencing or mapping errors drawn
//     independently read by read. Each read carries its own error rate, built
//     from its base quality, its mapping quality and, for an indel, the repeat
//     context the indel sits in. That gives genotype and somatic likelihoods,
//     and QLIK, which grows without bound with depth.
//  2. Systematic error: the alt reads agree because the error is shared --
//     strand-specific chemistry, a paralogous mapping, slippage in a long
//     repeat, a read population the phasing could not place. More reads do not
//     dilute it, so it cannot be a likelihood over independent reads. It is a
//     logistic model over features of the read population, fitted per
//     technology against truth sets by eval/qual_fit.py.
//
// The record is false if either explains it:
//     P(false) = 1 - (1 - P_lik(false)) * (1 - P_art)
// and QUAL = -10 log10 P(false). Where no artifact model is set, QUAL == QLIK.

#ifdef __cplusplus
extern "C" {
#endif

struct call_var_opt_t;
struct bam_chunk_t;
struct var1_t;

// Allele-fraction grid the per-read likelihoods are accumulated on. It holds
// every fraction the models evaluate -- 0, 1/2, 1, the somatic haplotype
// fractions and their halves -- so pooling slots or samples is an exact sum.
#define LONGCALLD_QGRID_N 22
extern const double LONGCALLD_QGRID[LONGCALLD_QGRID_N];

// Histogram bins for base and mapping quality; larger values land in the last
// bin. Rank-sum statistics are computed from these.
#define LONGCALLD_QHIST_N 64

// Profile cells whose base quality is not known: noisy-region alleles come from
// a re-alignment, not from a base in the read.
#define LONGCALLD_QUAL_UNKNOWN 0xFF

// Per-read error rates are binned to whole phred units: reads in one bin
// contribute identically, so the likelihoods are summed over occupied bins at
// scoring time rather than evaluated per read -- 22 log10 calls per read and
// variant were 14% of the run time.
#define LONGCALLD_QE_N 91 // phred 0..90

// Per-candidate accumulators, chunk-scoped. Built in make_variants() once read
// haplotypes are final, i.e. after stitching.
typedef struct var_qual_acc_t {
    uint16_t ne[2][3][2][LONGCALLD_QE_N]; // [sample][slot][ref/alt][error phred]: reads
    float ll[2][3][LONGCALLD_QGRID_N];    // [sample][slot][f]: sum over reads of log10 P(read | alt fraction f); filled by score_var_quality()
    int n_strand[2][2][2];                // [sample][ref/alt][fwd/rev]
    int n_unph[2][2];                     // [sample][ref/alt]: reads in slot 0 -- unassigned, or outside the site's block
    int n_alt_conflict[2];                // [sample]: alt reads whose own phasing has >= 2 conflicting SNPs
    uint16_t bq[2][2][LONGCALLD_QHIST_N]; // [sample][ref/alt][base quality]
    uint16_t mq[2][2][LONGCALLD_QHIST_N]; // [sample][ref/alt][mapping quality]
} var_qual_acc_t;

// Artifact-model features. The order is the weight-file order.
enum {
    VQF_BIAS = 0,        // intercept
    VQF_SOR,             // GATK StrandOddsRatio, over the alt-carrying sample
    VQF_ONE_STRAND,      // every alt read on one strand, given >= 4 alt reads
    VQF_MQ_RS,           // max(0, -z): alt mapping quality ranks below ref
    VQF_MQ_DIFF,         // max(0, mean ref MQ - mean alt MQ) / 10
    VQF_BQ_RS,           // max(0, -z): alt base quality ranks below ref
    VQF_BQ_LOW,          // max(0, 20 - median alt BQ) / 10
    VQF_HP_INDEL,        // indel in a homopolymer: min(run, 20)
    VQF_STR_INDEL,       // indel in a period 2-6 tandem repeat: min(tract bases, 40) / 2
    VQF_HP_SNV,          // substitution in or beside a homopolymer of >= 4: min(run, 15)
    VQF_UNPH_FRAC,       // unphased share of the alt-carrying sample's reads
    VQF_UNPH_ALT_EXCESS, // unphased share of its alt reads minus that of its ref reads, floored at 0
    VQF_LOG_UNPH,        // log1p(its unphased reads)
    VQF_ALT_MINOR_HAP,   // share of its phased alt reads on the haplotype holding fewer of them
    VQF_ALT_CONFLICT,    // share of its alt reads whose own phasing conflicts
    VQF_UNPHASED_CALL,   // paired: the verdict was reached on pooled counts
    VQF_DP_HIGH,         // max(0, log2(DP / chunk median DP))
    VQF_DP_LOW,          // max(0, -log2(DP / chunk median DP))
    VQF_LOG_ALT,         // log1p(alt reads in the alt-carrying sample)
    VQF_VAF,             // alt fraction in the alt-carrying sample
    VQF_N
};
extern const char *LONGCALLD_VQF_NAMES[VQF_N];

// Weight tables: [tech][class][feature], tech 0 = HiFi, 1 = ONT.
#define VQ_CLASS_GERM_SNV   0
#define VQ_CLASS_GERM_INDEL 1
#define VQ_CLASS_SOM_SNV    2
#define VQ_CLASS_SOM_INDEL  3
#define VQ_N_CLASS 4
#define VQ_N_TECH  2
extern const char *LONGCALLD_VQ_CLASS_NAMES[VQ_N_CLASS];
typedef struct var_qual_weights_t {
    double w[VQ_N_TECH][VQ_N_CLASS][VQF_N];
    uint8_t set[VQ_N_TECH][VQ_N_CLASS]; // 0: no artifact model for this cell, so QUAL == QLIK
} var_qual_weights_t;

// What a record carries into the VCF.
typedef struct var_qual_t {
    float qlik;              // likelihood-only phred
    float p_art;             // artifact-model probability
    float sor, mq_rs, bq_rs; // strand odds ratio; rank-sum z of alt against ref (negative: alt lower)
    float feat[VQF_N];       // artifact features, printed under LONGCALLD_QUALCHK
    int16_t pl[2][3];        // [sample] normalised PL for 0/0, 0/1, 1/1, capped at 9999
    int16_t n_dp[2], n_unph[2], n_unph_alt[2]; // [sample]
    uint8_t rep_period;      // 0 none, 1 homopolymer, 2-6 tandem repeat
    uint8_t rep_len;         // tract length in reference bases, capped at 255
    char rep_unit[7];        // repeat unit, NUL-terminated
    uint8_t vq_class;        // VQ_CLASS_*
    uint8_t is_set;          // scored; 0 leaves QUAL/GQ to the caller
} var_qual_t;

// Step 1, while digars exist: record each profile cell's base quality.
void fill_read_var_profile_qual(const struct call_var_opt_t *opt, struct bam_chunk_t *chunk);
// Step 2: accumulate per-candidate evidence, on the site's phase block as
// chunk->var_slots[].ps names it. Returns n_cand_vars accumulators, or NULL.
var_qual_acc_t *collect_var_qual_acc(const struct call_var_opt_t *opt, struct bam_chunk_t *chunk,
                                     const int *want);
// Repeat context at a reference position: an indel's event base, a substitution's base.
void ref_repeat_context(const char *ref_seq, hts_pos_t beg, hts_pos_t end, hts_pos_t pos, int is_snv,
                        uint8_t *period, uint8_t *tract_len, char *unit);
// Fill v->vq and set v->QUAL and v->GQ. med_dp is the chunk's median output depth.
void score_var_quality(const struct call_var_opt_t *opt, var_qual_acc_t *acc, int med_dp,
                       struct var1_t *v);
// GQ of a printed genotype carrying `dosage` alt alleles, from that sample's PL.
int var_qual_gq(const int16_t pl[3], int dosage);
void var_qual_weights_default(var_qual_weights_t *wt);
// Overlay weights from a file of "tech class feature weight" lines.
int var_qual_weights_load(var_qual_weights_t *wt, const char *fn);

#ifdef __cplusplus
}
#endif

#endif
