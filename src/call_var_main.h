#ifndef LONGCALLD_CALL_VAR_MAIN_H
#define LONGCALLD_CALL_VAR_MAIN_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "seq.h"
#include "kmer.h"
#include "cgranges.h"
#include "wavefront/wavefront_align.h"
#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"
#include "bam_utils.h"

#define CALL_VAR_PL_THREAD_N 2 // num of threads for pipeline
#define CALL_VAR_THREAD_N 8    // num of threads for variant calling

#define LONGCALLD_MIN_CAND_MQ 30 // ignore reads with MAPQ < 30
#define LONGCALLD_MIN_CAND_BQ 10 // ignore bases with BQ < 10
#define LONGCALLD_MIN_CAND_DP 5 // total depth < 5: skipped
#define LONGCALLD_MIN_ALT_DP 2 // alt depth < 2: skipped
#define LONGCALLD_MIN_CAND_AF 0.20 // AF < 0.20: not germline het.
#define LONGCALLD_MAX_CAND_AF 0.80 // AF > 0.80: not germline het.
// Paired-mode candidate and somatic-test defaults, tuned on SEQC2 HCC1395 chr1
// (PacBio HiFi) with som.py. They apply ONLY when -N is given, so single-sample
// calling is untouched; an explicit -a or --som-pval still wins.
#define LONGCALLD_PAIRED_MIN_CAND_AF 0.05
#define LONGCALLD_PAIRED_SOM_PVAL 1.0
// Within-slot alt fraction needed to assert an ALT allele that the pooled
// haplotype-consensus vote did not find. See make_variants(). Tuned on
// SEQC2 HCC1395 chr1: the gain is flat from 0.30 to 0.40 and collapses by 0.50.
#define LONGCALLD_PAIRED_MIN_HAP_VAF 0.35
// Paired-mode classifier thresholds. All of these are read only by
// classify_paired_var(), which runs only when a matched normal is given, so
// they need no single-sample default and cannot affect single-sample output.
// Tuned together on SEQC2 HCC1395 chr1; see the comments at each use.
#define LONGCALLD_MAX_UNPHASED_FRAC 0.02      // above this, the haplotype grid is not trusted
#define LONGCALLD_MAX_HAP_ALT_FRAC 0.20       // alt tolerated in a slot that must be clean
#define LONGCALLD_MIN_SOM_VAF 0.05            // tumour VAF floor, substitutions
#define LONGCALLD_MIN_SOM_INDEL_VAF 0.15      // tumour VAF floor, indels
#define LONGCALLD_MIN_SOM_REP_INDEL_VAF 0.30  // tumour VAF floor, indels in a repeat

// Artifact checks a paired somatic call is put through. Each bit is one check;
// --som-artifact-mask selects which are enforced, and LONGCALLD_ARTCHK=1
// reports every call's full mask so the set can be re-measured.
#define LONGCALLD_ART_STRAND      0x001 // alt seen on one strand only (ONT only)
#define LONGCALLD_ART_DIS_TO_VAR  0x002 // another variant within min_somatic_dis_to_var
#define LONGCALLD_ART_ALT_QUAL    0x004 // median alt base quality below the chunk median
#define LONGCALLD_ART_WIN_QUAL    0x008 // low-quality window around the alt base
#define LONGCALLD_ART_INDEL_ERR   0x010 // an indel error within min_somatic_dis_to_seq_error
#define LONGCALLD_ART_DENSE_DIFF  0x020 // alt reads carry dense mismatches
#define LONGCALLD_ART_LONG_CLIP   0x040 // a long soft clip nearby
#define LONGCALLD_ART_LOW_COMP    0x080 // error inside a low-complexity region
#define LONGCALLD_ART_HOMOPOLYMER 0x100 // homopolymer error
#define LONGCALLD_ART_ALL         0x1ff
// Enforced by default: the six checks that cost essentially no true positives
// on PacBio HiFi. Measured over 7260 paired somatic calls on chr1 (SEQC2
// HCC1395, som.py inside the HC BED; the instrumented run lost chr1:137.5-150Mb,
// so the counts below cover the rest), as the TP and FP each check would remove
// if enforced alone:
//
//   STRAND         0 TP /  0 FP   (ONT only, inert here)
//   ALT_QUAL       1 TP /  2 FP
//   WIN_QUAL       0 TP /  1 FP
//   INDEL_ERR      1 TP /  2 FP
//   DENSE_DIFF     2 TP /  2 FP
//   LONG_CLIP      0 TP /  1 FP
//   DIS_TO_VAR    67 TP /  8 FP   <- excluded
//   LOW_COMP     120 TP /  4 FP   <- excluded
//   HOMOPOLYMER 1281 TP / 17 FP   <- excluded
//
// The three excluded checks come from the single-sample mosaic path, where a
// candidate is a lone low-fraction allele and clustering with a repeat is good
// evidence of slippage. A paired somatic call already had to clear the
// four-slot test, so the same clustering mostly marks real variants sitting in
// repeats: enforcing HOMOPOLYMER alone takes SNV F1 from 0.9171 to 0.6607.
// The default set takes SNV F1 from 0.9171 to 0.9208 and indel F1 from 0.4978
// to 0.5045 on full chr1, removing 31 false positives for 4 true ones. Most of
// what it catches is a cluster of pooled-branch calls in the repeat-rich 1q21
// region, which is the shape DENSE_DIFF and LONG_CLIP exist to recognise.
#define LONGCALLD_ART_DEFAULT (LONGCALLD_ART_STRAND | LONGCALLD_ART_ALT_QUAL | \
                               LONGCALLD_ART_WIN_QUAL | LONGCALLD_ART_INDEL_ERR | \
                               LONGCALLD_ART_DENSE_DIFF | LONGCALLD_ART_LONG_CLIP)
#define LONGCALLD_DEF_PLOID 2 // diploid
#define LONGCALLD_REF_ALLELE 0
#define LONGCALLD_ALT_ALLELE1 1
#define LONGCALLD_ALT_ALLELE2 2

// define noisy region parameters:
#define LONGCALLD_NOISY_REG_MAX_XGAPS 5 //5 // or 10; dense X/gap region: more than n X/gap bases in a 100-bp window
#define LONGCALLD_NOISY_REG_HIFI_SLIDE_WIN 100 //200
#define LONGCALLD_NOISY_REG_ONT_SLIDE_WIN 25
#define LONGCALLD_MAX_NOISY_FRAC_PER_READ 0.5 // skip reads with more than 50% bases in noisy region
#define LONGCALLD_MAX_VAR_RATIO_PER_READ 0.05 // skip reads with n_var / ref_span > 5% 
#define LONGCALLD_MAX_READ_DEPTH 500 // vars with >500 reads will be skipped
#define LONGCALLD_MAX_NOISY_REG_COV 1000 // regions with >1000 reads will be skipped
#define LONGCALLD_NOISY_END_CLIP 30 // 100 // >= 100 bp clipping on both ends will be considered as long clipping
#define LONGCALLD_NOISY_END_CLIP_WIN 100 // 100 bp next to the long end-clipping will be considered as low-quality region
#define LONGCALLD_NOISY_REG_MERGE_DIS 500 // 500 bp, merge noisy regions within 500 bp
#define LONGCALLD_NOISY_REG_FLANK_LEN 10 // during re-alignment, include 10-bp flanking region for both ends of noisy region

#define LONGCALLD_MAX_NOISY_REG_LEN 50000 // >50kb noisy region will be skipped
#define LONGCALLD_NOISY_REG_READS 2 // >= 5 reads supporting noisy region
// #define LONGCALLD_NOISY_REG_RATIO 0.20 // >= 25% reads supporting noisy region

#define LONGCALLD_MIN_NOISY_REG_SIZE_TO_SAMPLE_READS 10000 // >=10kb noisy region, sample reads for re-alignment
#define LONGCALLD_PARTIAL_ALN_RATIO 1.1 // max length ratio for partial alignment, i.e. longer_aln_len / shorter_aln_len <= 1.1

#define LONGCALLD_MIN_HAP_FULL_READS 1 // >= 1full read supporting each haplotype
#define LONGCALLD_MIN_HAP_READS 2 // >= 2 reads supporting each haplotype, including partial/clipped reads, call consensus from >= 3 reads
// #define LONGCALLD_MIN_NO_HAP_FULL_READS 10 // >10 total full reads in noisy region
#define LONGCALLD_MIN_READ_TO_HAP_CONS_SIM 0.9 // for reads with >= 90% equal bases, assign haplotype
// SV related parameters
#define LONGCALLD_MIN_SV_LEN 30 // size >= 30 bp -> SV
#define LONGCALLD_MIN_TSD_LEN 2 // TSD >= 2 bp
#define LONGCALLD_MAX_TSD_LEN 100 // TSD <= 100 bp
#define LONGCALLD_MIN_POLYA_LEN 10 // polyA >= 10 bp
#define LONGCALLD_MIN_POLYA_RATIO 0.8 // polyA >= 80% of the total length
// SOMATIC/MOSAIC Var
#define LONGCALLD_MIN_SOMATIC_DIS_TO_VAR 5
#define LONGCALLD_MIN_SOMATIC_DIS_TO_HP_INDEL_ERROR 3
#define LONGCALLD_MIN_SOMATIC_DIS_TO_INDEL_VAR 5
#define LONGCALLD_MIN_SOMATIC_DIS_TO_SEQ_ERROR 5
// XXX MOSAIC (low) vs SOMATIC (high)
// #define LONGCALLD_MIN_SOMATIC_AF 0.01
// #define LONGCALLD_MAX_SOMATIC_ALT_AF 0.1
#define LONGCALLD_SOMATIC_WIN 1000
#define LONGCALLD_SOMATIC_WIN_MAX_VARS 5 // >= 2 somatic vars in a 1000-bp window, consider as artifact, tag both var & read as artifact
#define LONGCALLD_MIN_SOMATIC_HAP_READS 5 // >= 10 reads supporting each haplotype
#define LONGCALLD_MIN_SOMATIC_ALT_DP 2 // >= 2 reads supporting somatic variant
#define LONGCALLD_MIN_SOMATIC_TE_ALT_DP 1 // >= 1 read supporting somatic TE variant
// fisher is NOT used for PacBio-HiFi
#define LONGCALLD_MIN_SOMATIC_FISHER_PVAL 0.05
#define LONGCALLD_STRAND_BIAS_PVAL_ONT 0.01
// BLT50: 1/10/3 works better
// #define LONGCALLD_SOMATIC_BETA_ALPHA 2 // beta prior for somatic variant calling
// #define LONGCALLD_SOMATIC_BETA_BETA 10 // beta prior for somatic variant calling
// #define LONGCALLD_MIN_SOMATIC_LOG_BETA_BINOM -4 // min log10(p-value) for beta-binomial test


// for sdust
#define LONGCALLD_SDUST_T 5  // 10
#define LONGCALLD_SDUST_W 20 // 50

// for math_utils
#define LONGCALLD_LGAMMA_MAX_I 500




#ifdef __cplusplus
extern "C" {
#endif

// for each vartiant:
// #site-level:
// * QUAL = -10 * log(1-p) (int), site-wise, p = P(data|ref)
// * FILTER: PASS, LowQual, RefCall, NoCall
// * INFO: DP/AD/AF/END
// #sample-level:
// * FORMAT: GT/GQ/DP/AD/PL/PS 
//     GP = P(G|D)
//     GL = P(D|G)
//   * PL (int) = -10 * log10(P(D|G)), Phred-scaled likelihoods of the possible genotypes.
//        "Normalized": the PL of the most likely genotype is set to 0, and the rest are scaled relative to this.
//   * GQ (int) = PL_sec - PL_lowest (int)
//   * GT = argmin(PL)
typedef struct {
    uint8_t type; // BAM_CINS/BAM_CDEL/BAM_CDIFF
    hts_pos_t pos, PS; // phase set
    uint8_t *ref_bases; int ref_len;
    uint8_t is_somatic, is_sv, is_clean; // clean: SNP or simple indel in non-repetitive region
    uint8_t *tsd_seq; int tsd_len, polya_len; hts_pos_t tsd_pos1, tsd_pos2; // target site duplication, 2 TSDs for DEL
    int te_seq_i, te_is_rev; // char *rep_name, *rep_family, *rep_class; use te_seq_i to retrieve TE sequence info

    int n_alt_allele; // currently, only one alt allele
    uint8_t **alt_bases; int *alt_len;
    int DP, AD[2]; uint8_t GT[2]; // DP/AD/GT
    int QUAL, GQ, PL[6]; // phred-scaled, QUAL/FILTER/GQ/PL
    int *alt_read_i;

    // Paired tumour/normal evidence and verdict. Unset in single-sample mode.
    var_slot_counts_t slots;  // the per-(sample, haplotype) grid this was decided on
    uint8_t paired_class;     // LONGCALLD_PAIRED_*
    uint8_t in_repeat;        // indel sits in a homopolymer or short tandem repeat
    uint8_t som_hap;          // the tumour haplotype carrying a somatic alt; 0 if unphased
    uint8_t is_unphased_call; // decided on pooled per-sample counts, not on the slot grid
    float hap_vaf;            // alt fraction within that haplotype; estimates purity x cancer-cell fraction
    double paired_pval;       // Fisher exact, one slot against the other three
} var1_t;

typedef struct var_t {
    int n, m;
    var1_t *vars;
} var_t;

typedef struct call_var_opt_t {
    // input
    char *ref_fa_fn, *ref_fa_fai_fn;
    char *reg_bed_fn;
    char *sample_name;
    int n_in_bam_fn, m_in_bam_fn, input_is_list; char **in_bam_fns;
    // Paired tumour/normal input. The -N files are appended to in_bam_fns only
    // after set_bam_list_file() has finished reordering it, so sample identity
    // is simply "index >= n_tumor_bam_fn" and cannot be scrambled by that
    // reordering.
    int n_normal_bam_fn, m_normal_bam_fn; char **normal_bam_fns;
    int n_tumor_bam_fn;        // in_bam_fns[0 .. n_tumor_bam_fn-1] are the tumour/primary sample
    int n_samples;             // 1 normally, 2 once a normal is supplied
    int *bam_i_to_sample;      // size n_in_bam_fn
    char **sample_names;       // size n_samples
    char *tumor_name, *normal_name; // explicit VCF column names, if given
    uint8_t paired_mode;
    samFile **bams; bam_hdr_t **headers; hts_idx_t **idxs; // for write to output bam/cram
    uint8_t is_pb_hifi, is_ont; float strand_bias_pval; // for ONT reads
    // variant calling regions
    uint8_t only_autosome, only_autosome_XY;
    char **exc_tnames; int n_exc_tnames;
    // char *region_list; uint8_t region_is_file; // for -R/--region/--region-file option
    // filters for variant calling
    int max_ploid, min_mq, min_bq, min_dp, min_alt_dp;
    double min_af, max_af;
    uint8_t min_af_set, som_pval_set, min_hap_vaf_set; // user passed -a / --som-pval / --min-hap-vaf: do not override in paired mode
    double min_hap_vaf;               // within-slot alt fraction needed to assert an ALT the pooled vote missed
    double max_unphased_frac;         // share of unphased reads above which the haplotype grid is not trusted
    double max_hap_alt_frac;          // alt tolerated in a slot that must be clean, as a fraction of its depth
    double min_som_vaf;               // tumour VAF floor for a somatic substitution
    double min_som_indel_vaf;         // tumour VAF floor for a somatic indel
    double min_som_rep_indel_vaf;     // tumour VAF floor for a somatic indel in a homopolymer or short tandem repeat
    int som_artifact_mask;            // LONGCALLD_ART_* bits enforced on paired somatic calls
    int min_normal_phase_dp; // paired mode: min normal depth for a variant to drive phasing
    // Paired tumour/normal classification thresholds.
    int min_normal_dp, min_tumor_dp;         // per-sample depth floor before we say anything
    int min_normal_hap_dp, min_tumor_hap_dp; // per-haplotype depth floor for the four-slot test
    int max_hap_alt;                         // alt reads tolerated in a slot that should be clean
    double som_pval;                         // Fisher exact, one slot against the other three
    double min_normal_phase_af, max_normal_phase_af; // AF window a normal site must sit in to anchor phasing
    // somatic/mosaic variant
    int min_somatic_dis_to_var, min_somatic_dis_to_homopolymer_indel_error, min_somatic_dis_to_seq_error, min_somatic_alt_dp, min_somatic_te_dp, min_somatic_hap_dp;
    double min_somatic_log_beta_binom, min_somatic_fisher_pval; //, max_somatic_alt_af; // 0.01~0.1
    int somatic_beta_alpha, somatic_beta_beta; // beta prior for somatic variant calling
    int somatic_win, somatic_win_max_vars; // somatic variant calling window size, max vars in the window

    int noisy_reg_max_xgaps, noisy_reg_slide_win;
    int end_clip_reg, end_clip_reg_flank_win;
    int noisy_reg_merge_dis, noisy_reg_flank_len; // noisy_reg_merge_win; // for re-alignment
    // filters for noisy region, i.e., coverage/ratio
    int max_noisy_reg_len, max_noisy_reg_cov; //, min_noisy_reg_reads; 
    int64_t max_poa_cells; // largest abPOA DP matrix to attempt, in cells
    int64_t max_poa_msa_cells; // same, for the unbanded whole-MSA call site
    double max_var_ratio_per_read, max_noisy_frac_per_read; //, min_noisy_reg_ratio;
    int min_hap_full_reads, min_hap_reads; //, min_no_hap_full_reads;
    // alignment
    int match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2;
    int gap_aln; // default: 1: left (minimap2, abpoa), 2: right (wfa2)
    double min_read_to_hap_cons_sim, partial_aln_ratio;
    int min_noisy_reg_size_to_sample_reads; //, sampling_hap_read_count, sampling_non_hap_read_count;
    // int disable_read_sampling; // disable read-sampling in long noisy regions; by default read sampling is enabled; disable to capture mosaic variants in long noisy regions
    // TSD & polyA/T
    int min_tsd_len, max_tsd_len, min_polya_len; float min_polya_ratio;
    // Alu/L1/SVA sequences
    char *te_seq_fn; char **te_seq_names; int n_te_seqs;
    int te_kmer_len; kmer32_hash_t **te_for_h, **te_rev_h;
    int output_var_rnames, output_sv_rnames, output_somatic_var_rnames; // output supporting read IDs
    // general
    // int max_ploidy;
    int pl_threads, n_threads;
    // math utils
    double lgamma_cache[LONGCALLD_LGAMMA_MAX_I+1]; int min_lgamma_i, max_lgamma_i; // 0, 999

    // output
    int min_sv_len; // classify as SV if length >= min_sv_len (50)
    htsFile *out_aln_fp; uint8_t out_aln_is_cram; uint8_t refine_bam; // phased bam (tumour / single sample)
    htsFile *out_naln_fp; uint8_t out_naln_is_cram;                   // matched normal's phased bam
    // One output per sample, indexed by LONGCALLD_SAMPLE_*. NULL means that
    // sample's reads are not written; read_i must still advance over them.
    struct aln_out_t *aln_outs[LONGCALLD_MAX_SAMPLES];
    uint8_t no_reorder;        // write records in input order and declare SO:unsorted
    uint8_t out_aln_passthrough; // also emit reads on excluded contigs and unmapped reads
    htsFile *out_vcf; bcf_hdr_t *vcf_hdr; char *out_vcf_fn; char out_vcf_type; // u/b/v/z
    double p_error, log_p, log_1p, log_2; int max_gq; int max_qual;
    int8_t no_vcf_header, out_amb_base, out_somatic, out_methylation;
} call_var_opt_t;

typedef struct {
    int n_regions, m_regions;
    int *reg_tids; hts_pos_t *reg_begs, *reg_ends;
} reg_chunks_t;

typedef struct call_var_io_aux_t {
    faidx_t *fai;
    int n_bam;
    samFile **bams; bam_hdr_t **headers; hts_idx_t **idxs;
} call_var_io_aux_t; // per thread

// shared data for all threads
typedef struct call_var_pl_t {
    // input files
    call_var_io_aux_t *io_aux;
    // TODO mv aligner_obj (WFA/abPOA) here
    // parameters, output files
    struct call_var_opt_t *opt;
    // m-threads
    // struct bam_chunk_t *last_chunk; int n_last_chunk_reads, *last_chunk_read_i; hts_pos_t cur_active_reg_beg;
    // int max_reads_per_chunk, 
    int min_reg_chunks_per_run, max_reg_len_per_chunk;
    int reg_chunk_i, n_reg_chunks, m_reg_chunks; reg_chunks_t *reg_chunks;
    int n_threads;
    // Per-sample read tallies, summed in step 2 (which is serial, so no atomics).
    int64_t n_reads_per_sample[LONGCALLD_MAX_SAMPLES];
} call_var_pl_t;

struct bam_chunk_t;
// separately/parallelly processed
typedef struct call_var_step_t {
    struct call_var_pl_t *pl;
    int n_chunks, max_chunks;
    struct bam_chunk_t *chunks; // input, size: n_chunks
    var_t *vars; // output, size: n_chunks
} call_var_step_t;

int call_var_main(int argc, char *argv[]);
int paired_var_main(int argc, char *argv[]);
int call_var_main_mode(int argc, char *argv[], int paired_cli);
void var_free(var_t *v);
void var1_free(var1_t *v);

#ifdef __cplusplus
}
#endif

#endif // end of LONGCALLD_call_var_H