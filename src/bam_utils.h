#ifndef LONGCALLD_BAM_UTILS_H
#define LONGCALLD_BAM_UTILS_H

#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "cgranges.h"
#include "collect_var.h"

#define LONGCALLD_BAM_CHUNK_READ_COUNT 4000
#define LONGCALLD_BAM_CHUNK_REG_SIZE 500000 // 0.5M/1M
// region-based BAM_CHUNK: no split for a single chromosome
//                         multi regions may be merged into a single chunk if n_regs < 64
// #define LONGCALLD_BAM_MIN_REG_CHUNK_PER_RUN 64 // n_threads * 4


#define BAM_RECORD_LOW_QUAL 0x1
#define BAM_RECORD_WRONG_MAP 0x2
// #define BAM_RECORD_LARGE_CLIP 0x4

#define bam_bseq2base(bseq, qi) seq_nt16_str[bam_seqi(bseq, qi)]

#ifdef __cplusplus
extern "C" {
#endif

// read/base wise X/I/D operations from CIGAR
typedef struct digar1_t {
    // pos: 1-based ref position
    // qi: 0-based query position (for =XI); for DEL, qi is the first read base after the deletion
    hts_pos_t pos; int type, len, qi;
    uint8_t *alt_seq; // X/I: alt_seq; =/D: NULL
    uint8_t is_low_qual; // low-base-qual, effective for clean-region digars
} digar1_t;

typedef struct {
    hts_pos_t beg, end; // [beg, end]
    uint8_t is_rev;
    int n_digar, m_digar;
    digar1_t *digars;
    // read-wise noisy region: active region for re-alignment
    cgranges_t *noisy_regs; // merge low_qual digar1_t if they are next to each other
    int qlen; uint8_t *bseq, *qual; // copy from bam1_t
} digar_t; // detailed CIGAR for each read

// Per-sample, per-haplotype allele tallies for one candidate variant.
//
// With phasing shared across samples, haplotype index h means the same parental
// chromosome in both, so these form a 2x2 grid of comparable slots --
// (normal,H1), (normal,H2), (tumour,H1), (tumour,H2) -- plus an unphased bucket
// at h == 0. A somatic mutation arises on one chromosome copy in one clone, so
// it should occupy exactly one of those four slots, and a tumour one.
//
// Derived on demand from read_var_cr / read_var_profile rather than maintained
// beside the pooled counters in the hot counting path. That costs one pass over
// the reads overlapping the variant -- the same walk collect_somatic_var_aux_info()
// already does -- and in exchange the pooled single-sample numbers stay
// bit-identical by construction, because there is no second counting path to
// keep in step.
typedef struct var_slot_counts_t {
    int alt[LONGCALLD_MAX_SAMPLES][LONGCALLD_N_HAP_SLOTS]; // read_var_profile allele 1
    int ref[LONGCALLD_MAX_SAMPLES][LONGCALLD_N_HAP_SLOTS]; // allele 0
    // Reads that are neither: allele -1 (some third allele) and -2 (alt, but the
    // base failed the quality floor). They are kept apart rather than folded
    // into ref, both because calling a low-quality alt "reference" is actively
    // misleading and because alt+ref must stay the same denominator the pooled
    // total_cov uses, or the per-haplotype VAF would be computed against a
    // larger base than the pooled one.
    int oth[LONGCALLD_MAX_SAMPLES][LONGCALLD_N_HAP_SLOTS];
    // Phase block the contributing reads belong to. A variant excluded from
    // phasing -- which every somatic variant is, since anchoring is restricted
    // to sites that are germline het in the normal -- never gets a phase_set of
    // its own, yet it is still decided on that block's haplotypes and must
    // report its PS.
    hts_pos_t ps;
    int n_samples;
} var_slot_counts_t;

// Depth and alt count for one slot, and for a whole sample (summed over
// haplotypes, including the unphased bucket).
static inline int slot_dp(const var_slot_counts_t *c, int s, int h) { return c->alt[s][h] + c->ref[s][h]; }
static inline int smp_alt(const var_slot_counts_t *c, int s) { return c->alt[s][0] + c->alt[s][1] + c->alt[s][2]; }
static inline int smp_dp(const var_slot_counts_t *c, int s) { return slot_dp(c,s,0) + slot_dp(c,s,1) + slot_dp(c,s,2); }
static inline int slot_oth(const var_slot_counts_t *c, int s, int h) { return c->oth[s][h]; }

typedef struct bam_chunk_t {
    // input
    // tid = pl->reg_chunks[reg_chunk_i].reg_tids[reg_i] 
    // reg_beg = pl->reg_chunks[reg_chunk_i].reg_begs[reg_i]
    // reg_end = pl->reg_chunks[reg_chunk_i].reg_ends[reg_i]
    int reg_chunk_i, reg_i, tid; hts_pos_t reg_beg, reg_end;
    int *qual_counts; int min_qual, first_quar_qual, median_qual, third_quar_qual, max_qual; // for each read, count the number of bases with qual
    char *tname;
    // ref_seq: 
    //   reference sequence for this chunk, size: ref_end-ref_beg+1
    //   if reg_reg_cr contain >1 regions, ref_seq will be concatenated, including additonal gaps between regions
    char *ref_seq; hts_pos_t ref_beg, ref_end, whole_ref_len; // [ref_beg, ref_end]
    // reg_cr: 
    //   reference region for this chunk. usually just 1 region, but could be multiple regions when regions are small, or reads are long
    //   only variants within regions will be considered, variants outside will be skipped
    //   for variants across boundaries/spaning multiple regions, they will be processed during stitching
    //   ref_beg <= reg_beg < reg_end <= ref_end, ref_seq may include additional flanking regions (50kb)
    // should always be 1 region XXX
    cgranges_t *low_comp_cr; // tandem_rep_cr;
    int n_reads, m_reads; int *ordered_read_ids; // size: m_reads, for multiple input bams, merge sort reads by pos, end, NM, name
    uint64_t *read_sigs; // size: m_reads, identity of the record loaded at each read_i; see bam_record_sig()
    uint8_t *read_sample_idx; // size: m_reads, which input sample each read came from
    uint8_t *read_mapq, *read_is_rev; // size: m_reads; outlive the BAM records, for the quality model
    int *bam_read_off; // size: n_bam, the read_i at which each input file's reads start.
                       // The loader fills chunk->reads in BAM-major order, so the writer
                       // needs this to recover read_i while emitting in coordinate order.
    int n_samples; int *n_sample_reads; // size: n_samples
    uint8_t is_invalid;  // loading this chunk failed; step 2 must not emit anything for it
    char **read_names; // size: m_reads, for output
    int n_bam; // for each input bam, record the number of region-overlapping reads
    int *n_up_ovlp_reads, *n_down_ovlp_reads, *n_up_ovlp_skip_reads, *n_down_ovlp_skip_reads; // number of reads overlapping with up/downstream bam chunk
    int **up_ovlp_read_i, **down_ovlp_read_i;
    bam1_t **reads;
    // intermidiate
    int *n_clean_agree_snps, *n_clean_conflict_snps; // size: m_reads; XXX include both het and hom clean vars
    uint8_t *is_ont_palindrome; // size: m_reads, 1: palindromic read, 0: non-palindromic read
    digar_t *digars; uint8_t *is_skipped, *is_skipped_for_somatic; // size: m_reads, is_skipped: wrong mapping, low qual, etc.
    // variant-related
    // n_cand_vars: including candidate germline variants and somatic variants
    // for germline variants, including clean region and noisy region
    // for somatic, only non-noisy regions, we don't call somatic variants in noisy regions (for now)
    int n_cand_vars; cand_var_t *cand_vars;
    int *var_i_to_cate;
    var_slot_counts_t *var_slots; int n_var_slots; // size: n_cand_vars; paired mode only
    //   Artifact-check bitmasks for candidates that could be somatic, indexed by
    //   candidate. Gathered at the end of step 1 because the evidence they rest
    //   on -- read_var_cr, the digars, the read-variant profiles -- is freed by
    //   bam_chunk_mid_free() before step 2 reaches make_variants(). Two per
    //   candidate: one restricted to the haplotype carrying the alt, one over
    //   every overlapping read, for a call decided on pooled counts.
    int *som_art_mask_hap, *som_art_mask_pooled; int n_som_art; // size: n_cand_vars; paired mode only
    read_var_profile_t *read_var_profile; cgranges_t *read_var_cr;
    cgranges_t *var_noisy_read_cov_cr, *var_noisy_read_err_cr;
    int *var_noisy_read_marks, var_noisy_read_mark_id;
    // noisy regions
    // right now: not work with cooridinate > 2G (pow(2,31)), use noisy_reg_beg-reg_beg+1 instead for coordinates > 2G
    cgranges_t *chunk_noisy_regs; // merged noisy regions for all reads
    int **noisy_reg_to_reads, *noisy_reg_to_n_reads; // size: chunk->chunk_noisy_regs->n_r
    // output
    uint8_t flip_hap; // flip haplotype for this chunk, chromosome-wise global parameter
    hts_pos_t flip_pre_PS, flip_cur_PS; // so all the haplotypes need to be output sequentially to keep consistency
    // read-wise
    // phase_score: used to determine low-qual phased reads, which should not be used for somatic variant calling
    int *phase_scores, *haps; hts_pos_t *phase_sets; // size: m_reads 
} bam_chunk_t; // reg-based bam_chunk_t

struct call_var_pl_t;
struct call_var_opt_t;
struct cand_var_t;
struct var_site_t;
struct read_var_profile_t;
struct call_var_io_aux_t;

// static inline int double_check_digar(digar_t *digar) {
static inline int double_check_digar(digar1_t *digars, int n_digar) {
    if (n_digar == 0) return 0;
    for (int i = n_digar-1; i > 0; --i) {
        int qi, last_i = i-1; int last_qi = digars[last_i].qi;
        if (digars[last_i].type == BAM_CEQUAL ||
            digars[last_i].type == BAM_CMATCH ||
            digars[last_i].type == BAM_CDIFF ||
            digars[last_i].type == BAM_CINS ||
            digars[last_i].type == BAM_CSOFT_CLIP ||
            digars[last_i].type == BAM_CHARD_CLIP) {
                qi = last_qi + digars[last_i].len;
        } else qi = last_qi;
        if (qi != digars[i].qi) {
            // fprintf(stderr, "Error: digars[%d]=%d, digars[%d]=%d\n", i, digars[i].qi, last_i, qi);
            return 1;
        }
    }
    return 0;
}

static inline int digar2qlen(digar_t *digar) {
    if (digar->n_digar == 0) return 0;
    int qlen = digar->digars[digar->n_digar-1].qi;
    if (digar->digars[digar->n_digar-1].type == BAM_CEQUAL ||
        digar->digars[digar->n_digar-1].type == BAM_CMATCH ||
        digar->digars[digar->n_digar-1].type == BAM_CDIFF ||
        digar->digars[digar->n_digar-1].type == BAM_CINS ||
        digar->digars[digar->n_digar-1].type == BAM_CSOFT_CLIP ||
        digar->digars[digar->n_digar-1].type == BAM_CHARD_CLIP) {
        qlen += digar->digars[digar->n_digar-1].len;
    }
    return qlen;
}

static inline int is_in_noisy_reg(hts_pos_t pos, cgranges_t *noisy_regs) {
    int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
    ovlp_n = cr_overlap(noisy_regs, "cr", pos, pos+1, &ovlp_b, &max_b);
    free(ovlp_b);
    return (ovlp_n > 0);
}

static inline int is_overlap_cr(char *tname, hts_pos_t beg, hts_pos_t end, cgranges_t *reg_cr) {
    int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
    ovlp_n = cr_overlap(reg_cr, tname, beg-1, end, &ovlp_b, &max_b);
    free(ovlp_b);
    return (ovlp_n > 0);
}

static inline int is_overlap_reg(hts_pos_t beg, hts_pos_t end, hts_pos_t reg_beg, hts_pos_t reg_end) {
    if (beg > reg_end || end < reg_beg) return 0;
    return 1;
}

digar1_t *push_digar_alt_seq(digar1_t *digar, int *n_digar, int *m_digar, digar1_t d);
digar1_t *push_digar0(digar1_t *digar, int *n_digar, int *m_digar, digar1_t d);
void push_digar1(digar_t *digar, digar1_t d);
void free_digar1(digar1_t *digar1, int n_digar);
int get_aux_int_from_bam(bam1_t *b, const char *tag);
char *get_aux_str_from_bam(bam1_t *b, const char *tag);
void print_digar1(digar1_t *digar, int n_digar, FILE *fp);
void print_digar(digar_t *digar, FILE *fp);
int has_equal_X_in_bam_cigar(bam1_t *read);
int has_cs_in_bam(bam1_t *b);
int has_MD_in_bam(bam1_t *b);
int collect_reg_digars_var_seqs(bam_chunk_t *chunk, int read_i, hts_pos_t reg_beg, hts_pos_t reg_end, digar1_t *reg_digars, uint8_t **reg_var_seqs, int *fully_cover);
int collect_digar_from_eqx_cigar(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar);
int collect_digar_from_cs_tag(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar);
int collect_digar_from_MD_tag(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar);
int collect_digar_from_ref_seq(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar);
int update_cand_vars_from_digar(const struct call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_var_sites, struct var_site_t *var_sites, struct cand_var_t *cand_vars);
void update_read_var_profile_with_allele(int var_i, int allele_i, int alt_qi, read_var_profile_t *read_var_profile);
int update_read_vs_all_var_profile_from_digar(const struct call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_cand_vars, struct cand_var_t *cand_vars, int *var_i_to_cate, struct read_var_profile_t *read_var_profile);
int update_read_vs_somatic_var_profile_from_digar(const struct call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_cand_vars, struct cand_var_t *cand_vars, struct read_var_profile_t *read_var_profile);

int collect_ref_seq_bam_main(const struct call_var_pl_t *pl, struct call_var_io_aux_t *io_aux, int reg_chunk_i, int reg_i, bam_chunk_t *chunks);

// A phased-alignment output, with a bounded reorder buffer.
//
// --refine-aln rewrites a record's POS from its realigned CIGAR, which can move
// it LEFT. The writer visits records in input (coordinate) order, so the output
// stream is "nearly sorted": a record can only ever move earlier, and only by
// less than the length of the noisy region it was realigned inside. Buffering a
// window that wide and emitting in sorted order is therefore enough to keep the
// file genuinely coordinate-sorted, without sorting whole chunks (a 500 kb
// chunk at 100x is tens of Mbp of reads).
//
// Before this existed, --refine-aln emitted records out of order while @HD still
// claimed SO:coordinate -- 6 of 587 on the bundled ONT data -- so the file was
// silently invalid and samtools index refused it.
typedef struct aln_out_t {
    htsFile *fp;
    bam_hdr_t *hdr;          // header this file was opened with; write with THIS, not a per-input header
    uint8_t is_cram;
    uint8_t no_reorder;      // pass records straight through and declare SO:unsorted
    hts_pos_t window;        // upper bound on how far left a POS can move; 0 = no reordering needed
    bam1_t **buf; int n_buf, m_buf;
    int32_t last_tid; hts_pos_t last_pos;  // last record emitted, for the monotonicity check
    int64_t n_written;
    int64_t n_unrefined_supp; // supplementary records written without CIGAR refinement
    uint8_t seen_unplaced;    // an unmapped (tid < 0) record has been written; nothing placed may follow
} aln_out_t;

aln_out_t *aln_out_init(htsFile *fp, bam_hdr_t *hdr, int is_cram, int no_reorder, hts_pos_t window);
void aln_out_push(aln_out_t *o, bam1_t *b, const char *tname);
void aln_out_flush(aln_out_t *o);
void aln_out_free(aln_out_t *o);

void collect_all_var_slot_counts(bam_chunk_t *chunk);

// Does this run need chunk->read_names[] populated?
//
// Query names are needed to report supporting reads in the VCF, and also to let
// a supplementary alignment inherit its read's HP/PS when writing the phased
// alignment output. This predicate is used at all four sites that allocate,
// grow, fill and free the array; they must agree, or the result is a leak or a
// free of uninitialised pointers.
#define need_read_names(opt) ((opt)->output_var_rnames || (opt)->output_sv_rnames || \
                              (opt)->output_somatic_var_rnames || (opt)->out_aln_fp != NULL)

// Is this record one that the variant caller takes into chunk->reads[]?
//
// This predicate decides the mapping between BAM records and read_i, and it is
// evaluated TWICE on the same region: once in collect_ref_seq_bam_main() when
// the chunk is loaded, and again in write_read_to_bam() when the phased BAM is
// written. The two passes must agree exactly, or every read_i in the chunk is
// off and HP, PS and refined CIGARs get written onto the wrong records with no
// error. It lived as two copy-pasted expressions; keep it as one function.
static inline int bam_record_is_callable(const bam1_t *b, int min_mq) {
    if (b->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) return 0;
    if (b->core.qual < min_mq) return 0;
    return 1;
}

// A cheap identity for a BAM record: FNV-1a over the query name, mixed with
// the flag and position. Used to assert that the record the writer is about to
// tag really is the record that was loaded at that read_i.
static inline uint64_t bam_record_sig(const bam1_t *b) {
    uint64_t h = 1469598103934665603ULL; // FNV-1a offset basis
    const char *q = bam_get_qname(b);
    for (; *q; ++q) { h ^= (uint64_t)(unsigned char)(*q); h *= 1099511628211ULL; }
    h ^= (uint64_t)b->core.flag;      h *= 1099511628211ULL;
    h ^= (uint64_t)(b->core.pos + 1); h *= 1099511628211ULL;
    return h;
}

int write_read_to_bam(bam_chunk_t *chunk, const struct call_var_opt_t *opt, const struct call_var_io_aux_t *io_aux);
int write_untagged_read_to_bam(bam_chunk_t *chunk, aln_out_t *out, bam1_t *read, void *phase_by_qname);
void bam_chunk_mid_free(bam_chunk_t *chunk, const struct call_var_opt_t *opt);
void bam_chunks_mid_free(bam_chunk_t *chunks, int n_chunks, const struct call_var_opt_t *opt);
void bam_chunk_post_free(bam_chunk_t *chunk, const struct call_var_opt_t *opt);
void bam_chunks_post_free(bam_chunk_t *chunks, int n_chunks, const struct call_var_opt_t *opt);
struct read_var_profile_t *init_read_var_profile(int n_reads, int n_total_vars);
struct read_var_profile_t *init_read_var_profile_with_ids(int n_reads, int *read_ids, int n_total_vars);
struct read_var_profile_t *init_read_var_profile_banded(int n_reads, int *read_ids, const int *caps);
int read_var_band_width(struct cand_var_t *cand_vars, int n_cand_vars, const digar_t *digar);
void free_read_var_profile(struct read_var_profile_t *p, int n_reads);
char *extract_sample_name_from_bam_header(bam_hdr_t *header);

// int get_xid_from_eqx_cigar(bam1_t *read, const struct call_var_opt_t *opt, xid_t **mis_bases);

// int get_mis_bases_from_MD_tag(bam1_t *read, int min_bq, x_base_t **mis_bases);

// (beg, end]
static inline int is_in_reg(cgranges_t *reg_cr, char *tname, hts_pos_t beg, hts_pos_t end) {
    int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
    ovlp_n = cr_overlap(reg_cr, tname, beg, end, &ovlp_b, &max_b);
    free(ovlp_b);
    // fprintf(stderr, "is_in_reg: %s:%d-%d : %d\n", tname, beg, end, ovlp_n);
    // for (int i = 0; i < reg_cr->n_r; ++i) {
        // fprintf(stderr, "reg_cr: %s:%d-%d\n", reg_cr->ctg[cr_label(reg_cr, i)].name, cr_start(reg_cr, i), cr_end(reg_cr, i));
    // }
    return (ovlp_n > 0);
}

#ifdef __cplusplus
}
#endif

#endif // end of LONGCALLD_BAM_UTILS_H
