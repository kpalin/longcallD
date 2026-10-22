#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#ifdef __linux__
#include <sys/sysinfo.h>  // Linux-specific
#elif __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>    // macOS-specific
#include <mach/mach.h>     // For memory information
#else
#error "Unsupported platform"
#endif
#include "main.h"
#include "call_var_main.h"
#include "bam_utils.h"
#include "vcf_utils.h"
#include "utils.h"
#include "seq.h"
#include "collect_var.h"
#include "kthread.h"
#include "align.h"
#include "math_utils.h"
#include "kmer.h"

extern int LONGCALLD_VERBOSE;

const struct option call_var_opt [] = {
    // long options
    { "amb-base", 0, NULL, 0},
    { "hifi", 0, NULL, 0},
    { "ont", 0, NULL, 0},
    { "region-file", 1, NULL, 0},
    { "regions-file", 1, NULL, 0},
    { "all-ctg", 0, NULL, 0},
    { "autosome-XY", 0, NULL, 0},
    { "autosome", 0, NULL, 0},
    { "mosaic", 0, NULL, 0},
    { "refine-aln", 0, NULL, 0},
    { "out-var-rnames", 0, NULL, 0},
    { "out-sv-rnames", 0, NULL, 0},
    { "out-som-var-rnames", 0, NULL, 0},
    // { "alpha", 1, NULL, 0},
    // { "beta", 1, NULL, 0},
    { "max-somvar", 1, NULL, 0},
    { "som-alt", 1, NULL, 0},
    { "som-mei-alt", 1, NULL, 0},

    { "exclude-ctg", 1, NULL, 'E'},
    { "extra-bam", 1, NULL, 'X'},
    { "normal", 1, NULL, 'N'},
    { "input-is-list", 0, NULL, 'L'},

    // short options
    { "ref-idx", 1, NULL, 'r' },
    { "out-vcf", 1, NULL, 'o'},
    { "out-type", 1, NULL, 'O'},
    { "min-sv-len", 1, NULL, 'l'},
    { "out-sam", 1, NULL, 'S' },
    { "out-bam", 1, NULL, 'b' },
    { "out-cram", 1, NULL, 'C' },
    { "no-vcf-header", 0, NULL, 'H'},
    { "somatic", 0, NULL, 's'},
    { "methylation", 0, NULL, 'm'},
    { "gap-aln", 1, NULL, 'g'},
    // { "out-bam", 1, NULL, 'b' },
    { "sample-name", 1, NULL, 'n'},
    { "trans-elem", 1, NULL, 'T'},

    { "min-cov", 1, NULL, 'c' },
    { "alt-cov", 1, NULL, 'd' },
    { "alt-ratio", 1, NULL, 'a' },
    { "min-mapq", 1, NULL, 'M'}, // map quality
    { "min-bq", 1, NULL, 'B'}, // base quality
    { "max-cov", 1, NULL, 0}, // total coverage; long-only, 'C' is --out-cram
    // { "max-ploidy", 1, NULL, 'p' },

    { "max-xgap", 1, NULL, 'x' },
    { "win-size", 1, NULL, 'w'},
    // Long-only below: each of these previously carried a short code that either
    // collided with a different option ('c' with --min-cov, 'C' with --out-cram)
    // or was absent from opt_str entirely ('F'), so the option silently did the
    // wrong thing or fell through to `default`. --noisy-rat, --noisy-flank and
    // --full-read are gone: their backing fields were retired from
    // call_var_opt_t and their call sites now read opt->min_af / opt->min_dp,
    // so there is nothing left for them to set.
    { "merge-dis", 1, NULL, 0},
    { "end-clip", 1, NULL, 0},
    { "clip-flank", 1, NULL, 0},
    { "no-reorder", 0, NULL, 0},
    { "min-normal-phase-dp", 1, NULL, 0},
    { "min-normal-dp", 1, NULL, 0},
    { "min-tumor-dp", 1, NULL, 0},
    { "min-normal-hap-dp", 1, NULL, 0},
    { "min-tumor-hap-dp", 1, NULL, 0},
    { "max-hap-alt", 1, NULL, 0},
    { "max-poa-cells", 1, NULL, 0},
    { "som-pval", 1, NULL, 0},
    { "no-out-aln-passthrough", 0, NULL, 0},
    { "out-normal-bam", 1, NULL, 0},
    { "out-normal-cram", 1, NULL, 0},
    { "out-normal-sam", 1, NULL, 0},
    { "tumor-name", 1, NULL, 0},
    { "tumour-name", 1, NULL, 0},
    { "normal-name", 1, NULL, 0},
    { "hap-read", 1, NULL, 'p' },

    { "threads", 1, NULL, 't' },
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "verbose", 1, NULL, 'V' },
    { 0, 0, 0, 0}
};

int get_num_processors() {
#ifdef __linux__
    return get_nprocs();
#elif __APPLE__
    int nm[2] = {CTL_HW, HW_AVAILCPU};
    int ncpu;
    size_t len = sizeof(ncpu);

    if (sysctl(nm, 2, &ncpu, &len, NULL, 0) == -1 || ncpu < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &ncpu, &len, NULL, 0);
        if (ncpu < 1) ncpu = 1;
    }
    return ncpu;
#else
    #error "Unsupported platform"
#endif
}

void set_hifi_opt(call_var_opt_t *opt) {
    opt->is_pb_hifi = 1; opt->is_ont = 0;
    // HiFi: 10/200
    opt->noisy_reg_max_xgaps = LONGCALLD_NOISY_REG_MAX_XGAPS; //10
    opt->noisy_reg_slide_win = LONGCALLD_NOISY_REG_HIFI_SLIDE_WIN; // 200bp slide window
}

void set_ont_opt(call_var_opt_t *opt) {
    opt->is_ont = 1; opt->is_pb_hifi = 0;
    opt->strand_bias_pval = LONGCALLD_STRAND_BIAS_PVAL_ONT;
    // ONT: 10/50
    opt->noisy_reg_max_xgaps = LONGCALLD_NOISY_REG_MAX_XGAPS; // 10
    opt->noisy_reg_slide_win = LONGCALLD_NOISY_REG_ONT_SLIDE_WIN; // 50bp slide window
}

// extreme low allele frequency < 5% or even 1%
// void set_mosaic_opt(call_var_opt_t *opt) {
//     opt->somatic_beta_alpha = 1; opt->somatic_beta_beta = 10;
//     opt->min_somatic_log_beta_binom = -3.0; // log10(0.001)
// }

// // moderate low allele frequency, > 5%, (tumor purity: 40~100% -> HET var: 20%~50% allele frequency)
// void set_tumor_somatic_opt(call_var_opt_t *opt) {
//     opt->somatic_beta_alpha = 2; opt->somatic_beta_beta = 5;
//     opt->min_somatic_log_beta_binom = -4.0; // log10(0.0001) ???
// }

call_var_opt_t *call_var_init_para(void) {
    call_var_opt_t *opt = (call_var_opt_t*)_err_malloc(sizeof(call_var_opt_t));

    opt->sample_name = NULL;
    opt->ref_fa_fn = NULL; opt->ref_fa_fai_fn = NULL; opt->reg_bed_fn = NULL;
    opt->input_is_list = 0; opt->n_in_bam_fn = 1; opt->m_in_bam_fn = 4; opt->in_bam_fns = (char**)malloc(4 * sizeof(char*));
    for (int i = 0; i < 4; i++) opt->in_bam_fns[i] = NULL;
    opt->n_normal_bam_fn = 0; opt->m_normal_bam_fn = 0; opt->normal_bam_fns = NULL;
    opt->n_tumor_bam_fn = 0; opt->n_samples = 1; opt->bam_i_to_sample = NULL;
    opt->sample_names = NULL; opt->paired_mode = 0;
    opt->tumor_name = NULL; opt->normal_name = NULL;

    opt->is_pb_hifi = 1; opt->is_ont = 0;
    opt->only_autosome = 0; opt->only_autosome_XY = 1; 
    opt->n_exc_tnames = 0; opt->exc_tnames = NULL;

    // opt->region_list = NULL; opt->region_is_file = 0;
    opt->max_ploid = LONGCALLD_DEF_PLOID;
    opt->min_mq = LONGCALLD_MIN_CAND_MQ;
    opt->min_bq = LONGCALLD_MIN_CAND_BQ;
    opt->min_dp = LONGCALLD_MIN_CAND_DP;
    opt->min_alt_dp = LONGCALLD_MIN_ALT_DP;
    opt->min_af = LONGCALLD_MIN_CAND_AF;
    opt->max_af = LONGCALLD_MAX_CAND_AF;
    opt->min_normal_phase_dp = 4;
    opt->min_normal_dp = 8; opt->min_tumor_dp = 8;
    opt->min_normal_hap_dp = 4; opt->min_tumor_hap_dp = 4;
    opt->max_hap_alt = 0;
    // The inherited mosaic default of 0.05 is far too permissive at these
    // depths; retune against the E3 precision/recall curve.
    opt->som_pval = 1e-3;
    opt->max_noisy_reg_cov = LONGCALLD_MAX_NOISY_REG_COV;
    // somatic/mosaic var
    opt->min_somatic_dis_to_var = LONGCALLD_MIN_SOMATIC_DIS_TO_VAR;
    opt->min_somatic_dis_to_homopolymer_indel_error = LONGCALLD_MIN_SOMATIC_DIS_TO_HP_INDEL_ERROR;
    opt->min_somatic_dis_to_seq_error = LONGCALLD_MIN_SOMATIC_DIS_TO_SEQ_ERROR;
    opt->min_somatic_fisher_pval = LONGCALLD_MIN_SOMATIC_FISHER_PVAL;
    opt->min_somatic_alt_dp = LONGCALLD_MIN_SOMATIC_ALT_DP;
    opt->min_somatic_hap_dp = LONGCALLD_MIN_SOMATIC_HAP_READS;
    opt->min_somatic_te_dp = LONGCALLD_MIN_SOMATIC_TE_ALT_DP;
    // opt->min_somatic_af = LONGCALLD_MIN_SOMATIC_AF;
    // opt->max_somatic_alt_af = LONGCALLD_MAX_SOMATIC_ALT_AF;
    // opt->min_somatic_log_beta_binom = LONGCALLD_MIN_SOMATIC_LOG_BETA_BINOM;
    // opt->somatic_beta_alpha = LONGCALLD_SOMATIC_BETA_ALPHA;
    // opt->somatic_beta_beta = LONGCALLD_SOMATIC_BETA_BETA;
    opt->somatic_win = LONGCALLD_SOMATIC_WIN;
    opt->somatic_win_max_vars = LONGCALLD_SOMATIC_WIN_MAX_VARS;

    opt->noisy_reg_merge_dis = LONGCALLD_NOISY_REG_MERGE_DIS;
    opt->noisy_reg_flank_len = LONGCALLD_NOISY_REG_FLANK_LEN;

    opt->noisy_reg_max_xgaps = LONGCALLD_NOISY_REG_MAX_XGAPS;
    opt->noisy_reg_slide_win = LONGCALLD_NOISY_REG_HIFI_SLIDE_WIN;

    opt->end_clip_reg = LONGCALLD_NOISY_END_CLIP;
    opt->end_clip_reg_flank_win = LONGCALLD_NOISY_END_CLIP_WIN;

    opt->max_var_ratio_per_read = LONGCALLD_MAX_VAR_RATIO_PER_READ;
    opt->max_noisy_reg_len  = LONGCALLD_MAX_NOISY_REG_LEN;
    // Between the largest POA a normal region needs (6.3e8 cells, under 1 GB)
    // and the satellite call that asked for 16 GiB in one go (3.2e10).
    opt->max_poa_cells = 2000000000;
    // opt->disable_read_sampling = 0; // by default read-sampling is enabled for long noisy regions (10kb+)
    // opt->min_noisy_reg_reads = LONGCALLD_NOISY_REG_READS;
    // opt->min_noisy_reg_ratio = LONGCALLD_NOISY_REG_RATIO;
    opt->max_noisy_frac_per_read = LONGCALLD_MAX_NOISY_FRAC_PER_READ;
    opt->min_hap_full_reads = LONGCALLD_MIN_HAP_FULL_READS;
    opt->min_hap_reads = LONGCALLD_MIN_HAP_READS;
    // opt->min_no_hap_full_reads = LONGCALLD_MIN_NO_HAP_FULL_READS;

    opt->match = LONGCALLD_MATCH_SCORE;
    opt->mismatch = LONGCALLD_MISMATCH_SCORE;
    opt->gap_open1 = LONGCALLD_GAP_OPEN1_SCORE;
    opt->gap_ext1 = LONGCALLD_GAP_EXT1_SCORE;
    opt->gap_open2 = LONGCALLD_GAP_OPEN2_SCORE;
    opt->gap_ext2 = LONGCALLD_GAP_EXT2_SCORE;
    opt->gap_aln = LONGCALLD_GAP_LEFT_ALN;
    opt->partial_aln_ratio = LONGCALLD_PARTIAL_ALN_RATIO;
    opt->min_noisy_reg_size_to_sample_reads = LONGCALLD_MIN_NOISY_REG_SIZE_TO_SAMPLE_READS;

    opt->pl_threads = MIN_OF_TWO(CALL_VAR_PL_THREAD_N, get_num_processors());
    opt->n_threads = MIN_OF_TWO(CALL_VAR_THREAD_N, get_num_processors());

    opt->min_sv_len = LONGCALLD_MIN_SV_LEN;
    opt->min_tsd_len = LONGCALLD_MIN_TSD_LEN; opt->max_tsd_len = LONGCALLD_MAX_TSD_LEN;
    opt->min_polya_len = LONGCALLD_MIN_POLYA_LEN; opt->min_polya_ratio = LONGCALLD_MIN_POLYA_RATIO;
    opt->output_sv_rnames = 0; opt->output_somatic_var_rnames = 0; opt->output_var_rnames = 0;

    initialize_lgamma_cache(opt);
    opt->te_seq_fn = NULL; opt->n_te_seqs = 0; opt->te_kmer_len = 15; opt->te_for_h = NULL; opt->te_rev_h = NULL;

    opt->p_error = 0.001; opt->log_p = -3.0; opt->log_1p = log10(1-opt->p_error); opt->log_2 = 0.301023;
    opt->max_gq = 60; opt->max_qual = 60;
    opt->out_vcf = NULL; opt->vcf_hdr = NULL; opt->out_vcf_fn = NULL; opt->out_vcf_type = 'v'; opt->no_vcf_header = 0; opt->out_amb_base = 0;
    opt->out_aln_fp = NULL; opt->out_aln_is_cram = 0; opt->refine_bam = 0;
    opt->out_naln_fp = NULL; opt->out_naln_is_cram = 0;
    for (int i = 0; i < LONGCALLD_MAX_SAMPLES; ++i) opt->aln_outs[i] = NULL;
    opt->no_reorder = 0; opt->out_aln_passthrough = 1;
    opt->out_somatic = 0; opt->out_methylation = 0;
    // opt->verbose = 0;
    return opt;
}

void call_var_free_para(call_var_opt_t *opt) {
    if (opt->ref_fa_fn) free(opt->ref_fa_fn);
    for (int i = 0; i < opt->n_in_bam_fn; ++i) {
        if (opt->in_bam_fns[i] != NULL) free(opt->in_bam_fns[i]);
    }
    free(opt->in_bam_fns);
    if (opt->normal_bam_fns != NULL) {
        for (int i = 0; i < opt->n_normal_bam_fn; ++i) free(opt->normal_bam_fns[i]); // NULL once moved
        free(opt->normal_bam_fns);
    }
    free(opt->bam_i_to_sample);
    if (opt->sample_names != NULL) {
        for (int i = 0; i < opt->n_samples; ++i) free(opt->sample_names[i]);
        free(opt->sample_names);
    }
    if (opt->bams != NULL) {
        for (int i = 0; i < opt->n_in_bam_fn; ++i) {
            if (opt->bams[i]) {
                if (opt->idxs[i]) hts_idx_destroy(opt->idxs[i]);
                bam_hdr_destroy(opt->headers[i]);
                sam_close(opt->bams[i]);
            }
        }
        free(opt->bams); free(opt->idxs); free(opt->headers);
    }
    if (opt->sample_name) free(opt->sample_name);
    free(opt->tumor_name); free(opt->normal_name);
    // if (opt->region_list) free(opt->region_list);
    if (opt->n_exc_tnames > 0) {
        for (int i = 0; i < opt->n_exc_tnames; ++i) free(opt->exc_tnames[i]);
        free(opt->exc_tnames);
    }
    if (opt->out_vcf_fn != NULL) free(opt->out_vcf_fn);
    if (opt->te_seq_fn != NULL) {
        free(opt->te_seq_fn);
        for (int i = 0; i < opt->n_te_seqs; ++i) {
            free(opt->te_seq_names[i]);
            kh_destroy(kmer32, opt->te_for_h[i]);
            kh_destroy(kmer32, opt->te_rev_h[i]);
        }
        free(opt->te_seq_names); free(opt->te_for_h); free(opt->te_rev_h);
    }
    free(opt);
}

void reg_chunks_free(reg_chunks_t *reg_chunks, int m_reg_chunks) {
    for (int i = 0; i < m_reg_chunks; ++i) {
        if (reg_chunks[i].reg_tids) free(reg_chunks[i].reg_tids);
        if (reg_chunks[i].reg_begs) free(reg_chunks[i].reg_begs);
        if (reg_chunks[i].reg_ends) free(reg_chunks[i].reg_ends);
    }
    free(reg_chunks);
}

void call_var_io_aux_free(call_var_io_aux_t *aux, int n) {
    for (int i = 0; i < n; ++i) {
        if (aux[i].fai) fai_destroy(aux[i].fai);
        if (aux[i].n_bam > 0) {
            for (int j = 0; j < aux[i].n_bam; ++j) {
                if (aux[i].bams[j]) {
                    if (aux[i].idxs[j]) hts_idx_destroy(aux[i].idxs[j]);
                    bam_hdr_destroy(aux[i].headers[j]);
                    sam_close(aux[i].bams[j]);
                }
            }
            free(aux[i].bams); free(aux[i].idxs); free(aux[i].headers);
        }
    }
    free(aux);
}
void call_var_free_pl(call_var_pl_t pl) {
    call_var_io_aux_free(pl.io_aux, pl.n_threads);
    reg_chunks_free(pl.reg_chunks, pl.m_reg_chunks);
}

void var_free(var_t *v) {
    if (v->vars) {
        for (int i = 0; i < v->n; ++i) {
            if (v->vars[i].ref_bases) free(v->vars[i].ref_bases);
            if (v->vars[i].alt_len) free(v->vars[i].alt_len);
            if (v->vars[i].alt_bases) {
                for (int j = 0; j < v->vars[i].n_alt_allele; ++j) free(v->vars[i].alt_bases[j]);
                free(v->vars[i].alt_bases);
            }
            if (v->vars[i].tsd_len > 0) {
                free(v->vars[i].tsd_seq);
            }
            // alt_read_i is allocated in make_variants for every variant with
            // AD[1] > 0, not just SVs, and AD[1] may be lowered afterwards -- so
            // free it unconditionally. make_variants sets it to NULL otherwise,
            // and every element in [0, v->n) passes through that assignment.
            free(v->vars[i].alt_read_i);
        }
        free(v->vars);
    }
}

void var1_free(var1_t *v) {
    free(v->alt_read_i);
    if (v->ref_bases) free(v->ref_bases);
    if (v->alt_len) free(v->alt_len);
    if (v->alt_bases) {
        for (int j = 0; j < v->n_alt_allele; ++j) free(v->alt_bases[j]);
        free(v->alt_bases);
    }
    if (v->tsd_len > 0) {
        free(v->tsd_seq);
    }
}

static void collect_bam_call_var_worker_for(void *_data, long ii, int tid) {
    call_var_step_t *step = (call_var_step_t*)_data;
    bam_chunk_t *c = step->chunks + ii;
    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "[%s] thread-id: %d, chunk: %ld (%d) ... \n", __func__, tid, ii, step->n_chunks);
    if (collect_ref_seq_bam_main(step->pl, step->pl->io_aux+tid, step->pl->reg_chunk_i, ii, c) >= 0 && !c->is_invalid) {
        collect_var_main(step->pl, c);
    } else {
        _err_warning("Skipping chunk %s:%" PRId64 "-%" PRId64 ": failed to load reads.\n", c->tname, c->reg_beg, c->reg_end);
        c->is_invalid = 1;
    }
    bam_chunk_mid_free(c, step->pl->opt);
    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "[%s] thread-id: %d, chunk: %ld (%d) ... done\n", __func__, tid, ii, step->n_chunks);
}

// merge variants from two chunks
static void make_var_worker_for(void *_data, long ii, int tid) {
    call_var_step_t *step = (call_var_step_t*)_data;
    bam_chunk_t *c = step->chunks+ii; var_t *var = step->vars+ii;
    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "[%s] thread-id: %d, chunk: %ld (%d), n_reads: %d\n", __func__, tid, ii, step->n_chunks, c->n_reads);
    if (c->is_invalid) return; // step->vars is calloc'd, so an untouched var_t frees cleanly
    make_var_main(step, c, var, ii);
}

static void reg_chunks_realloc(call_var_pl_t *pl) {
    if (pl->n_reg_chunks >= pl->m_reg_chunks) {
        pl->m_reg_chunks *= 2;
        pl->reg_chunks = (reg_chunks_t*)realloc(pl->reg_chunks, pl->m_reg_chunks * sizeof(reg_chunks_t));
        for (int k = pl->n_reg_chunks; k < pl->m_reg_chunks; k++) {
            pl->reg_chunks[k].n_regions = 0;
            pl->reg_chunks[k].m_regions = 512;
            pl->reg_chunks[k].reg_tids = (int*)malloc(pl->reg_chunks[k].m_regions * sizeof(int));
            pl->reg_chunks[k].reg_begs = (hts_pos_t*)malloc(pl->reg_chunks[k].m_regions * sizeof(hts_pos_t));
            pl->reg_chunks[k].reg_ends = (hts_pos_t*)malloc(pl->reg_chunks[k].m_regions * sizeof(hts_pos_t));
        }
    }
}

void set_exclude_ctg(call_var_opt_t *opt, char *ctg) {
    opt->exc_tnames = (char**)realloc(opt->exc_tnames, (opt->n_exc_tnames+1) * sizeof(char*));
    opt->exc_tnames[opt->n_exc_tnames] = strdup(ctg);
    opt->n_exc_tnames++;
}

void set_extra_input_bam(call_var_opt_t *opt, char *fn) {
    if (opt->n_in_bam_fn >= opt->m_in_bam_fn) {
        opt->m_in_bam_fn *= 2;
        opt->in_bam_fns = (char**)realloc(opt->in_bam_fns, opt->m_in_bam_fn * sizeof(char*));
    }
    opt->in_bam_fns[opt->n_in_bam_fn] = retrieve_full_url(fn);
    opt->n_in_bam_fn++;
}

void set_normal_input_bam(call_var_opt_t *opt, char *fn) {
    if (opt->n_normal_bam_fn >= opt->m_normal_bam_fn) {
        opt->m_normal_bam_fn = opt->m_normal_bam_fn ? opt->m_normal_bam_fn * 2 : 4;
        opt->normal_bam_fns = (char**)realloc(opt->normal_bam_fns, opt->m_normal_bam_fn * sizeof(char*));
    }
    opt->normal_bam_fns[opt->n_normal_bam_fn++] = retrieve_full_url(fn);
}

// Append the -N files to in_bam_fns and build bam_i_to_sample.
//
// Must run AFTER set_bam_list_file(), which permutes in_bam_fns: doing it
// earlier would scramble which index belongs to which sample. Appending the
// normals last makes the mapping a single comparison and keeps it that way.
void finalise_sample_layout(call_var_opt_t *opt) {
    opt->n_tumor_bam_fn = opt->n_in_bam_fn;
    for (int i = 0; i < opt->n_normal_bam_fn; ++i) {
        if (opt->n_in_bam_fn >= opt->m_in_bam_fn) {
            opt->m_in_bam_fn *= 2;
            opt->in_bam_fns = (char**)realloc(opt->in_bam_fns, opt->m_in_bam_fn * sizeof(char*));
        }
        opt->in_bam_fns[opt->n_in_bam_fn++] = opt->normal_bam_fns[i];
        opt->normal_bam_fns[i] = NULL; // ownership moved
    }
    opt->n_samples = opt->n_normal_bam_fn > 0 ? 2 : 1;
    opt->bam_i_to_sample = (int*)malloc(opt->n_in_bam_fn * sizeof(int));
    for (int i = 0; i < opt->n_in_bam_fn; ++i)
        opt->bam_i_to_sample[i] = (i < opt->n_tumor_bam_fn) ? LONGCALLD_SAMPLE_TUMOR : LONGCALLD_SAMPLE_NORMAL;
    if (opt->n_samples > 1) {
        opt->paired_mode = 1;
        // Paired mode implies the low-allele-fraction candidate path. Candidate
        // selection measures allele fraction on the POOLED tumour+normal depth,
        // and a somatic variant is low-fraction there by construction: a clonal
        // het somatic variant at 100% purity with equal tumour and normal depth
        // sits at about 0.25 pooled, and a subclonal one well below that --
        // under the 0.20 germline threshold either way. Without this, paired
        // mode would systematically miss exactly the variants it exists to
        // find. Measured: a synthetic somatic SNV at tumour VAF 0.25 (pooled
        // 0.138) is not even a candidate without it.
        opt->out_somatic = 1;
        _err_info("Paired input: %d tumour file(s), %d normal file(s)\n", opt->n_tumor_bam_fn, opt->n_normal_bam_fn);
    }
}

void set_bam_list_file(call_var_opt_t *opt, char *fn) {
    // open bam list file
    FILE *fp = fopen(fn, "r");
    if (fp == NULL) _err_error_exit("Cannot open bam list file: %s\n", fn);
    char line[1024];
    int n_file = 0, n_extra_file = opt->n_in_bam_fn-1;
    while (fgets(line, sizeof(line), fp)) {
        char *bam_fn = strtok(line, "\n\r");
        if (bam_fn != NULL) {
            if (n_file == 0) {
                opt->in_bam_fns[0] = retrieve_full_url(bam_fn);
            } else {
                set_extra_input_bam(opt, bam_fn);
            }
            n_file++;
        }
    }
    // change [1 .. n_extra_file .. n_file-1] to [1 .. n_file-1 .. n_extra_file]
    if (n_extra_file > 0) {
        char **extra_fns = (char**)malloc(n_extra_file * sizeof(char*));
        for (int i = 0; i < n_extra_file; ++i) {
            extra_fns[i] = opt->in_bam_fns[1 + i];
        }
        for (int i = 1; i < n_file; ++i) {
            opt->in_bam_fns[i] = opt->in_bam_fns[n_extra_file + i];
        }
        for (int i = 0; i < n_extra_file; ++i) {
            opt->in_bam_fns[n_file + i] = extra_fns[i];
        }
        free(extra_fns);
    }
    fclose(fp);
}

// Enum for chromosome types
typedef enum {
    AUTOSOME,        // Chromosomes 1-22
    SEX_CHROMOSOME,  // X or Y
    OTHER            // Mitochondrial (MT/M) or unplaced/unusual chromosomes
} ChromosomeType;

ChromosomeType classify_chromosome(const char *chr) {
    char chr_copy[1000];  // Buffer to hold the extracted chromosome name
    strncpy(chr_copy, chr, sizeof(chr_copy) - 1);
    chr_copy[sizeof(chr_copy) - 1] = '\0';

    // Find first occurrence of ':', indicating a genomic region
    char *colon_pos = strchr(chr_copy, ':');
    if (colon_pos != NULL) {
        *colon_pos = '\0';  // Truncate at ':'
    }

    // Remove "chr" prefix if present
    if (strncmp(chr_copy, "chr", 3) == 0) {
        memmove(chr_copy, chr_copy + 3, strlen(chr_copy) - 2);  // Shift left
    }

    // Check for sex chromosomes
    if (strcmp(chr_copy, "X") == 0 || strcmp(chr_copy, "Y") == 0) {
        return SEX_CHROMOSOME;
    }

    // Check for other chromosomes (e.g., Mitochondrial DNA)
    if (strcmp(chr_copy, "MT") == 0 || strcmp(chr_copy, "M") == 0) {
        return OTHER;
    }

    // Convert string to integer and check for autosome range
    char *endptr;
    long num = strtol(chr_copy, &endptr, 10);

    // if (*endptr == '\0' && num >= 1 && num <= 22) {
    if (*endptr == '\0' && num >= 1) {
        return AUTOSOME;
    }
    return OTHER; // Anything else falls under "OTHER"
}

static void reg_chunks_region_realloc(reg_chunks_t *reg_chunks) {
    if (reg_chunks->n_regions >= reg_chunks->m_regions) {
        reg_chunks->m_regions *= 2;
        reg_chunks->reg_tids = (int*)realloc(reg_chunks->reg_tids, reg_chunks->m_regions * sizeof(int));
        reg_chunks->reg_begs = (hts_pos_t*)realloc(reg_chunks->reg_begs, reg_chunks->m_regions * sizeof(hts_pos_t));
        reg_chunks->reg_ends = (hts_pos_t*)realloc(reg_chunks->reg_ends, reg_chunks->m_regions * sizeof(hts_pos_t));
    }
}

static int exclude_target_ctg(call_var_opt_t *opt, char *tname) {
    for (int i = 0; i < opt->n_exc_tnames; ++i) {
        if (strcmp(tname, opt->exc_tnames[i]) == 0) return 1;
    }
    return 0;
}

static int skip_target_region(call_var_opt_t *opt, char *tname) {
    ChromosomeType t = classify_chromosome(tname);
    if (opt->only_autosome && t != AUTOSOME) return 1;
    if (opt->only_autosome_XY && t != AUTOSOME && t != SEX_CHROMOSOME) return 1;
    if (exclude_target_ctg(opt, tname)) return 1;
    return 0;
}

static int collect_regions_from_region_list(call_var_opt_t *opt, call_var_pl_t *pl, hts_itr_t *iter, int n_regions, char **regions) {
    int last_tid = -1; bam_hdr_t *hdr = pl->io_aux[0].headers[0];
    int min_reg_chunks_per_run = pl->min_reg_chunks_per_run, max_reg_len_per_chunk = pl->max_reg_len_per_chunk;
    for (int i = 0; i < iter->n_reg; ++i) {
        for (int j = 0; j < iter->reg_list[i].count; ++j) { // ACGT:1234, (beg, end]:(0,4]
            int tid = iter->reg_list[i].tid; hts_pos_t chr_len = hdr->target_len[tid]; char *tname = hdr->target_name[tid];
            if (skip_target_region(opt, tname)) continue;
            // fprintf(stderr, "Region: %s:%" PRId64 "-%" PRId64 "\n", hdr->target_name[tid], iter->reg_list[i].intervals[j].beg, iter->reg_list[i].intervals[j].end);
            if (last_tid != -1 && tid != last_tid && pl->reg_chunks[pl->n_reg_chunks].n_regions >= min_reg_chunks_per_run)
                pl->n_reg_chunks++;
            reg_chunks_realloc(pl);
            hts_pos_t reg_beg = MAX_OF_TWO(1, iter->reg_list[i].intervals[j].beg + 1); // 0-base
            hts_pos_t reg_end = MIN_OF_TWO(iter->reg_list[i].intervals[j].end, chr_len);
            hts_pos_t n_regions = (reg_end - reg_beg + max_reg_len_per_chunk) / max_reg_len_per_chunk;
            reg_chunks_t *reg_chunks = &pl->reg_chunks[pl->n_reg_chunks];
            for (int reg_i = 0; reg_i < n_regions; ++ reg_i) {
                hts_pos_t beg = (hts_pos_t) reg_i * max_reg_len_per_chunk + reg_beg;
                hts_pos_t end = MIN_OF_TWO((hts_pos_t) (reg_i + 1) * max_reg_len_per_chunk + reg_beg-1, reg_end);
                // add region
                reg_chunks_region_realloc(reg_chunks);
                reg_chunks->reg_tids[reg_chunks->n_regions] = tid;
                reg_chunks->reg_begs[reg_chunks->n_regions] = beg; // [beg, end]
                reg_chunks->reg_ends[reg_chunks->n_regions] = end;
                reg_chunks->n_regions++;
            }
            last_tid = tid;
        }
    }
    // n_reg_chunks may have reached m_reg_chunks on the last iteration, so make
    // sure the slot we are about to inspect exists before reading n_regions.
    reg_chunks_realloc(pl);
    if (pl->reg_chunks[pl->n_reg_chunks].n_regions > 0) pl->n_reg_chunks++;
    hts_itr_destroy(iter);
    // Print the region_chunks
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < pl->n_reg_chunks; i++) {
            fprintf(stderr, "Chunk %d:\n", i);
            for (int j = 0; j < pl->reg_chunks[i].n_regions; j++) {
                fprintf(stderr, "\tRegion %d %s:%" PRIi64 "-%" PRIi64 "\n", j, pl->io_aux[0].headers[0]->target_name[pl->reg_chunks[i].reg_tids[j]], pl->reg_chunks[i].reg_begs[j], pl->reg_chunks[i].reg_ends[j]);
            }
        }
    }
    return pl->n_reg_chunks;
}

static int collect_regions_from_bed_file(call_var_opt_t *opt, call_var_pl_t *pl) {
    FILE *fp = fopen(opt->reg_bed_fn, "r");
    if (fp == NULL) {
        _err_error("Failed to open region bed file: %s\n", opt->reg_bed_fn);
        return -1;
    }
    char line[1024];
    // collect regions
    bam_hdr_t *hdr = pl->io_aux[0].headers[0];
    int n_regions = 0, m_regions = 1024;
    char **regions = (char**)malloc(m_regions * sizeof(char*));
    for (int i = 0; i < m_regions; ++i) regions[i] = (char*)malloc(1024);
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue; // skip header/comment
        // exclude \n in line if necessary
        if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
        char *token;
        token = strtok(line, "\t");
        if (token == NULL) continue;
        char *tname = token; int tid = bam_name2id(hdr, tname);
        if (tid == -1) {
            _err_warning("Skip unknonw region: %s\n", line);
            continue;
        }
        if (skip_target_region(opt, tname)) continue;
        hts_pos_t beg = 1, end = hdr->target_len[tid];
        token = strtok(NULL, "\t");
        if (token != NULL) {
            beg = atoi(token)+1;
            token = strtok(NULL, "\t");
            if (token != NULL) end = atoi(token);
        }
        if (beg > end || beg <= 0 || end <= 0) continue;
        if (n_regions == m_regions) {
            m_regions *= 2; // should be m, not n_regions
            regions = (char**)realloc(regions, m_regions * sizeof(char*));
            for (int i = n_regions; i < m_regions; ++i) regions[i] = (char*)malloc(1024);
        }
        sprintf(regions[n_regions], "%s:%" PRIi64 "-%" PRIi64 "", tname, beg, end);
        n_regions++;
    }
    hts_idx_t *idx = pl->io_aux[0].idxs[0];
    hts_itr_t *iter = sam_itr_regarray(idx, hdr, regions, n_regions);
    if (iter != NULL) return collect_regions_from_region_list(opt, pl, iter, n_regions, regions);
    else return -1;
}

static int collect_regions(call_var_pl_t *pl, call_var_opt_t *opt, int n_regions, char **regions) {
    int min_reg_chunks_per_run = pl->min_reg_chunks_per_run;
    int max_reg_len_per_chunk = pl->max_reg_len_per_chunk;
    pl->n_reg_chunks = 0; pl->m_reg_chunks = 100;
    pl->reg_chunks = (reg_chunks_t*)malloc(pl->m_reg_chunks * sizeof(reg_chunks_t));
    for (int i = 0; i < pl->m_reg_chunks; i++) {
        pl->reg_chunks[i].n_regions = 0;
        pl->reg_chunks[i].m_regions = 512;
        pl->reg_chunks[i].reg_tids = (int*)malloc(pl->reg_chunks[i].m_regions * sizeof(int));
        pl->reg_chunks[i].reg_begs = (hts_pos_t*)malloc(pl->reg_chunks[i].m_regions * sizeof(hts_pos_t));
        pl->reg_chunks[i].reg_ends = (hts_pos_t*)malloc(pl->reg_chunks[i].m_regions * sizeof(hts_pos_t));
    }
    // check if region(s) is provided
    if (n_regions > 0) {
        opt->only_autosome = 0, opt->only_autosome_XY=0;
        if (opt->reg_bed_fn != NULL) {
            _err_error("Both region(s) and region bed file are provided. Only region(s) will be used.\n");
        }
        hts_idx_t *idx = pl->io_aux[0].idxs[0]; bam_hdr_t *hdr = pl->io_aux[0].headers[0];
        hts_itr_t *iter = sam_itr_regarray(idx, hdr, regions, n_regions);
        if (iter == NULL) {
            _err_error("Failed to parse provided regions: ");
            for (int i = 0; i < n_regions; ++i) fprintf(stderr, "%s ", regions[i]);
            fprintf(stderr, "\n");
            _err_error("Will process the whole alignment file.\n");
        } else {
            return collect_regions_from_region_list(opt, pl, iter, n_regions, regions);
        }
    } else if (opt->reg_bed_fn != NULL) { // check if region bed file is provided
        opt->only_autosome = 0, opt->only_autosome_XY=0;
        collect_regions_from_bed_file(opt, pl);
        if (pl->n_reg_chunks == 0) {
            _err_error("Failed to parse provided region bed file: %s\n", opt->reg_bed_fn);
            _err_error("Will process the whole alignment file.\n");
        } else return pl->n_reg_chunks;
    }
    // whole genome, plit chromosomes into region_chunks
    // collect_regions_from_whole_genome(opt, pl, min_reg_chunks_per_run, max_reg_len_per_chunk);
    for (int i = 0; i < pl->io_aux[0].headers[0]->n_targets; i++) {
        reg_chunks_realloc(pl);
        // all regions in a chromosome will be in one reg_chunk
        char *tname = pl->io_aux[0].headers[0]->target_name[i];
        if (skip_target_region(opt, tname)) continue;
        int tid = i; hts_pos_t chr_len = pl->io_aux[0].headers[0]->target_len[i];
        int n_regions_chr = (chr_len + max_reg_len_per_chunk - 1) / max_reg_len_per_chunk;
        reg_chunks_t *reg_chunks = &pl->reg_chunks[pl->n_reg_chunks];
        for (int reg_i = 0; reg_i < n_regions_chr; reg_i++) {
            hts_pos_t beg = (hts_pos_t) reg_i * max_reg_len_per_chunk + 1;
            hts_pos_t end = MIN_OF_TWO((hts_pos_t) (reg_i + 1) * max_reg_len_per_chunk, chr_len);
            // add region
            reg_chunks_region_realloc(reg_chunks);
            reg_chunks->reg_tids[reg_chunks->n_regions] = tid;
            reg_chunks->reg_begs[reg_chunks->n_regions] = beg;
            reg_chunks->reg_ends[reg_chunks->n_regions] = end;
            reg_chunks->n_regions++;
        }
        if (reg_chunks->n_regions >= min_reg_chunks_per_run) {
            pl->n_reg_chunks++;
        }
    }
    // n_reg_chunks may have reached m_reg_chunks on the last iteration, so make
    // sure the slot we are about to inspect exists before reading n_regions.
    reg_chunks_realloc(pl);
    if (pl->reg_chunks[pl->n_reg_chunks].n_regions > 0) pl->n_reg_chunks++;

    // Print the region_chunks
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "Collected %d region chunks:\n", pl->n_reg_chunks);
        for (int i = 0; i < pl->n_reg_chunks; i++) {
            fprintf(stderr, "Chunk %d:\n", i);
            for (int j = 0; j < pl->reg_chunks[i].n_regions; j++) {
                fprintf(stderr, "\tRegion %d %s:%" PRIi64 "-%" PRIi64 "\n", j, pl->io_aux[0].headers[0]->target_name[pl->reg_chunks[i].reg_tids[j]], pl->reg_chunks[i].reg_begs[j], pl->reg_chunks[i].reg_ends[j]);
            }
        }
    }
    return pl->n_reg_chunks; // num of collected chromosomes/contigs
}

// only works with sorted index BAM/CRAM
static void call_var_pl_open_fa_bam(call_var_opt_t *opt, call_var_pl_t *pl, char **regions, int n_regions) {
    // input BAM file
    if (opt->n_in_bam_fn <= 0) _err_error_exit("No input BAM/CRAM files provided.\n");
    for (int j = 0; j < opt->n_in_bam_fn; ++j) {
        if (strcmp(opt->in_bam_fns[j], "-") == 0)
            _err_error_exit("Input from pipe/stdin is not supported\n");
        _err_info("Opening alignment file: %s\n", opt->in_bam_fns[j]);
    }
    pl->io_aux = (call_var_io_aux_t*)malloc(pl->n_threads * sizeof(call_var_io_aux_t));
    // Checked once the headers are read, below: every input must declare the
    // same contigs in the same order. tid is just an index into that list, so
    // two inputs that disagree would silently put records on the wrong
    // chromosome -- harmless to ignore with one input and one output, fatal
    // once tumour and normal are written to separate files.

    // multi-threading
    for (int i = 0; i < pl->n_threads; ++i){
        call_var_io_aux_t *aux = &pl->io_aux[i];
        aux->n_bam = opt->n_in_bam_fn;
        aux->bams = (samFile **)calloc(aux->n_bam, sizeof(samFile *));
        aux->headers = (bam_hdr_t **)calloc(aux->n_bam, sizeof(bam_hdr_t *));
        aux->idxs = (hts_idx_t **)calloc(aux->n_bam, sizeof(hts_idx_t *));
        if (!aux->bams || !aux->headers || !aux->idxs) _err_error_exit("Memory allocation failure\n");
        aux->fai = fai_load3(opt->ref_fa_fn, opt->ref_fa_fai_fn, NULL, FAI_CREATE);
        if (aux->fai == NULL)
            _err_error_exit("Failed to load/build reference fasta index: %s\n", opt->ref_fa_fn);

        for (int j = 0; j < aux->n_bam; ++j) {
            aux->bams[j] = sam_open(opt->in_bam_fns[j], "r"); if (aux->bams[j] == NULL) _err_error_exit("Failed to open alignment file \'%s\'\n", opt->in_bam_fns[j]);
            const htsFormat *fmt = hts_get_format(aux->bams[j]); if (!fmt) _err_error_exit("Failed to get format of alignment file \'%s\'\n", opt->in_bam_fns[j]);
            if (fmt->format != bam && fmt->format != cram) {
                sam_close(aux->bams[j]); _err_error_exit("Input alignment file must be BAM or CRAM format.\n");
            }
            aux->headers[j] = sam_hdr_read(aux->bams[j]);
            if (aux->headers[j] == NULL) {
                sam_close(aux->bams[j]);
                _err_error_exit("Failed to read alignment file header \'%s\'\n", opt->in_bam_fns[j]);
            }
            if (fmt->format == cram) {
                if (hts_set_fai_filename(aux->bams[j], opt->ref_fa_fn) != 0) {
                    fai_destroy(aux->fai); sam_hdr_destroy(aux->headers[j]); sam_close(aux->bams[j]);
                    _err_error_exit("Failed to set reference file for CRAM decoding: %s %s\n", opt->in_bam_fns[j], opt->ref_fa_fn);
                }
            }
            aux->idxs[j] = sam_index_load(aux->bams[j], opt->in_bam_fns[j]);
            if (aux->idxs[j] == NULL) { // attempt to create index if not exist
                _err_warning("Index not found for \'%s\', creating now ...\n", opt->in_bam_fns[j]);
                if (sam_index_build(opt->in_bam_fns[j], 0) < 0) {
                    fai_destroy(aux->fai); sam_hdr_destroy(aux->headers[j]); sam_close(aux->bams[j]);
                    _err_error_exit("Failed to build index for \'%s\'\n", opt->in_bam_fns[j]);
                } else {
                    aux->idxs[j] = sam_index_load(aux->bams[j], opt->in_bam_fns[j]);
                    if (aux->idxs[j] == NULL) {
                        fai_destroy(aux->fai); sam_hdr_destroy(aux->headers[j]); sam_close(aux->bams[j]);
                        _err_error_exit("Failed to load index for \'%s\'\n", opt->in_bam_fns[j]);
                    }
                }
            }
        }
    }
    if (opt->out_aln_fp != NULL) {
        opt->bams = (samFile **)calloc(opt->n_in_bam_fn, sizeof(samFile *));
        opt->headers = (bam_hdr_t **)calloc(opt->n_in_bam_fn, sizeof(bam_hdr_t *));
        opt->idxs = (hts_idx_t **)calloc(opt->n_in_bam_fn, sizeof(hts_idx_t *));
        if (!opt->bams || !opt->headers || !opt->idxs) _err_error_exit("Memory allocation failure\n");
        for (int j = 0; j < opt->n_in_bam_fn; ++j) {
            opt->bams[j] = sam_open(opt->in_bam_fns[j], "r"); if (opt->bams[j] == NULL) _err_error_exit("Failed to open alignment file \'%s\'\n", opt->in_bam_fns[j]);
            const htsFormat *fmt = hts_get_format(opt->bams[j]); if (!fmt) _err_error_exit("Failed to get format of alignment file \'%s\'\n", opt->in_bam_fns[j]);
            if (fmt->format != bam && fmt->format != cram) {
                sam_close(opt->bams[j]); _err_error_exit("Input alignment file must be BAM or CRAM format.\n");
            }
            opt->headers[j] = sam_hdr_read(opt->bams[j]);
            if (opt->headers[j] == NULL) {
                sam_close(opt->bams[j]);
                _err_error_exit("Failed to read alignment file header \'%s\'\n", opt->in_bam_fns[j]);
            }
            if (fmt->format == cram) {
                if (hts_set_fai_filename(opt->bams[j], opt->ref_fa_fn) != 0) {
                    sam_hdr_destroy(opt->headers[j]); sam_close(opt->bams[j]);
                    _err_error_exit("Failed to set reference file for CRAM decoding: %s %s\n", opt->in_bam_fns[j], opt->ref_fa_fn);
                }
            }
            // hts_set_threads(opt->bams[j], opt->n_threads);
            opt->idxs[j] = sam_index_load(opt->bams[j], opt->in_bam_fns[j]);
            if (opt->idxs[j] == NULL) { // attempt to create index if not exist
                _err_warning("Index not found for \'%s\', creating now ...\n", opt->in_bam_fns[j]);
                if (sam_index_build(opt->in_bam_fns[j], 0) < 0) {
                    sam_hdr_destroy(opt->headers[j]); sam_close(opt->bams[j]);
                    _err_error_exit("Failed to build index for \'%s\'\n", opt->in_bam_fns[j]);
                } else {
                    opt->idxs[j] = sam_index_load(opt->bams[j], opt->in_bam_fns[j]);
                    if (opt->idxs[j] == NULL) {
                        sam_hdr_destroy(opt->headers[j]); sam_close(opt->bams[j]);
                        _err_error_exit("Failed to load index for \'%s\'\n", opt->in_bam_fns[j]);
                    }
                }
            }
        }
    } else {
        opt->bams = NULL; opt->headers = NULL; opt->idxs = NULL;
    }
    // collect sample name (SM) from BAM header
    if (opt->sample_name == NULL) opt->sample_name = extract_sample_name_from_bam_header(pl->io_aux[0].headers[0]); // extract sample names from first BAM header
    if (opt->sample_name == NULL) {
        opt->sample_name = strdup(opt->in_bam_fns[0]);
        if (opt->n_in_bam_fn > 1) { // cat all input bam file names
            for (int i = 1; i < opt->n_in_bam_fn; ++i) {
                opt->sample_name = realloc(opt->sample_name, strlen(opt->sample_name) + strlen(opt->in_bam_fns[i]) + 2);
                sprintf(opt->sample_name + strlen(opt->sample_name), ",%s", opt->in_bam_fns[i]);
            }
        }
    }
    // Per-sample VCF column names. Each sample's name comes from its own
    // first input's @RG SM, not from the first input overall, or both columns
    // would be labelled with the tumour's name. A duplicate-sample VCF is
    // rejected by bcftools outright, so a collision is resolved rather than
    // emitted.
    if (opt->n_samples > 1) {
        opt->sample_names = (char**)calloc(opt->n_samples, sizeof(char*));
        for (int smp = 0; smp < opt->n_samples; ++smp) {
            if (opt->sample_names[smp] != NULL) continue;
            for (int j = 0; j < opt->n_in_bam_fn; ++j) {
                if (opt->bam_i_to_sample[j] != smp) continue;
                opt->sample_names[smp] = extract_sample_name_from_bam_header(pl->io_aux[0].headers[j]);
                break;
            }
            if (opt->sample_names[smp] == NULL)
                opt->sample_names[smp] = strdup(smp == LONGCALLD_SAMPLE_NORMAL ? "NORMAL" : "TUMOR");
        }
        if (opt->tumor_name != NULL) {
            free(opt->sample_names[LONGCALLD_SAMPLE_TUMOR]);
            opt->sample_names[LONGCALLD_SAMPLE_TUMOR] = strdup(opt->tumor_name);
        }
        if (opt->normal_name != NULL) {
            free(opt->sample_names[LONGCALLD_SAMPLE_NORMAL]);
            opt->sample_names[LONGCALLD_SAMPLE_NORMAL] = strdup(opt->normal_name);
        }
        if (strcmp(opt->sample_names[LONGCALLD_SAMPLE_TUMOR], opt->sample_names[LONGCALLD_SAMPLE_NORMAL]) == 0) {
            _err_warning("Tumour and normal report the same sample name '%s'; suffixing them _T and _N. "
                         "Use --tumor-name/--normal-name to choose your own.\n", opt->sample_names[LONGCALLD_SAMPLE_TUMOR]);
            for (int smp = 0; smp < 2; ++smp) {
                char *old_nm = opt->sample_names[smp];
                char *nm = (char*)malloc(strlen(old_nm) + 3);
                sprintf(nm, "%s_%c", old_nm, smp == LONGCALLD_SAMPLE_NORMAL ? 'N' : 'T');
                free(old_nm); opt->sample_names[smp] = nm;
            }
        }
        _err_info("VCF sample columns: %s (normal), %s (tumour)\n",
                  opt->sample_names[LONGCALLD_SAMPLE_NORMAL], opt->sample_names[LONGCALLD_SAMPLE_TUMOR]);
    }
    // collect regions
    if (collect_regions(pl, opt, n_regions, regions) <= 0) {
        if (n_regions > 0) _err_warning("No valid regions found in the alignment file. Processing the entire alignment file.\n");
        else _err_warning("No autosomes or sex chromosomes (chr{1-22,XY}) found in the alignment file. Processing the entire alignment file.\n");
        opt->only_autosome = 0, opt->only_autosome_XY=0;
        collect_regions(pl, opt, 0, NULL);
    }
}

// Is this @PG ID already used by some line in the header?
static int pg_id_is_taken(bam_hdr_t *header, int n_pg, const char *id) {
    for (int i = 0; i < n_pg; ++i) {
        const char *nm = sam_hdr_line_name(header, "PG", i);
        if (nm != NULL && strcmp(nm, id) == 0) return 1;
    }
    return 0;
}

// Add exactly one @PG record describing this run.
//
// sam_hdr_add_pg() appends one @PG per *leaf* of the header's existing PP
// chain, so an input header whose chain has 87 leaves gained 87 identical
// longcallD records (measured on the bundled ONT test data). Instead, chain a
// single record to the last @PG line present.
static void add_longcalld_pg_line(bam_hdr_t *header) {
    int n_pg = sam_hdr_count_lines(header, "PG");
    if (n_pg < 0) n_pg = 0;

    // Pick a collision-free ID: "longcallD", else "longcallD.1", ...
    char id[32];
    snprintf(id, sizeof(id), "%s", PROG);
    for (int i = 1; i < 100000 && pg_id_is_taken(header, n_pg, id); ++i)
        snprintf(id, sizeof(id), "%s.%d", PROG, i);

    // Keep CL bounded; an unbounded command line makes the header unpleasant
    // for downstream tools and some parsers cope badly with very long records.
    char cl_buf[4096]; const char *cl = CMD;
    if (cl != NULL && strlen(cl) >= sizeof(cl_buf)) {
        snprintf(cl_buf, sizeof(cl_buf), "%.*s...", (int)(sizeof(cl_buf) - 4), cl);
        cl = cl_buf;
    }

    const char *pp = n_pg > 0 ? sam_hdr_line_name(header, "PG", n_pg - 1) : NULL;
    int r;
    if (pp != NULL && cl != NULL)
        r = sam_hdr_add_line(header, "PG", "ID", id, "PN", PROG, "VN", LONGCALLD_VERSION, "CL", cl, "PP", pp, NULL);
    else if (pp != NULL)
        r = sam_hdr_add_line(header, "PG", "ID", id, "PN", PROG, "VN", LONGCALLD_VERSION, "PP", pp, NULL);
    else if (cl != NULL)
        r = sam_hdr_add_line(header, "PG", "ID", id, "PN", PROG, "VN", LONGCALLD_VERSION, "CL", cl, NULL);
    else
        r = sam_hdr_add_line(header, "PG", "ID", id, "PN", PROG, "VN", LONGCALLD_VERSION, NULL);
    if (r < 0) _err_error_exit("Fail to add @PG line to bam header.\n");
}

// Is the header partitioned so that every contig we skip sorts after every
// contig we process? If so, the skipped contigs' records can simply be appended
// after the pipeline's output and the file stays coordinate-sorted. GRCh38
// headers are laid out this way (chr1..chrY, then chrM, randoms, decoys), but
// nothing guarantees it.
static int passthrough_keeps_sort_order(call_var_opt_t *opt, bam_hdr_t *header) {
    int last_kept = -1;
    for (int i = 0; i < header->n_targets; ++i) {
        if (!skip_target_region(opt, header->target_name[i])) last_kept = i;
    }
    for (int i = 0; i < last_kept; ++i) {
        if (skip_target_region(opt, header->target_name[i])) return 0; // skipped contig sorts before a kept one
    }
    return 1;
}

// Stream the records the region pipeline never visits: those on contigs the
// target filter excludes, and unmapped reads. Without this the output is a
// region-restricted subset rather than a drop-in replacement for the input,
// which is not what "output the realigned reads" should mean.
//
// Our HP/PS are stripped: these reads were never phased, so carrying a tag over
// from the input would misattribute someone else's phasing to us.
static int64_t write_passthrough_records(call_var_opt_t *opt, bam_hdr_t *header) {
    int64_t n = 0;
    bam1_t *b = bam_init1();
    for (int i = 0; i < opt->n_in_bam_fn; ++i) {
        samFile *in_bam = opt->bams[i]; hts_idx_t *idx = opt->idxs[i];
        if (in_bam == NULL || idx == NULL) continue;
        int smp = opt->bam_i_to_sample != NULL ? opt->bam_i_to_sample[i] : 0;
        aln_out_t *out = (smp >= 0 && smp < LONGCALLD_MAX_SAMPLES) ? opt->aln_outs[smp] : NULL;
        if (out == NULL) continue; // this sample has no alignment output
        for (int t = 0; t < header->n_targets; ++t) {
            if (!skip_target_region(opt, header->target_name[t])) continue;
            hts_itr_t *iter = sam_itr_queryi(idx, t, 0, header->target_len[t]);
            if (iter == NULL) continue;
            while (sam_itr_next(in_bam, iter, b) >= 0) {
                uint8_t *hp = bam_aux_get(b, "HP"); if (hp != NULL) bam_aux_del(b, hp);
                uint8_t *ps = bam_aux_get(b, "PS"); if (ps != NULL) bam_aux_del(b, ps);
                aln_out_push(out, b, header->target_name[t]);
                n++;
            }
            hts_itr_destroy(iter);
        }
        // Unmapped reads with no coordinate; these sort after everything placed.
        hts_itr_t *iter = sam_itr_queryi(idx, HTS_IDX_NOCOOR, 0, 0);
        if (iter != NULL) {
            while (sam_itr_next(in_bam, iter, b) >= 0) {
                uint8_t *hp = bam_aux_get(b, "HP"); if (hp != NULL) bam_aux_del(b, hp);
                uint8_t *ps = bam_aux_get(b, "PS"); if (ps != NULL) bam_aux_del(b, ps);
                aln_out_push(out, b, "*");
                n++;
            }
            hts_itr_destroy(iter);
        }
    }
    bam_destroy1(b);
    return n;
}

// Open one alignment output: its own @PG, its own @HD sort order, its own
// reorder buffer. The header is shared between outputs, which is safe only
// because call_var_pl_open_fa_bam() has already rejected inputs whose contig
// lists disagree -- otherwise a record's tid would mean different things in the
// two files.
static struct aln_out_t *open_one_aln_out(call_var_opt_t *opt, htsFile *fp, int is_cram,
                                          bam_hdr_t *header, hts_pos_t window, int n_threads) {
    if (fp == NULL) return NULL;
    if (is_cram) {
        if (hts_set_fai_filename(fp, opt->ref_fa_fn) != 0)
            _err_error_exit("Failed to set reference file for output CRAM encoding: %s\n", opt->ref_fa_fn);
    }
    hts_set_threads(fp, n_threads);
    if (sam_hdr_write(fp, header) < 0) _err_error_exit("Failed to write BAM header.\n");
    return aln_out_init(fp, header, is_cram, opt->no_reorder, window);
}

static void call_var_pl_write_bam_header(call_var_opt_t *opt, bam_hdr_t *header) {
    add_longcalld_pg_line(header);

    // Say what we actually produce. With the reorder buffer the output really is
    // coordinate-sorted; with --no-reorder it is not, and claiming otherwise is
    // how --refine-aln used to emit a file that samtools index rejected.
    if (opt->out_aln_passthrough && !opt->no_reorder && !passthrough_keeps_sort_order(opt, header)) {
        _err_warning("This header interleaves excluded contigs with processed ones, so appending the "
                     "excluded records keeps every read but cannot keep the file sorted; declaring "
                     "SO:unsorted. Use --all-ctg to process everything, or "
                     "--no-out-aln-passthrough to drop those records instead.\n");
        opt->no_reorder = 1;
    }
    if (opt->no_reorder) {
        if (sam_hdr_update_hd(header, "SO", "unsorted", "GO", "none") < 0)
            _err_warning("Failed to set @HD SO:unsorted in the output header.\n");
    } else {
        if (sam_hdr_update_hd(header, "SO", "coordinate") < 0)
            _err_warning("Failed to set @HD SO:coordinate in the output header.\n");
    }
    // Realignment can only move a record's POS left, and only within the noisy
    // region it was realigned inside -- regions longer than max_noisy_reg_len are
    // skipped entirely -- so that length (plus both flanks) bounds the shift.
    // Measured worst case on the bundled data is far smaller: 617 bp (ONT),
    // 500 bp (HiFi). Without --refine-aln no POS changes at all, so the window
    // is zero and records stream straight through.
    hts_pos_t window = opt->refine_bam ? (hts_pos_t)opt->max_noisy_reg_len + 2 * opt->noisy_reg_flank_len : 0;
    // Thread budget is shared between the outputs.
    int n_out = (opt->out_aln_fp != NULL) + (opt->out_naln_fp != NULL);
    int per_out = n_out > 1 ? (opt->n_threads + 1) / 2 : opt->n_threads;
    opt->aln_outs[LONGCALLD_SAMPLE_TUMOR] =
        open_one_aln_out(opt, opt->out_aln_fp, opt->out_aln_is_cram, header, window, per_out);
    opt->aln_outs[LONGCALLD_SAMPLE_NORMAL] =
        open_one_aln_out(opt, opt->out_naln_fp, opt->out_naln_is_cram, header, window, per_out);
}

// work with sorted SAM/BAM/CRAM
static void *call_var_worker_pipeline(void *shared, int step, void *in) { // kt_pipeline() callback
    call_var_pl_t *pl = (call_var_pl_t*)shared;
    if (step == 0) { // step 0: read bam records into BAM chunks, call variants
        if (pl->reg_chunk_i >= pl->n_reg_chunks) return 0;
        if (LONGCALLD_VERBOSE >= 2) _err_info("Step 0: load BAM & call variants.\n");
        // _err_info("Processing contig(s): %d/%d ...\n", pl->reg_chunk_i+1, pl->n_reg_chunks);
        call_var_step_t *s = calloc(1, sizeof(call_var_step_t));
        s->pl = pl;
        s->n_chunks = pl->reg_chunks[pl->reg_chunk_i].n_regions;
        s->chunks = calloc(s->n_chunks, sizeof(bam_chunk_t));
        // _err_info("Processing ctg-id: %d, n_chunks: %d (%d/%d)\n", pl->reg_chunks[pl->reg_chunk_i].reg_tids[0], s->n_chunks, pl->reg_chunk_i+1, pl->n_reg_chunks);
        kt_for(pl->n_threads, collect_bam_call_var_worker_for, s, s->n_chunks);
        pl->reg_chunk_i++;
        return s;
    } else if (step == 1) { // step 1: stitch the results
                            //         merge variants from different chunks
                            //         1) merge overlapping variants
                            //         2) update phase set and haplotype (flip if needed)
        if (LONGCALLD_VERBOSE >= 2) _err_info("Step 2: stitch & make variants\n");
        if (((call_var_step_t*)in)->n_chunks > 0) {
            ((call_var_step_t*)in)->vars = calloc(((call_var_step_t*)in)->n_chunks, sizeof(var_t));
            stitch_var_main((call_var_step_t*)in, ((call_var_step_t*)in)->chunks);
            kt_for(pl->n_threads, make_var_worker_for, in, ((call_var_step_t*)in)->n_chunks);
        }
        int64_t n_processed_reads = 0;
        for (int i = 0; i < ((call_var_step_t*)in)->n_chunks; ++i) {
            int n_up_ovlp_reads = 0;
            for (int j = 0; j < pl->opt->n_in_bam_fn; ++j) n_up_ovlp_reads += ((call_var_step_t*)in)->chunks[i].n_up_ovlp_reads[j]; 
            n_processed_reads += (((call_var_step_t*)in)->chunks[i].n_reads - n_up_ovlp_reads);
        }
        if (n_processed_reads > 0) _err_info("Processed %" PRIi64 " reads, %d/%d chunks\n", n_processed_reads, pl->reg_chunk_i, pl->n_reg_chunks);
        return in;
    } else if (step == 2) { // step 3: write the buffer to output
        if (LONGCALLD_VERBOSE >= 2) _err_info("Step 3: output variants (& phased bam)\n");
        call_var_step_t *s = (call_var_step_t*)in;
        int n_out_vars = 0, n_out_reads = 0;
        for (int i = 0; i < s->n_chunks; ++i) {
            bam_chunk_t *c = s->chunks + i;
            if (c->is_invalid) { var_free(s->vars + i); continue; }
            for (int smp = 0; smp < c->n_samples && smp < LONGCALLD_MAX_SAMPLES; ++smp)
                pl->n_reads_per_sample[smp] += c->n_sample_reads[smp];
            n_out_vars += write_var_to_vcf(s->vars+i, pl->opt, c);
            if (pl->opt->aln_outs[LONGCALLD_SAMPLE_TUMOR] != NULL || pl->opt->aln_outs[LONGCALLD_SAMPLE_NORMAL] != NULL)
                n_out_reads += write_read_to_bam(c, pl->opt, pl->io_aux);
            var_free(s->vars + i);  // free output
        }
        if (n_out_vars > 0) _err_info("Output %d variants to VCF\n", n_out_vars);
        if (n_out_reads > 0) {
            if (pl->opt->out_aln_is_cram) _err_info("Output %d reads to CRAM\n", n_out_reads);
            else _err_info("Output %d reads to BAM\n", n_out_reads);
        }
        bam_chunks_post_free(s->chunks, s->n_chunks, pl->opt); // free input
        free(s->vars); free(s);
    }
    return 0;
}

// 0 = `call`, 1 = `paired`. One option table and one usage function serve
// both; the subcommand only changes the synopsis line and which options are
// worth showing first.
static int LONGCALLD_CLI_PAIRED = 0;

static void call_var_usage(void) {//main usage
    fprintf(stderr, "\n");
    // fprintf(stderr, "Program: %s (%s)\n", PROG, DESCRIP);
    // fprintf(stderr, "Version: %s\tContact: %s\n\n", LONGCALLD_VERSION, CONTACT); 

    if (LONGCALLD_CLI_PAIRED)
        fprintf(stderr, "Usage: %s paired [options] <ref.fa> <tumor.bam/cram> <normal.bam/cram> [region ...] > paired.vcf\n", PROG);
    else
        fprintf(stderr, "Usage: %s call [options] <ref.fa> <input.bam/cram> [region ...] > phased.vcf\n", PROG);
    fprintf(stderr, "\n");
    fprintf(stderr, "Note: \'ref.fa\' should contain the same contig/chromosome names as \'input.bam/cram\'\n");
    fprintf(stderr, "      \'input.bam/cram' should be sorted and indexed\n");
    // fprintf(stderr, "\n");

    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  Input:\n");
    fprintf(stderr, "    --hifi                HiFi reads [default]\n");
    fprintf(stderr, "    --ont                 ONT reads \n");
    fprintf(stderr, "    --region-file   FILE  region bed file []\n");
    fprintf(stderr, "                          each line is a region, e.g., chr1         (whole chromosome)\n");
    fprintf(stderr, "                                                       chr1 100     (chr1:100-END)\n");
    fprintf(stderr, "                                                       chr1 100 200 (chr1:100-200)\n");
    fprintf(stderr, "    -L --input-is-list    input file is a list of BAM/CRAM files of the same sample\n");
    fprintf(stderr, "                          each line is a file path, e.g., longcallD call -L ref.fa bam_list.txt\n");
    fprintf(stderr, "    -X --extra-bam  FILE  extra input BAM/CRAM file(s) of the same sample for variant calling []\n");
    fprintf(stderr, "    -N --normal     FILE  matched normal BAM/CRAM file(s); the positional input is then\n");
    fprintf(stderr, "                          treated as the tumour. Can be used multiple times. []\n");
    fprintf(stderr, "       --min-normal-phase-dp INT\n");
    fprintf(stderr, "                          with -N, min. normal-sample depth for a variant to be used for\n");
    fprintf(stderr, "                          phasing [%d]\n", 4);
    fprintf(stderr, "       --min-normal-dp INT   min. normal depth to classify a variant at all [8]\n");
    fprintf(stderr, "       --min-tumor-dp  INT   min. tumour depth to classify a variant at all [8]\n");
    fprintf(stderr, "       --min-normal-hap-dp INT\n");
    fprintf(stderr, "                          min. depth on EACH normal haplotype for the four-slot test [4]\n");
    fprintf(stderr, "       --min-tumor-hap-dp INT\n");
    fprintf(stderr, "                          min. depth on the candidate tumour haplotype [4]\n");
    fprintf(stderr, "       --max-hap-alt   INT   alt reads tolerated in a slot that should be clean [0]\n");
    fprintf(stderr, "       --som-pval    FLOAT   max. Fisher exact p for a somatic call [1e-3]\n");
    fprintf(stderr, "                          can be used multiple times, e.g., -X bam1 -X bam2\n");
    fprintf(stderr, "    --autosome-XY         only call variants on autosomes and sex chromosomes, e.g., (chr){1-22,XY} [default]\n");
    // fprintf(stderr, "                          if not found, all chromosomes/contigs will be processed automatically\n");
    fprintf(stderr, "    --autosome            only call variants on autosomes, e.g., (chr){1-22}\n");
    fprintf(stderr, "    --all-ctg             call variants on all chromosomes/contigs\n");
    fprintf(stderr, "    -E --exclude-ctg STR  exclude contig/chromosome []\n");
    fprintf(stderr, "                          can be used multiple times, e.g., -E hs37d5 -E chrM\n");
    fprintf(stderr, "    -r --ref-idx    FILE  .fai index file for reference FASTA file, detect automaticaly if not provided []\n");
    fprintf(stderr, "    -T --trans-elem FILE  transposable element sequence FASTA file []\n");
    // fprintf(stderr, "    --STR           FILE  short tandem repeat annotation file []\n");
    // fprintf(stderr, "                          required 4 columns: chrom, start, end, motif\n");
    // fprintf(stderr, "                          optional 5th column: ID of STR; use chrom_start_end as ID if not provided\n"); 
    fprintf(stderr, "  Output:\n");
    fprintf(stderr, "    -n --sample-name STR  sample name in VCF output. RG/SM tag or BAM file name will be used if not provided []\n");
    fprintf(stderr, "    -o --out-vcf    FILE  output phased VCF file [stdout]\n");
    fprintf(stderr, "    -O --out-type    STR  v/z: un/compressed VCF [v]\n");
    fprintf(stderr, "    -l --min-sv-len  INT  min. length to be considered as SV [%d]\n", LONGCALLD_MIN_SV_LEN);
    fprintf(stderr, "                          SV-related information will be added to INFO field, e.g., SVLEN/SVTYPE/TSD\n");
    fprintf(stderr, "    -H --no-vcf-header    suppress the header in VCF output\n");
    fprintf(stderr, "    -s --mosaic/somatic   output low allele-frequency mosaic/somatic variants\n");
    // fprintf(stderr, "    -m --methylation    output methylation site information [False]\n");
    // fprintf(stderr, "                          if present, MM/ML tags will be used to calculate methylation level\n");
    fprintf(stderr, "       --amb-base         output variant with ambiguous base (N)\n");
    fprintf(stderr, "    -S/b/C --out-sam/bam/cram  FILE\n");
    fprintf(stderr, "                          output phased SAM/BAM/CRAM file []\n");
    fprintf(stderr, "                          note: multiple input BAM/CRAM files will be merged in SAM/BAM/CRAM output\n");
    fprintf(stderr, "    --refine-aln          refine alignment in SAM/BAM/CRAM output\n");
    fprintf(stderr, "                          output stays coordinate-sorted and indexable\n");
    fprintf(stderr, "    --out-normal-bam/-cram/-sam FILE\n");
    fprintf(stderr, "                          with -N, write the matched normal's phased reads here;\n");
    fprintf(stderr, "                          -b/-C/-S then carry the tumour's []\n");
    fprintf(stderr, "    --no-out-aln-passthrough\n");
    fprintf(stderr, "                          omit reads on excluded contigs and unmapped reads from the\n");
    fprintf(stderr, "                          alignment output (they are included by default)\n");
    fprintf(stderr, "    --no-reorder          do not re-sort refined records; output is written in input\n");
    fprintf(stderr, "                          order and the header declares SO:unsorted\n");
    fprintf(stderr, "    --out-var-rnames      output names of supporting reads in VCF FORMAT field for all variants [False]\n");
    fprintf(stderr, "    --out-som-var-rnames  output names of supporting reads in VCF FORMAT field for somatic variants [False]\n");
    fprintf(stderr, "    --out-sv-rnames       output names of supporting reads in VCF FORMAT field for SVs [False]\n");
    // fprintf(stderr, "\n");
    fprintf(stderr, "  Variant calling:\n");
    fprintf(stderr, "    -c --min-cov     INT  min. total read coverage for candidate variant [%d]\n", LONGCALLD_MIN_CAND_DP);
    fprintf(stderr, "    -d --alt-cov     INT  min. alt. read coverage for candidate variant [%d]\n", LONGCALLD_MIN_ALT_DP);
    fprintf(stderr, "    -a --alt-ratio FLOAT  min. alt. read ratio for candidate variant [%.2f]\n", LONGCALLD_MIN_CAND_AF);
    fprintf(stderr, "    -M --min-mapq    INT  min. mapping quality score to be used [%d]\n", LONGCALLD_MIN_CAND_MQ);
    fprintf(stderr, "    -B --min-bq      INT  min. base quality score to be used [%d]\n", LONGCALLD_MIN_CAND_BQ);
    fprintf(stderr, "       --max-cov     INT  max. total read coverage for candidate variant [%d]\n", LONGCALLD_MAX_NOISY_REG_COV);
    // fprintf(stderr, "    -p --max-ploidy  INT  max. ploidy [%d]\n", LONGCALLD_DEF_PLOID);
    fprintf(stderr, "  Low allele-frequency mosaic/somatic variant calling: (effective when -s/--mosaic/--somatic is used)\n");
    // fprintf(stderr, "    --alpha          INT  alpha value of beta-binomial distribution [%d]\n", LONGCALLD_SOMATIC_BETA_ALPHA);
    // fprintf(stderr, "    --beta           INT  beta value of beta-binomial distribution [%d]\n", LONGCALLD_SOMATIC_BETA_BETA);
    fprintf(stderr, "    --som-alt        INT  min. alt. read coverage for somatic variant [%d]\n", LONGCALLD_MIN_SOMATIC_ALT_DP);
    fprintf(stderr, "    --som-mei-alt    INT  min. alt. read coverage for somatic MEI variant [%d]\n", LONGCALLD_MIN_SOMATIC_TE_ALT_DP);
    fprintf(stderr, "    --max-somvar INT,INT  max. number of mosaic/somatic variants allowed in a window (m,w) [%d,%d]\n", LONGCALLD_SOMATIC_WIN_MAX_VARS, LONGCALLD_SOMATIC_WIN);
    // fprintf(stderr, "  Variant calling in noisy regions:\n");
    // fprintf(stderr, "    -x --max-xgap    INT  max. number of allowed substitutions/gap-bases in a sliding window(-w/--win-size) [%d]\n", LONGCALLD_NOISY_REG_MAX_XGAPS);
    // fprintf(stderr, "                          window with more than -x subs/gap-bases will be considered as noisy region\n");
    // fprintf(stderr, "    -w --win-size    INT  window size for searching noisy region [%d]\n", LONGCALLD_NOISY_REG_HIFI_SLIDE_WIN);
    // fprintf(stderr, "    -D --merge-dis   INT  default max. distance to merge two potential noisy SV regions [%d]\n", LONGCALLD_NOISY_REG_MERGE_DIS);
    // fprintf(stderr, "    -j --noisy-rat FLOAT  min. ratio of noisy reads in a window to call a noisy region [%.2f]\n", LONGCALLD_NOISY_REG_RATIO);
    // fprintf(stderr, "    -f --noisy-flank INT  flanking mask window size for noisy region [%d]\n", LONGCALLD_DENSE_FLANK_WIN);
    // fprintf(stderr, "    -c --end-clip    INT  max. number of clipping bases on both ends [%d]\n", LONGCALLD_NOISY_END_CLIP);
    // fprintf(stderr, "                          end-clipping region with more than -c bases will be considered as noisy clipping region\n");
    // fprintf(stderr, "    -F --clip-flank  INT  flanking mask window size for noisy clipping region [%d]\n", LONGCALLD_NOISY_END_CLIP_WIN);
    // fprintf(stderr, "\n");
    // fprintf(stderr, "    -p --hap-read    INT  when haplotype is available, min. number of full-spanning reads for each haplotype in noisy region to call a variant [%d]\n", LONGCALLD_MIN_HAP_FULL_READS);
    // fprintf(stderr, "    -f --full-read   INT  when haplotype is not available, min. number of full-spanning reads in noisy region to call a variant [%d]\n", LONGCALLD_MIN_NO_HAP_FULL_READS);

    // fprintf(stderr, "\n");
    fprintf(stderr, "  General:\n");
    fprintf(stderr, "    -t --threads     INT  number of threads to use [%d]\n", MIN_OF_TWO(CALL_VAR_THREAD_N, get_num_processors()));
    // fprintf(stderr, "    -h --help             print this help usage\n");
    fprintf(stderr, "    -v --version          print version number\n");
    // fprintf(stderr, "    -V --verbose     INT  verbose level (0-2). 0: none, 1: information, 2: debug [0]\n");
    // fprintf(stderr, "\n");
    exit(1);
}

int call_var_main(int argc, char *argv[]) { return call_var_main_mode(argc, argv, 0); }
int paired_var_main(int argc, char *argv[]) { return call_var_main_mode(argc, argv, 1); }

int call_var_main_mode(int argc, char *argv[], int paired_cli) {
    LONGCALLD_CLI_PAIRED = paired_cli;
    // _err_cmd("%s\n", CMD);
    // 'L' takes no argument (--input-is-list); it was "L:" here while the option
    // table said no_argument, so `-L ref.fa bam_list.txt` ate ref.fa as -L's
    // operand. 'D', 'j' and 'f' are gone along with the options that used them.
    const char *opt_str = "r:T:E:X:N:o:O:ml:Hb:C:S:sc:d:M:B:a:n:x:w:Lp:g:t:hvV:";
    int c, op_idx; call_var_opt_t *opt = call_var_init_para(); char *s;
    double realtime0 = realtime();
    // first round of getopt_long() to parse preset options
    while ((c = getopt_long(argc, argv, opt_str, call_var_opt, &op_idx)) >= 0) {
        switch(c) {
            case 0: 
                if (strcmp(call_var_opt[op_idx].name, "hifi") == 0) set_hifi_opt(opt);
                else if (strcmp(call_var_opt[op_idx].name, "ont") == 0) set_ont_opt(opt);
                break;
        }
    }
    optind = 1; // reset optind for next getopt_long() call
    // second round of getopt_long() to parse other options
    while ((c = getopt_long(argc, argv, opt_str, call_var_opt, &op_idx)) >= 0) {
        switch(c) {
            case 'r': opt->ref_fa_fai_fn = strdup(optarg); break;
            // case 'b': cgp->var_block_size = atoi(optarg); break;
            case 'o': opt->out_vcf_fn = strdup(optarg); break; // = fopen(optarg, "w"); break;
            case 'O': if (strcmp(optarg, "v") == 0 || strcmp(optarg, "z") == 0) opt->out_vcf_type = optarg[0];
                      else _err_error_exit("\'-O/--out-type\' can only be \'v\' or \'z\'\n");
                      break;
            case 'H': opt->no_vcf_header = 1; break;
            case 'l': opt->min_sv_len = atoi(optarg); break;
            case 'T': opt->te_seq_fn = strdup(optarg); make_te_kmer_idx(opt); break;
            case 0: if (strcmp(call_var_opt[op_idx].name, "amb-base") == 0) opt->out_amb_base = 1; 
                    // else if (strcmp(call_var_opt[op_idx].name, "hifi") == 0) set_hifi_opt(opt);
                    // else if (strcmp(call_var_opt[op_idx].name, "ont") == 0) set_ont_opt(opt);
                    else if (strcmp(call_var_opt[op_idx].name, "region-file") == 0 || 
                             strcmp(call_var_opt[op_idx].name, "regions-file") == 0) opt->reg_bed_fn = strdup(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "all-ctg") == 0) opt->only_autosome = 0, opt->only_autosome_XY=0;
                    else if (strcmp(call_var_opt[op_idx].name, "autosome") == 0) opt->only_autosome = 1;
                    else if (strcmp(call_var_opt[op_idx].name, "autosome-XY") == 0) opt->only_autosome_XY = 1;
                    else if (strcmp(call_var_opt[op_idx].name, "mosaic") == 0) opt->out_somatic = 1;
                    // else if (strcmp(call_var_opt[op_idx].name, "alpha") == 0) opt->somatic_beta_alpha = atoi(optarg);
                    // else if (strcmp(call_var_opt[op_idx].name, "beta") == 0) opt->somatic_beta_beta = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "som-alt") == 0) opt->min_somatic_alt_dp = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "som-mei-alt") == 0) opt->min_somatic_te_dp = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "max-somvar") == 0) {
                        opt->somatic_win_max_vars = strtol(optarg, &s, 10); 
                        if (*s == ',') opt->somatic_win = strtol(s+1, &s, 10);
                    } else if (strcmp(call_var_opt[op_idx].name, "refine-aln") == 0) opt->refine_bam = 1;
                    else if (strcmp(call_var_opt[op_idx].name, "out-var-rnames") == 0) opt->output_var_rnames = 1;
                    else if (strcmp(call_var_opt[op_idx].name, "out-sv-rnames") == 0) opt->output_sv_rnames = 1;
                    else if (strcmp(call_var_opt[op_idx].name, "out-som-var-rnames") == 0) opt->output_somatic_var_rnames = 1;
                    else if (strcmp(call_var_opt[op_idx].name, "max-cov") == 0) opt->max_noisy_reg_cov = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "merge-dis") == 0) opt->noisy_reg_merge_dis = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "end-clip") == 0) opt->end_clip_reg = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "clip-flank") == 0) opt->end_clip_reg_flank_win = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "no-reorder") == 0) opt->no_reorder = 1;
                    else if (strcmp(call_var_opt[op_idx].name, "min-normal-phase-dp") == 0) opt->min_normal_phase_dp = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "min-normal-dp") == 0) opt->min_normal_dp = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "min-tumor-dp") == 0) opt->min_tumor_dp = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "min-normal-hap-dp") == 0) opt->min_normal_hap_dp = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "min-tumor-hap-dp") == 0) opt->min_tumor_hap_dp = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "max-hap-alt") == 0) opt->max_hap_alt = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "max-poa-cells") == 0) opt->max_poa_cells = strtoll(optarg, NULL, 10);
                    else if (strcmp(call_var_opt[op_idx].name, "som-pval") == 0) opt->som_pval = atof(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "no-out-aln-passthrough") == 0) opt->out_aln_passthrough = 0;
                    else if (strcmp(call_var_opt[op_idx].name, "out-normal-bam") == 0)  { opt->out_naln_fp = hts_open(optarg, "wb"); opt->out_naln_is_cram = 0; }
                    else if (strcmp(call_var_opt[op_idx].name, "out-normal-cram") == 0) { opt->out_naln_fp = hts_open(optarg, "wc"); opt->out_naln_is_cram = 1; }
                    else if (strcmp(call_var_opt[op_idx].name, "out-normal-sam") == 0)  { opt->out_naln_fp = hts_open(optarg, "w");  opt->out_naln_is_cram = 0; }
                    else if (strcmp(call_var_opt[op_idx].name, "tumor-name") == 0 ||
                             strcmp(call_var_opt[op_idx].name, "tumour-name") == 0) { free(opt->tumor_name); opt->tumor_name = strdup(optarg); }
                    else if (strcmp(call_var_opt[op_idx].name, "normal-name") == 0) { free(opt->normal_name); opt->normal_name = strdup(optarg); }
                    break;
            case 's': opt->out_somatic = 1; break;
            case 'm': opt->out_methylation = 1; break;
            case 'E': set_exclude_ctg(opt, optarg); break;
            case 'X': set_extra_input_bam(opt, optarg); break;
            case 'N': set_normal_input_bam(opt, optarg); break;
            case 'b': opt->out_aln_fp = hts_open(optarg, "wb"); opt->out_aln_is_cram = 0; break;
            case 'C': opt->out_aln_fp = hts_open(optarg, "wc"); opt->out_aln_is_cram = 1; break;
            case 'S': opt->out_aln_fp = hts_open(optarg, "w"); opt->out_aln_is_cram = 0; break;
            case 'c': opt->min_dp = atoi(optarg); break;
            case 'd': opt->min_alt_dp = atoi(optarg); break;
            case 'a': opt->min_af = atof(optarg); break;
            case 'M': opt->min_mq = atoi(optarg); break;
            case 'B': opt->min_bq = atoi(optarg); break;
            case 'n': opt->sample_name = strdup(optarg); break;
            case 'L': opt->input_is_list = 1; break;
            case 'x': opt->noisy_reg_max_xgaps = atoi(optarg); break;
            case 'w': opt->noisy_reg_slide_win = atoi(optarg); break;
            case 'p': opt->min_hap_full_reads = atoi(optarg); break;
            // case 'f': opt->min_no_hap_full_reads = atoi(optarg); break;
            // case 'D': opt->noisy_reg_merge_dis = atoi(optarg); break;
            // case 'j': opt->min_noisy_reg_ratio = atof(optarg); break;
            // case 'L': opt->max_noisy_reg_len = atoi(optarg); break;
            // case 'c': opt->end_clip_reg = atoi(optarg); break;
            // case 'F': opt->end_clip_reg_flank_win = atoi(optarg); break;
            case 'g': if (strcmp(optarg, "right") == 0 || strcmp(optarg, "r") == 0) opt->gap_aln = LONGCALLD_GAP_RIGHT_ALN;
                      else if (strcmp(optarg, "left") == 0 || strcmp(optarg, "l") == 0) opt->gap_aln = LONGCALLD_GAP_LEFT_ALN;
                      else _err_error_exit("\'-g/--gap-aln\' can only be \'left\'/\'l\' or \'right\'/\'r\'\n"); // call_var_usage();
                      break; // without this, --gap-aln fell through to -t and set n_threads = atoi("left") = 0
            case 't': opt->n_threads = atoi(optarg); break;
            case 'h': call_var_free_para(opt); call_var_usage();
            case 'v': fprintf(stdout, "%s\n", LONGCALLD_VERSION); call_var_free_para(opt); return 0;
            case 'V': LONGCALLD_VERBOSE = atoi(optarg); break;
            // getopt_long has already printed what went wrong; exit non-zero so a
            // mistyped option cannot look like a successful run that produced no VCF.
            default: call_var_free_para(opt); return 1;
        }
    }

    if (argc - optind < 2) {
        call_var_free_para(opt); call_var_usage();
    } else {
        opt->ref_fa_fn = retrieve_full_url(argv[optind++]);
        if (opt->input_is_list) set_bam_list_file(opt, argv[optind++]);
        else opt->in_bam_fns[0] = retrieve_full_url(argv[optind++]);
        // `paired` takes the matched normal as its third operand. -N still
        // works and wins, so the remaining operands stay available as regions.
        if (paired_cli && opt->n_normal_bam_fn == 0) {
            if (argc - optind < 1) {
                call_var_free_para(opt);
                _err_error("paired needs a matched normal: either a third operand or -N/--normal.\n");
                call_var_usage();
            }
            set_normal_input_bam(opt, argv[optind++]);
        }
    }
    finalise_sample_layout(opt);
    if (opt->out_naln_fp != NULL && opt->n_samples < 2)
        _err_error_exit("--out-normal-bam/-cram/-sam needs a matched normal; supply one with -N/--normal.\n");
    // check if hifi and ont are both set
    if (opt->is_pb_hifi && opt->is_ont) _err_error_exit("Cannot set both --hifi and --ont\n");
    if (opt->only_autosome && opt->only_autosome_XY) _err_error_exit("Cannot set both --autosome and --autosome-XY\n");
    // threads
    if (opt->n_threads <= 0) {
        opt->n_threads = 1; opt->pl_threads = 1;
    } else {
        opt->n_threads = MIN_OF_TWO(opt->n_threads, get_num_processors());
        opt->pl_threads = MIN_OF_TWO(2, opt->n_threads); // MIN_OF_TWO(opt->n_threads, CALL_VAR_PL_THREAD_N);
    }
    if (opt->out_vcf_fn == NULL) opt->out_vcf_fn = strdup("-");
    if (opt->out_vcf_type == 'z') opt->out_vcf = hts_open(opt->out_vcf_fn, "wz");
    else opt->out_vcf = hts_open(opt->out_vcf_fn, "w");
    // set up pipeline for multi-threading
    call_var_pl_t pl;
    memset(&pl, 0, sizeof(call_var_pl_t));
    pl.max_reg_len_per_chunk = LONGCALLD_BAM_CHUNK_REG_SIZE; // pl.max_reads_per_chunk = LONGCALLD_BAM_CHUNK_READ_COUNT; 
    pl.n_threads = opt->n_threads; pl.min_reg_chunks_per_run = opt->n_threads * 4;
    // open BAM file & reference genome
    call_var_pl_open_fa_bam(opt, &pl, argv+optind, argc-optind);

    // write VCF/BAM header
    if (opt->out_aln_fp != NULL || opt->out_naln_fp != NULL) call_var_pl_write_bam_header(opt, pl.io_aux[0].headers[0]);
    if (opt->out_vcf != NULL && (opt->out_vcf_type == 'z' || opt->no_vcf_header == 0)) write_vcf_header(pl.io_aux[0].headers[0], opt);
    // start to work !!!
    pl.opt = opt;
    memset(pl.n_reads_per_sample, 0, sizeof(pl.n_reads_per_sample));
    kt_pipeline(opt->pl_threads, call_var_worker_pipeline, &pl, 3);
    if (opt->n_samples > 1) {
        _err_info("Reads used for calling: tumour %" PRId64 ", normal %" PRId64 "\n",
                  pl.n_reads_per_sample[LONGCALLD_SAMPLE_TUMOR],
                  pl.n_reads_per_sample[LONGCALLD_SAMPLE_NORMAL]);
    }
    int any_aln_out = 0;
    for (int i = 0; i < LONGCALLD_MAX_SAMPLES; ++i) if (opt->aln_outs[i] != NULL) any_aln_out = 1;
    if (any_aln_out && opt->out_aln_passthrough) {
        int64_t n_pass = write_passthrough_records(opt, pl.io_aux[0].headers[0]);
        if (n_pass > 0) _err_info("Passed through %" PRId64 " reads on excluded contigs or without a mapping\n", n_pass);
    }
    int64_t n_unrefined = 0;
    for (int i = 0; i < LONGCALLD_MAX_SAMPLES; ++i) {
        if (opt->aln_outs[i] == NULL) continue;
        n_unrefined += opt->aln_outs[i]->n_unrefined_supp;
        aln_out_free(opt->aln_outs[i]); opt->aln_outs[i] = NULL; // flush before closing
    }
    if (opt->refine_bam && n_unrefined > 0)
        _err_warning("%" PRId64 " supplementary alignment(s) were written with their original CIGAR: "
                     "--refine-aln rebuilds a CIGAR from the primary alignment's digar, which is not valid "
                     "for a supplementary.\n", n_unrefined);
    if (opt->out_aln_fp != NULL) hts_close(opt->out_aln_fp);
    if (opt->out_naln_fp != NULL) hts_close(opt->out_naln_fp);
    if (opt->out_vcf != NULL) {
        if (opt->vcf_hdr != NULL) bcf_hdr_destroy(opt->vcf_hdr);
        hts_close(opt->out_vcf);
    }
    call_var_free_pl(pl); call_var_free_para(opt); 
    // finish
    _err_info("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.\n", realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    _err_success("%s\n", CMD);
    return 0;
}
