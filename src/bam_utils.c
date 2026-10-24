#include <ctype.h>
#include <string.h>
#include "bam_utils.h"
#include "sam_internal.h"
#include "utils.h"
#include "collect_var.h"
#include "call_var_main.h"
#include "sdust.h"
#include "khash.h"

// qname -> the phase of that read's primary alignment, so a supplementary can
// inherit it. KHASH_MAP_INIT_STR does not copy keys, and the keys here are
// chunk->read_names[] entries, which outlive the map.
typedef struct { int hap; hts_pos_t ps; } read_phase_t;
KHASH_MAP_INIT_STR(rphase, read_phase_t)

extern int LONGCALLD_VERBOSE;


static read_var_profile_t *init_read_var_profile_inner(int n_reads, int *read_ids, int n_total_vars) {
    size_t n_profile_bytes = (size_t)n_reads * sizeof(read_var_profile_t);
    size_t n_var_values = (size_t)n_reads * (size_t)n_total_vars;
    size_t n_value_bytes = n_var_values * sizeof(int);
    uint8_t *raw = (uint8_t*)malloc(n_profile_bytes + n_value_bytes * 2);
    read_var_profile_t *p = (read_var_profile_t*)raw;
    int *alleles = NULL, *alt_qi = NULL;

    if (p == NULL) return NULL;

    alleles = (int*)(raw + n_profile_bytes);
    alt_qi = alleles + n_var_values;
    memset(alleles, 0xFF, n_value_bytes);
    memset(alt_qi, 0xFF, n_value_bytes);

    for (int i = 0; i < n_reads; ++i) {
        p[i].read_id = (read_ids == NULL ? i : read_ids[i]);
        p[i].start_var_idx = -1;
        p[i].end_var_idx = -2;
        p[i].alleles = alleles + (size_t)i * (size_t)n_total_vars;
        p[i].alt_qi = alt_qi + (size_t)i * (size_t)n_total_vars;
    }
    return p;
}

read_var_profile_t *init_read_var_profile(int n_reads, int n_total_vars) {
    return init_read_var_profile_inner(n_reads, NULL, n_total_vars);
}

read_var_profile_t *init_read_var_profile_with_ids(int n_reads, int *read_ids, int n_total_vars) {
    return init_read_var_profile_inner(n_reads, read_ids, n_total_vars);
}

void free_read_var_profile(read_var_profile_t *p, int n_reads) {
    (void)n_reads;
    free(p);
}

int has_equal_X_in_bam_cigar(bam1_t *read) {
    // Get the CIGAR string for the read
    const uint32_t *cigar = bam_get_cigar(read);
    int n_cigar = read->core.n_cigar;
    if (n_cigar <= 0) return 0;
    // Check if '=' or 'X' is in the CIGAR operations
    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            return 1;
        } else if (op == BAM_CMATCH) {
            return 0;
        }
    }
    return 0;
}

int has_cs_in_bam(bam1_t *b) {
    uint8_t *s = bam_aux_get(b, "cs");
    return s != NULL;
}

int has_MD_in_bam(bam1_t *b) {
    uint8_t *s = bam_aux_get(b, "MD");
    return s != NULL;
}

int get_aux_int_from_bam(bam1_t *b, const char *tag) {
    uint8_t *s = bam_aux_get(b, tag);
    if (s == NULL) return -1;
    return bam_aux2i(s);
}

char *get_aux_str_from_bam(bam1_t *b, const char *tag) {
    uint8_t *s = bam_aux_get(b, tag);
    if (s == NULL) return NULL;
    return (char*)s;
}

static void longcalld_copy_digar_read_buffers(bam_chunk_t *chunk, bam1_t *read, digar_t *digar) {
    uint32_t qlen = read->core.l_qseq;
    size_t packed_bseq_bytes = ((size_t)qlen + 1) / 2;

    digar->qlen = qlen;
    digar->bseq = packed_bseq_bytes > 0 ? (uint8_t *)malloc(packed_bseq_bytes) : NULL;
    if (packed_bseq_bytes > 0) memcpy(digar->bseq, bam_get_seq(read), packed_bseq_bytes);

    digar->qual = qlen > 0 ? (uint8_t *)malloc((size_t)qlen) : NULL;
    if (qlen > 0) memcpy(digar->qual, bam_get_qual(read), (size_t)qlen);
    for (uint32_t i = 0; i < qlen; ++i) {
        chunk->qual_counts[digar->qual[i]]++;
    }
}

void print_digar(digar_t *digar, FILE *fp) {
    fprintf(fp, "pos\ttype\tlen\tqi\tis_low_qual\n");
    for (int i = 0; i < digar->n_digar; ++i) {
        fprintf(fp, "%" PRId64 "\t%c\t%d\t%d\t%d\n", digar->digars[i].pos, BAM_CIGAR_STR[digar->digars[i].type], digar->digars[i].len, digar->digars[i].qi, digar->digars[i].is_low_qual);
    }
}

static inline uint32_t *longcalld_push_cigar(int *n_cigar, int *m_cigar, uint32_t *cigar, uint32_t op, int len) {
	if (*n_cigar == 0 || op != (cigar[(*n_cigar) - 1]&0xf)) {
		if (*n_cigar == *m_cigar) {
			*m_cigar = *m_cigar? (*m_cigar)<<1 : 4;
			cigar = (uint32_t*)realloc(cigar, (*m_cigar) << 2);
		}
		cigar[(*n_cigar)++] = len<<4 | op;
	} else cigar[(*n_cigar)-1] += len<<4;
	return cigar;
}

typedef struct {
    hts_pos_t *pos; int *lens;
    int *counts, *is_dense;
    int front, rear, count;
    int max_s, win; // max_s: max sub/gaps in win
    int size_m;
} xid_queue_t;

xid_queue_t *init_xid_queue(int max_sites, int max_s, int win) {
    xid_queue_t *q = (xid_queue_t*)malloc(sizeof(xid_queue_t));
    q->pos = (hts_pos_t*)malloc(max_sites * sizeof(hts_pos_t));
    q->counts = (int*)malloc(max_sites * sizeof(int));
    q->lens = (int*)malloc(max_sites * sizeof(int));
    q->is_dense = (int*)calloc(max_sites, sizeof(int));
    q->size_m = max_sites;
    q->front = 0; q->rear = -1; q->count = 0;
    q->max_s = max_s; q->win = win;
    return q;
}

void realloc_xid_queue(xid_queue_t *q) {
    if (q->rear+1 >= q->size_m) {
        q->size_m *= 2;
        q->pos = (hts_pos_t*)realloc(q->pos, q->size_m * sizeof(hts_pos_t));
        q->counts = (int*)realloc(q->counts, q->size_m * sizeof(int));
        q->lens = (int*)realloc(q->lens, q->size_m * sizeof(int));
        q->is_dense = (int*)realloc(q->is_dense, q->size_m * sizeof(int));
        for (int i = q->rear+1; i < q->size_m; ++i) q->is_dense[i] = 0;
    }
}

void free_xid_queue(xid_queue_t *q) {
    free(q->is_dense);
    free(q->pos); free(q->counts); 
    free(q->lens);
    free(q);
}

void push_xid_size_queue_win(xid_queue_t *q, hts_pos_t pos, int len, int count,
                             cgranges_t *cr, hts_pos_t *cr_cur_start, hts_pos_t *cr_cur_end,
                             int *cr_q_start, int *cr_q_end) {
    realloc_xid_queue(q);
    q->pos[++q->rear] = pos; // left-most pos of the region
    q->lens[q->rear] = len; q->counts[q->rear] = count;
    q->count += count;

    hts_pos_t noisy_start = -1, noisy_end = -1;
    while (q->pos[q->front]+q->lens[q->front]-1 <= pos - q->win) {
        q->count -= q->counts[q->front]; // q->lens[q->front];
        q->front++;
    }
    // fprintf(stderr, "pos: %" PRIi64 ", len: %d, count: %d, front: %d, rear: %d, count: %d\n", pos, len, count, q->front, q->rear, q->count);
    if (count > 0) {
        if (q->count > q->max_s) {
            for (int j = q->front; j <= q->rear; j++)
                q->is_dense[j] = 1; // mark all dense regions
            noisy_start = q->pos[q->front]; noisy_end = q->pos[q->rear]+q->lens[q->rear];
            if (*cr_cur_start == -1) {
                *cr_cur_start = noisy_start; *cr_cur_end = noisy_end;
                *cr_q_start = q->front; *cr_q_end = q->rear;
            } else {
                if (noisy_start <= *cr_cur_end) { // merge
                    *cr_cur_end = noisy_end; // *cr_qi stays the same
                    *cr_q_end = q->rear; // update end
                } else {
                    // add the previous region
                    int var_size = 0;
                    for (int i = *cr_q_start; i <= *cr_q_end; i++) var_size += q->counts[i];
                    var_size = var_size > (*cr_cur_end - *cr_cur_start + 1) ? var_size : (*cr_cur_end - *cr_cur_start + 1);
                    cr_add(cr, "cr", *cr_cur_start-1, *cr_cur_end, var_size);
                    *cr_cur_start = noisy_start; *cr_cur_end = noisy_end;
                    *cr_q_start = q->front; *cr_q_end = q->rear;
                }
            }
            // fprintf(stderr, "noisy region: %" PRId64 "-%" PRId64 "\n", noisy_start, noisy_end);
        }
    }
}

int get_var_start(cand_var_t *var_sites, int cur_site_i, int n_total_pos, hts_pos_t start) {
    hts_pos_t target = start > 0 ? start - 1 : start;
    int left = cur_site_i, right = n_total_pos;
    while (left < right) {
        int mid = left + (right - left) / 2;
        hts_pos_t mid_pos = var_sites[mid].var_type == BAM_CDIFF ? var_sites[mid].pos : var_sites[mid].pos - 1;
        if (mid_pos < target) left = mid + 1;
        else right = mid;
    }
    while (left < n_total_pos && var_sites[left].pos < start) left++;
    return left;
}

int get_var_site_start(var_site_t *var_sites, int cur_site_i, int n_total_pos, hts_pos_t start) {
    hts_pos_t target = start > 0 ? start - 1 : start;
    int left = cur_site_i, right = n_total_pos;
    while (left < right) {
        int mid = left + (right - left) / 2;
        hts_pos_t mid_pos = var_sites[mid].var_type == BAM_CDIFF ? var_sites[mid].pos : var_sites[mid].pos - 1;
        if (mid_pos < target) left = mid + 1;
        else right = mid;
    }
    while (left < n_total_pos && var_sites[left].pos < start) left++;
    return left;
}

// only for mismatch and insertion
uint8_t *make_alt_seq(digar1_t *digar, const uint8_t *bseq) {
    if (digar->type == BAM_CDEL) return NULL;
    uint8_t *alt_seq = (uint8_t*)malloc((digar->len) * sizeof(uint8_t));
    for (int i = 0; i < digar->len; ++i) {
        alt_seq[i] = seq_nt16_int[bam_seqi(bseq, digar->qi+i)];
    }
    return alt_seq;
}

void update_var_site_with_allele(cand_var_t *cand_var, int is_low_qual, uint8_t strand, int allele_i) {
    if (is_low_qual) { // low qual base quality, allele_i considered as -1
        cand_var->low_qual_cov++;
        return;
    }
    cand_var->total_cov++;
    cand_var->alle_covs[allele_i] += 1; // ref/alt
    cand_var->strand_to_alle_covs[strand][allele_i] += 1; // strand-wise: 0:forward/1:reverse -> allele_i -> read count
}

void update_read_var_profile_with_allele(int var_i, int allele_i, int alt_qi, read_var_profile_t *read_var_profile) {
    // even allele_i is -1/-2, we still set it XXX
    if (read_var_profile->start_var_idx == -1) read_var_profile->start_var_idx = var_i;
    read_var_profile->end_var_idx = var_i;
    int _var_i = var_i - read_var_profile->start_var_idx;
    read_var_profile->alleles[_var_i] = allele_i; // ref:0, alt:1~n, or -1: not ref/alt
    read_var_profile->alt_qi[_var_i] = alt_qi; // position of alt base in read, -1 if not alt
}


int get_digar_ave_qual(digar1_t *digar, const uint8_t *qual) {
    if (digar->is_low_qual) return 0;
    if (digar->qi < 0) return 0;
    int q_start, q_end;
    if (digar->type == BAM_CDEL) {
        if (digar->qi == 0) {
            q_start = 0;
            q_end = 0;
        } else {
            q_start = digar->qi - 1;
            q_end = digar->qi;
        }
    } else {
        q_start = digar->qi;
        q_end = digar->qi + digar->len - 1;
    }
    int ave_qual = 0, qlen = q_end - q_start + 1;
    for (int i = q_start; i <= q_end; ++i) {
        ave_qual += qual[i];
    }
    return ave_qual / qlen;
}

// cand_vars:
// X: 1 -> M SNPs
// I: 0 -> M INSs
// D: M -> 1 DEL
// XXX no minor_alt_allele here, will be used in read_var_profile
// here, minor_alt_allele will be simply treated as non-alt, i.e., ref allele
int update_cand_vars_from_digar(const call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_var_sites, var_site_t *var_sites, cand_var_t *cand_vars) {
    int site_i, digar_i = 0;
    // hts_pos_t pos_start = read->core.pos+1, pos_end = bam_endpos(read);
    hts_pos_t pos_start = digar->beg, pos_end = digar->end;
    digar1_t *digars = digar->digars; uint8_t strand = digar->is_rev;

    // int cur_start_i = site_i = get_var_site_start(var_sites, start_i, n_var_sites, pos_start);
    site_i = get_var_site_start(var_sites, 0, n_var_sites, pos_start);
    int tid = chunk->tid; // read->core.tid;
    // compare sorted var_sites and digars
    // both var_sites and digars are sorted by pos/type/ref_len
    for (; site_i < n_var_sites && digar_i < digar->n_digar; ) {
        if (digars[digar_i].type == BAM_CEQUAL) {
            digar_i++; continue;
        }
        var_site_t digar_var_site = make_var_site_from_digar(tid, digar->digars+digar_i);
        int ave_qual = get_digar_ave_qual(digar->digars+digar_i, digar->qual);
        // allow in-exact match for large INSs
        int ret = exact_comp_var_site_ins(opt, var_sites+site_i, &digar_var_site); // this is consistent with <collect_all_cand_var_sites> in collect_var.c
        // int ret, is_ovlp;
        // ret = fuzzy_comp_ovlp_var_site(opt, var_sites+site_i, &digar_var_site, &is_ovlp);
        if (ret < 0) { // var_site < digar_var_site
            // update_var_site_with_allele(cand_vars+site_i, is_in_noisy_reg(cand_vars[site_i].pos, digar->noisy_regs), strand, 0);
            update_var_site_with_allele(cand_vars+site_i, 0, strand, 0);
            site_i++;
        } else if (ret == 0) { // merge together, update/add alt allele
            update_var_site_with_allele(cand_vars+site_i, (digar->digars[digar_i].is_low_qual) || (ave_qual < opt->min_bq), strand, 1);
            // update_var_site_with_allele(cand_vars+site_i, (digar->digars[digar_i].is_low_qual) || (ave_qual < chunk->first_quar_qual), strand, 1);
            site_i++;
        } else { // ret > 0, var_site > digar_var_site
            digar_i++;
        }
    }
    for (; site_i < n_var_sites; ++site_i) {
        if (var_sites[site_i].pos > pos_end) break;
        // update_var_site_with_allele(cand_vars+site_i, is_in_noisy_reg(cand_vars[site_i].pos, digar->noisy_regs), strand, 0);
        update_var_site_with_allele(cand_vars+site_i, 0, strand, 0);
    }
    // return cur_start_i;
    return 0;
}

// int get_alt_allele_idx(cand_var_t *cand_var, digar1_t *digar, const uint8_t *bseq) {
    // if (digar->is_low_qual) return -2; // low qual base quality
    // else return 1; // alt allele
    // if (cand_var->var_type == BAM_CDEL) return 1; // only 0/1 for deletion
    // uint8_t *alt_seq = make_alt_seq(digar, bseq);
    // int alt_len = digar->len;
    // int allele_i = -1;
    // for (int i = 0; i < cand_var->n_uniq_alles-1; ++i) {
    //     if (cand_var->alt_lens[i] == alt_len && memcmp(cand_var->alt_seqs[i], alt_seq, alt_len) == 0) {
    //         allele_i = i+1;
    //         break;
    //     }
    // }
    // free(alt_seq);
    // return allele_i;
    // if (allele_i == -1) {
    //     // _err_error_exit("Alt allele not found: %d\n", allele_i);
    //     allele_i = cand_var->n_uniq_alles;
    // }
    // return allele_i;
// }

// update both cand_vars' total_cov/allele_cov and read_var_profile
int update_read_vs_somatic_var_profile_from_digar(const call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_cand_vars, cand_var_t *cand_vars, read_var_profile_t *read_var_profile) {
    int var_i, digar_i = 0, var_read_pos = -1, allele_i=-1;
    // hts_pos_t pos_start = read->core.pos+1, pos_end = bam_endpos(read);
    hts_pos_t pos_start = digar->beg, pos_end = digar->end;
    digar1_t *digar1 = digar->digars;
    var_i = get_var_start(cand_vars, 0, n_cand_vars, pos_start);
    int tid = chunk->tid; //read->core.tid;
    // compare sorted cand_vars and digars
    // both cand_vars and digars are sorted by pos/type/ref_len
    for (; var_i < n_cand_vars && digar_i < digar->n_digar; ) {
        if (digar1[digar_i].type == BAM_CEQUAL) {
            digar_i++; continue;
        }
        var_site_t var_site0 = make_var_site_from_cand_var(cand_vars+var_i);
        hts_pos_t var_pos_beg = var_site0.pos, var_pos_end = var_site0.pos + var_site0.ref_len - 1;
        var_site_t digar_var_site = make_var_site_from_digar(tid, digar1+digar_i);
        int ave_qual = get_digar_ave_qual(digar1+digar_i, digar->qual);
        var_read_pos = digar1[digar_i].qi; // read position of the digar
        int is_ovlp, ret;
        ret = fuzzy_comp_ovlp_var_site(opt, &var_site0, &digar_var_site, &is_ovlp);
        if (var_site0.var_type == BAM_CINS && var_site0.alt_len >= opt->min_sv_len) {
            // left side
            if (is_ovlp == 0 || ret != 0) { // check if nearby INS digar_var_site is a match
                int _is_ovlp, _ret, _digar_i;
                for (_digar_i = digar_i-1; _digar_i >= 0; --_digar_i) {
                    if (digar1[_digar_i].type != BAM_CINS || digar1[_digar_i].len < opt->min_sv_len) continue;
                    if (digar1[_digar_i].pos < var_site0.pos - var_site0.alt_len) break; // no overlap
                    var_site_t _digar_var_site = make_var_site_from_digar(tid, digar1+_digar_i);
                    _ret = fuzzy_comp_ovlp_var_site(opt, &var_site0, &_digar_var_site, &_is_ovlp);
                    if (_is_ovlp == 1 && _ret == 0) {
                        is_ovlp = 1, ret = 0; // we found a match
                        var_read_pos = digar1[_digar_i].qi; // update read position
                        break;
                    }
                }
            }
            // right side
            if (is_ovlp == 0 || ret != 0) { // check if nearby INS digar_var_site is a match
                int _is_ovlp, _ret, _digar_i;
                for (_digar_i = digar_i+1; _digar_i < digar->n_digar; ++_digar_i) {
                    if (digar1[_digar_i].type != BAM_CINS || digar1[_digar_i].len < opt->min_sv_len) continue;
                    if (digar1[_digar_i].pos > var_site0.pos + var_site0.alt_len) break; // no overlap
                    var_site_t _digar_var_site = make_var_site_from_digar(tid, digar1+_digar_i);
                    _ret = fuzzy_comp_ovlp_var_site(opt, &var_site0, &_digar_var_site, &_is_ovlp);
                    if (_is_ovlp == 1 && _ret == 0) {
                        is_ovlp = 1, ret = 0; // we found a match
                        var_read_pos = digar1[_digar_i].qi; // update read position
                        break;
                    }
                }
            }
        }
        if (is_ovlp == 0) { // no overlapping
            if (ret < 0) { // var_site < digar_var_site: ref_allele
                cand_vars[var_i].total_cov++;
                cand_vars[var_i].alle_covs[0] += 1; // ref allele
                update_read_var_profile_with_allele(var_i, 0, -1, read_var_profile);
                var_i++;
            } else if (ret > 0) { // var_site > digar_var_site
                digar_i++;
            } else { // var_site == digar_var_site
                _err_error("Unexpected case: is_ovlp == 0 && var_site == digar_var_site, %d-%c-%d-%d\n", cand_vars[var_i].pos, BAM_CIGAR_STR[cand_vars[var_i].var_type], digar_var_site.pos, BAM_CIGAR_STR[digar_var_site.var_type]);
                var_i++; digar_i++;
            }
        } else { // overlap
            if (ret == 0) { // exact the same with var_site
                allele_i = 1; // get_alt_allele_idx(cand_vars+var_i, digar->digars+digar_i, digar->bseq);
                if (ave_qual < opt->min_bq) allele_i = -2;
                // if (ave_qual < chunk->first_quar_qual) allele_i = -2;
                update_read_var_profile_with_allele(var_i, allele_i, var_read_pos, read_var_profile);
                cand_vars[var_i].total_cov++;
                if (allele_i >= 0) cand_vars[var_i].alle_covs[allele_i] += 1; 
                var_i++;
            } else {
                if (var_pos_beg >= pos_start && var_pos_end <= pos_end) {
                    cand_vars[var_i].alle_covs[0] += 1; 
                    cand_vars[var_i].total_cov++;
                }
                update_read_var_profile_with_allele(var_i, -1, -1, read_var_profile);
                var_i++;
            }
        }
    }
    for (; var_i < n_cand_vars; ++var_i) {
        if (cand_vars[var_i].pos > pos_end) break;
        cand_vars[var_i].total_cov++;
        cand_vars[var_i].alle_covs[0] += 1; // ref allele
        update_read_var_profile_with_allele(var_i, 0, -1, read_var_profile);
    }
    return 0;
}

// for germline and somatic: for SNP/small indels in clean regions, require exact match
// for somatic: for large indels, allow in-exact match
int update_read_vs_all_var_profile_from_digar(const call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_cand_vars, cand_var_t *cand_vars, int *var_i_to_cate, read_var_profile_t *read_var_profile) {
    int var_i, digar_i = 0, var_read_pos = -1, allele_i=-1;
    hts_pos_t pos_start = digar->beg, pos_end = digar->end;
    digar1_t *digar1 = digar->digars;
    // int cur_start_i = var_i = get_var_start(cand_vars, start_var_i, n_cand_vars, pos_start);
    var_i = get_var_start(cand_vars, 0, n_cand_vars, pos_start);
    int tid = chunk->tid; //read->core.tid;
    // compare sorted cand_vars and digars
    // both cand_vars and digars are sorted by pos/type/ref_len
    for (; var_i < n_cand_vars && digar_i < digar->n_digar; ) {
        if (var_i_to_cate[var_i] == LONGCALLD_NON_VAR) {
            var_i++; continue;
        }
        if (digar1[digar_i].type == BAM_CEQUAL) {
            digar_i++; continue;
        }
        var_site_t var_site0 = make_var_site_from_cand_var(cand_vars+var_i);
        hts_pos_t var_pos_beg = var_site0.pos, var_pos_end = var_site0.pos + var_site0.ref_len - 1;
        var_site_t digar_var_site = make_var_site_from_digar(tid, digar1+digar_i);
        int ave_qual = get_digar_ave_qual(digar1+digar_i, digar->qual);
        var_read_pos = digar1[digar_i].qi; // read position of the digar
        int is_ovlp, ret;
        if (var_i_to_cate[var_i] == LONGCALLD_CAND_SOMATIC_VAR) {
            ret = fuzzy_comp_ovlp_var_site(opt, &var_site0, &digar_var_site, &is_ovlp);
            if (var_site0.var_type == BAM_CINS && var_site0.alt_len >= opt->min_sv_len) {
                // left side
                if (is_ovlp == 0 || ret != 0) { // check if nearby INS digar_var_site is a match
                    int _is_ovlp, _ret, _digar_i;
                    for (_digar_i = digar_i-1; _digar_i >= 0; --_digar_i) {
                        if (digar1[_digar_i].type != BAM_CINS || digar1[_digar_i].len < opt->min_sv_len) continue;
                        if (digar1[_digar_i].pos < var_site0.pos - var_site0.alt_len) break; // no overlap
                        var_site_t _digar_var_site = make_var_site_from_digar(tid, digar1+_digar_i);
                        _ret = fuzzy_comp_ovlp_var_site(opt, &var_site0, &_digar_var_site, &_is_ovlp);
                        if (_is_ovlp == 1 && _ret == 0) {
                            is_ovlp = 1, ret = 0; // we found a match
                            var_read_pos = digar1[_digar_i].qi; // update read position
                            break;
                        }
                    }
                }
                // right side
                if (is_ovlp == 0 || ret != 0) { // check if nearby INS digar_var_site is a match
                    int _is_ovlp, _ret, _digar_i;
                    for (_digar_i = digar_i+1; _digar_i < digar->n_digar; ++_digar_i) {
                        if (digar1[_digar_i].type != BAM_CINS || digar1[_digar_i].len < opt->min_sv_len) continue;
                        if (digar1[_digar_i].pos > var_site0.pos + var_site0.alt_len) break; // no overlap
                        var_site_t _digar_var_site = make_var_site_from_digar(tid, digar1+_digar_i);
                        _ret = fuzzy_comp_ovlp_var_site(opt, &var_site0, &_digar_var_site, &_is_ovlp);
                        if (_is_ovlp == 1 && _ret == 0) {
                            is_ovlp = 1, ret = 0; // we found a match
                            var_read_pos = digar1[_digar_i].qi; // update read position
                            break;
                        }
                    }
                }
            }
        } else ret = comp_ovlp_var_site(opt, &var_site0, &digar_var_site, &is_ovlp);
        if (is_ovlp == 0) { // no overlapping
            if (ret < 0) { // var_site < digar_var_site: ref_allele
                if (var_i_to_cate[var_i] == LONGCALLD_CAND_SOMATIC_VAR) {
                    cand_vars[var_i].total_cov++;
                    cand_vars[var_i].alle_covs[0] += 1; // ref allele
                }
                update_read_var_profile_with_allele(var_i, 0, -1, read_var_profile);
                var_i++;
            } else if (ret > 0) { // var_site > digar_var_site
                digar_i++;
            } else { // var_site == digar_var_site
                _err_error("Unexpected case: is_ovlp == 0 && var_site == digar_var_site, %d-%c-%d-%d\n", cand_vars[var_i].pos, BAM_CIGAR_STR[cand_vars[var_i].var_type], digar_var_site.pos, BAM_CIGAR_STR[digar_var_site.var_type]);
                var_i++; digar_i++;
            }
        } else { // overlap
            if (ret == 0) { // exact the same with var_site
                allele_i = 1; // get_alt_allele_idx(cand_vars+var_i, digar->digars+digar_i, digar->bseq);
                if (ave_qual < opt->min_bq) allele_i = -2;
                // if (ave_qual < chunk->first_quar_qual) allele_i = -2;
                update_read_var_profile_with_allele(var_i, allele_i, var_read_pos, read_var_profile);
                if (var_i_to_cate[var_i] == LONGCALLD_CAND_SOMATIC_VAR) {
                    cand_vars[var_i].total_cov++;
                    if (allele_i >= 0) cand_vars[var_i].alle_covs[allele_i] += 1; 
                }
                var_i++;
            } else {
                if (var_i_to_cate[var_i] == LONGCALLD_CAND_SOMATIC_VAR) {
                    if (var_pos_beg >= pos_start && var_pos_end <= pos_end) {
                        cand_vars[var_i].alle_covs[0] += 1; 
                        cand_vars[var_i].total_cov++;
                    }
                }
                update_read_var_profile_with_allele(var_i, -1, -1, read_var_profile);
                var_i++;
            }
        }
    }
    for (; var_i < n_cand_vars; ++var_i) {
        if (cand_vars[var_i].pos > pos_end) break;
        if (is_in_noisy_reg(cand_vars[var_i].pos, digar->noisy_regs)) continue;
        if (var_i_to_cate[var_i] == LONGCALLD_CAND_SOMATIC_VAR) {
            cand_vars[var_i].total_cov++;
            cand_vars[var_i].alle_covs[0] += 1; // ref allele
        }
        update_read_var_profile_with_allele(var_i, 0, -1, read_var_profile);
    }
    // return cur_start_i;
    return 0;
}

void set_digar(digar1_t *digar, hts_pos_t pos, int type, int len, int qi, int is_low_qual, uint8_t *alt_seq) {
    digar->pos = pos; digar->type = type; digar->len = len; digar->qi = qi; digar->alt_seq = alt_seq; digar->is_low_qual = is_low_qual;
}

int same_digar1(digar1_t d1, digar1_t d2) {
    if (d1.type != BAM_CEQUAL && d1.type != BAM_CINS && d1.type != BAM_CDEL) return 0; // not ins or del
    if (d1.type == d2.type && d1.is_low_qual == d2.is_low_qual) {
        return 1; // same
    } return 0; // diff
}

// allocate new memory for alt_seq, instead of using the pointer directly
digar1_t *push_digar_alt_seq(digar1_t *digar, int *n_digar, int *m_digar, digar1_t d) {
    if (d.len <= 0) return digar; // no need to push
    if (*n_digar == 0 || same_digar1(digar[*n_digar-1], d) == 0) { // not the same as the last digar
        if (*n_digar == *m_digar) {
            *m_digar = *m_digar ? (*m_digar << 1) : 4;
            digar = (digar1_t*)realloc(digar, *m_digar * sizeof(digar1_t));
        }
        digar[*n_digar].pos = d.pos; digar[*n_digar].type = d.type;
        digar[*n_digar].len = d.len; digar[*n_digar].alt_seq = NULL;
        if (d.type == BAM_CDIFF || d.type == BAM_CINS) {
            digar[*n_digar].alt_seq = (uint8_t*)malloc(d.len * sizeof(uint8_t));
            for (int i = 0; i < d.len; ++i) digar[*n_digar].alt_seq[i] =  d.alt_seq[i];
        }
        digar[*n_digar].qi = d.qi; digar[*n_digar].is_low_qual = d.is_low_qual;
        (*n_digar)++;
    } else {
        if (digar[*n_digar-1].type == BAM_CINS) {
            uint8_t *alt_seq = (uint8_t*)realloc(digar[*n_digar-1].alt_seq, (digar[*n_digar-1].len + d.len) * sizeof(uint8_t));
            for (int i = digar[*n_digar-1].len; i < digar[*n_digar-1].len + d.len; ++i)
                alt_seq[i] = d.alt_seq[i - digar[*n_digar-1].len];
            digar[*n_digar-1].alt_seq = alt_seq;
        }
        digar[*n_digar-1].len += d.len; // merge INS/DEL with the last digar
    }
    return digar; // return the new digar pointer
}

digar1_t *push_digar0(digar1_t *digar, int *n_digar, int *m_digar, digar1_t d) {
    if (d.len <= 0) return digar;
    if (*n_digar == 0 || same_digar1(digar[*n_digar-1], d) == 0) { // not the same as the last digar
        if (*n_digar == *m_digar) {
            *m_digar = *m_digar ? (*m_digar << 1) : 4;
            digar = (digar1_t*)realloc(digar, *m_digar * sizeof(digar1_t));
        }
        digar[*n_digar].pos = d.pos; digar[*n_digar].type = d.type;
        digar[*n_digar].len = d.len; digar[*n_digar].alt_seq = d.alt_seq;
        digar[*n_digar].qi = d.qi; digar[*n_digar].is_low_qual = d.is_low_qual;
        (*n_digar)++;
    } else {
        if (digar[*n_digar-1].type == BAM_CINS) {
            uint8_t *alt_seq = (uint8_t*)realloc(digar[*n_digar-1].alt_seq, (digar[*n_digar-1].len + d.len) * sizeof(uint8_t));
            for (int i = digar[*n_digar-1].len; i < digar[*n_digar-1].len + d.len; ++i) 
                alt_seq[i] = d.alt_seq[i - digar[*n_digar-1].len];
            digar[*n_digar-1].alt_seq = alt_seq;
            free(d.alt_seq); // free the new alt_seq
        }
        digar[*n_digar-1].len += d.len; // merge INS/DEL with the last digar
    }
    return digar; // return the new digar pointer
}

void push_digar1(digar_t *digar, digar1_t d) {
    _uni_realloc(digar->digars, digar->n_digar, digar->m_digar, digar1_t);
    digar->digars[digar->n_digar].pos = d.pos; digar->digars[digar->n_digar].type = d.type;
    digar->digars[digar->n_digar].len = d.len; digar->digars[digar->n_digar].alt_seq = d.alt_seq;
    digar->digars[digar->n_digar].qi = d.qi; digar->digars[digar->n_digar].is_low_qual = d.is_low_qual;
    digar->n_digar++;
}

void print_digar1(digar1_t *digars, int n_digar, FILE *fp) {
    fprintf(fp, "pos\ttype\tlen\tqi\tbase\tis_low_qual\n");
    for (int i = 0; i < n_digar; ++i) {
        fprintf(fp, "%" PRId64 "\t%c\t%d\t%d\t%d\n", digars[i].pos, BAM_CIGAR_STR[digars[i].type], digars[i].len, digars[i].qi, digars[i].is_low_qual);
    }
}

int collect_noisy_region_len(cgranges_t *noisy_reg) {
    int len = 0;
    for (int i = 0; i < noisy_reg->n_r; ++i) {
        hts_pos_t noisy_reg_beg = cr_start(noisy_reg, i), noisy_reg_end = cr_end(noisy_reg, i);
        len += (noisy_reg_end - noisy_reg_beg + 1);
    }
    return len;
}

// XXX for ONT palindrome, ignore the long clipping
// if [sa_pos, sa_end] overlaps with [primary_pos, primary_end] by over 90% of the length, return 1, else return 0
int check_ont_palindrome(hts_pos_t primary_pos, hts_pos_t primary_end, hts_pos_t sa_pos, hts_pos_t sa_end) {
    hts_pos_t primary_len = primary_end - primary_pos + 1, sa_len = sa_end - sa_pos + 1;
    hts_pos_t overlap_len = 0;
    if (sa_pos <= primary_pos) {
        if (sa_end >= primary_end) overlap_len = primary_len;
        else if (sa_end >= primary_pos) overlap_len = sa_end - primary_pos + 1;
    } else if (sa_pos <= primary_end) {
        if (sa_end >= primary_end) overlap_len = primary_end - sa_pos + 1;
        else overlap_len = sa_len;
    }
    if (overlap_len >= primary_len * 0.9) return 1;
    return 0;
}

// return 
//   1: palindrome
//   0: not
int is_ont_palindrome_clip(const call_var_opt_t *opt, bam1_t *read) {
    int res = 0; char *sa_value = NULL; char *sa_saveptr = NULL;
    if (opt->is_ont)  { // check SA tag
        hts_pos_t primary_pos = read->core.pos+1, primary_end = bam_endpos(read);
        uint8_t *sa_tag = bam_aux_get(read, "SA");  // Get SA tag
        if (sa_tag) { // collect mapping chrom, pos, strand, end pos, MAPQ, NM
            // fprintf(stderr, "is_ont_palindrome_clip: %s\n", bam_get_qname(read));
            char *sa = bam_aux2Z(sa_tag);
            sa_value = strdup(sa);
            char *sa_entry = strtok_r(sa_value, ";", &sa_saveptr);
            while (sa_entry) {
                int sa_entry_len = strlen(sa_entry);
                char *rname = (char*)malloc((sa_entry_len+1) * sizeof(char));
                char *cigar = (char*)malloc((sa_entry_len+1) * sizeof(char));
                char strand; int pos, mapq, flag;
                sscanf(sa_entry, "%[^,],%d,%c,%[^,],%d,%d", rname, &pos, &strand, cigar, &mapq, &flag);
                // Estimate end position from CIGAR
                int32_t sa_end_pos = pos;
                char *cptr = cigar;
                while (*cptr) {
                    int len = strtol(cptr, &cptr, 10);
                    if (*cptr == 'M' || *cptr == 'D' || *cptr == '=' || *cptr == 'X') {
                        sa_end_pos += len;
                    }
                    cptr++;
                }
                free(cigar); free(rname);
                sa_entry = strtok_r(NULL, ";", &sa_saveptr);
                // if any SA entry is a palindrome, then it is not a noisy region
                if (check_ont_palindrome(primary_pos, primary_end, pos, sa_end_pos)) {
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "palindrome: %s %" PRIi64 "-%" PRIi64 ", sa: %d-%d\n", bam_get_qname(read), primary_pos, primary_end, pos, sa_end_pos);
                    res = 1;
                    break;
                }
            }
        }
    }
    if (sa_value != NULL) free(sa_value);
    return res;
}

// is_low_qual: low base quality
int collect_digar_from_eqx_cigar(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar) {
    bam1_t *read = chunk->reads[read_i];
    hts_pos_t pos = read->core.pos+1, qi = 0;
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;
    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    digar->beg = pos; digar->end = bam_endpos(read); digar->is_rev = bam_is_rev(read);
    uint32_t qlen = read->core.l_qseq;
    longcalld_copy_digar_read_buffers(chunk, read, digar);
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = chunk->whole_ref_len;
    xid_queue_t *q = init_xid_queue(rlen, max_s, win); // for noisy region
    hts_pos_t noisy_start = -1, noisy_end = -1; int cr_q_start = -1, cr_q_end = -1;
    int n_total_cand_vars = 0;

    int is_palindrome = is_ont_palindrome_clip(opt, read); int left_clip_is_palindrome = 0, right_clip_is_palindrome = 0;
    if (is_palindrome) {
        chunk->is_ont_palindrome[read_i] = 1;
        if (digar->is_rev) left_clip_is_palindrome = 1;
        else right_clip_is_palindrome = 1;
    }
    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]), len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CDIFF) {
            for (int j = 0; j < len; ++j) {
                _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                uint8_t *x_seq = (uint8_t*)malloc(sizeof(uint8_t));
                x_seq[0] = seq_nt16_int[bam_seqi(digar->bseq, qi)];
                if (digar->qual[qi] >= opt->min_bq) { // skip low-quality bases for noisy regions
                    push_xid_size_queue_win(q, pos, 1, 1, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                    set_digar(_digars+_n_digar, pos, op, 1, qi, 0, x_seq);
                } else set_digar(_digars+_n_digar, pos, op, 1, qi, 1, x_seq);
                _n_digar++; // push_xid_queue(q, pos, 1, 1);
                n_total_cand_vars++;
                pos++; qi++;
            }
        } else if (op == BAM_CEQUAL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, op, len, qi, 0, NULL);
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            // push_xid_queue_win(q, pos, len, 0, digar->noisy_regs, &noisy_start, &noisy_end);
            pos += len; qi += len;
        } else if (op == BAM_CDEL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((qi == 0 || digar->qual[qi-1] >= opt->min_bq) && digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, len, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, op, len, qi, 0, NULL);
            } else set_digar(_digars+_n_digar, pos, op, len, qi, 1, NULL);
            _n_digar++; // push_xid_queue(q, pos, len, 1);
            n_total_cand_vars++;
            pos += len;
        } else if (op == BAM_CINS) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            uint8_t *ins_seq = (uint8_t*)malloc(len * sizeof(uint8_t));
            for (int j = 0; j < len; ++j) ins_seq[j] = seq_nt16_int[bam_seqi(digar->bseq, qi+j)];
            int is_low_qual = 1;
            for (int _i = 0; _i < len; ++_i) {
                if (digar->qual[qi+_i] >= opt->min_bq) {
                    is_low_qual = 0; break;
                }
            }
            if (!is_low_qual) push_xid_size_queue_win(q, pos, 0, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
            set_digar(_digars+_n_digar, pos, op, len, qi, is_low_qual, ins_seq); // insertion
            _n_digar++; // push_xid_queue(q, pos, 0, 1);
            n_total_cand_vars++;
            qi += len;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((i == 0 && left_clip_is_palindrome) || (i != 0 && right_clip_is_palindrome))
                set_digar(_digars+_n_digar, pos, BAM_CHARD_CLIP, len, qi, 0, NULL); // palindrome
            else set_digar(_digars+_n_digar, pos, op, len, qi, 0, NULL); // clipping
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            // skip checking noisy region for clipping close to start/end of chromosome/contig: within 10bp
            if ((i == 0 && pos > 10) || (i != 0 && pos < tlen - 10)) {
                if (len > end_clip_reg) {
                    if (i == 0 && !left_clip_is_palindrome) {
                        if (pos > 1) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
                        n_total_cand_vars++;
                    } else if (i != 0 && !right_clip_is_palindrome) {
                        if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                        n_total_cand_vars++;
                    }
                }
            }
            if (op == BAM_CSOFT_CLIP) qi += len;
        } else if (op == BAM_CREF_SKIP) { // XXX no action for N op
            pos += len;
        } else if (op == BAM_CMATCH) {
            _err_error_exit("CIGAR operation 'M' is not expected in EQX CIGAR: %s\n", bam_get_qname(read));
        }
    }
    for (int i = 0; i < _n_digar; ++i) { // XXX not set noisy-region for digars
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
        var_size = var_size > (noisy_end - noisy_start + 1) ? var_size : (noisy_end - noisy_start + 1);
        // fprintf(stderr, "var_size: %d, noisy_start: %" PRIi64 ", noisy_end: %" PRIi64 "\n", var_size, noisy_start, noisy_end);
        cr_add(digar->noisy_regs, "cr", noisy_start-1, noisy_end, var_size);
    }
    cr_index(digar->noisy_regs);
    int skip = 0;
    // if total noisy region length > mapped length * 0.8, skip the read
    int total_noisy_reg_len = collect_noisy_region_len(digar->noisy_regs);
    int mapped_len = digar->end - digar->beg + 1;
    if (total_noisy_reg_len > mapped_len * max_noisy_frac_per_read || n_total_cand_vars > mapped_len * max_var_ratio_per_read) {
        // fprintf(stderr, "SkipRead: %s %d-noisy %d-var %d-ref_span\n", bam_get_qname(read), total_noisy_reg_len, n_total_cand_vars, mapped_len);
        skip = 1;
    } else {
        for (int i = 0; i < digar->noisy_regs->n_r; ++i) { // we may double-count reads with multiple nearby noisy regions 
                                                        // downside is minor, as we just consider an easy region as noisy, variant calling result should not be affected
                                                        // read1:  ---- [noisy] ---- [noisy] ----
                                                        // read2:  --------- [  noisy  ] --------
            if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "NoisyRead %s %s:%d-%d %d\n", bam_get_qname(read), chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
            if (is_overlap_reg(cr_start(digar->noisy_regs, i)+1, cr_end(digar->noisy_regs, i), reg_beg, reg_end)) {
                // set cr_label as the largest insertion in noisy region
                // the insertion size is used to determine the merge_distance:
                // ins_size: merge_dis
                //   <= 100: 100
                //    < 500: ins_size
                //   >= 500: 500
                cr_add(chunk_noisy_regs, "cr", cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
            }
        }
        if (LONGCALLD_VERBOSE >= 3) {
            fprintf(stderr, "DIGAR1: %s\n", bam_get_qname(read));
            print_digar(digar, stderr);
            if (digar->noisy_regs->n_r > 0) fprintf(stderr, "Noisy-%s\n", bam_get_qname(read));
            for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
                fprintf(stderr, "Noisy: %s:%d-%d\n", chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i));
            }
        }
    }
    free(_digars); free_xid_queue(q);
    return skip == 1 ? -1 : 0;
}

// XXX TODO make it consistent with collect_digar_from_eqx_cigar
int collect_digar_from_cs_tag(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar) {
    bam1_t *read = chunk->reads[read_i];
    uint8_t *cs_tag = bam_aux_get(read, "cs");
    if (cs_tag == NULL) { _err_error_exit("cs tag not found in the BAM file: %s", bam_get_qname(read)); }
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    if (n_cigar <= 0) { _err_error_exit("CIGAR not found in the BAM file: %s", bam_get_qname(read)); }
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    hts_pos_t pos = read->core.pos+1, qi = 0;
    digar->beg = pos; digar->end = bam_endpos(read); digar->is_rev = bam_is_rev(read);
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;

    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    uint32_t qlen = read->core.l_qseq;
    longcalld_copy_digar_read_buffers(chunk, read, digar);
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    char *cs = bam_aux2Z(cs_tag);
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = chunk->whole_ref_len;
    xid_queue_t *q = init_xid_queue(rlen, max_s, win);
    hts_pos_t noisy_start = -1, noisy_end = -1; int cr_q_start = -1, cr_q_end = -1;
    int n_total_cand_vars = 0;

    int is_palindrome = is_ont_palindrome_clip(opt, read); int left_clip_is_palindrome = 0, right_clip_is_palindrome = 0;
    if (is_palindrome) {
        chunk->is_ont_palindrome[read_i] = 1;
        if (digar->is_rev) left_clip_is_palindrome = 1;
        else right_clip_is_palindrome = 1;
    }

    // left-end clipping
    if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP || bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP) {
        int len = bam_cigar_oplen(cigar[0]);
        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
        if (left_clip_is_palindrome) set_digar(_digars+_n_digar, pos, BAM_CHARD_CLIP, len, qi, 0, NULL); // palindrome
        else set_digar(_digars+_n_digar, pos, bam_cigar_op(cigar[0]), len, qi, 0, NULL); // clipping
        _n_digar++; // push_xid_queue(q, pos, 0, 0);
        if (len > end_clip_reg && !left_clip_is_palindrome) {
            if (pos > 10) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
            // if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
            n_total_cand_vars++;
        }
        if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) qi += len;
    }
    while (*cs) {
        if (*cs == ':') { // identical seq length
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            int len = strtol(&cs[1], &cs, 10);
            set_digar(_digars+_n_digar, pos, BAM_CEQUAL, len, qi, 0, NULL);
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            pos += len; qi += len;
        } else if (*cs == '=') { // identical sequence
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            cs++; int len = 0;
            while (isalpha(*cs)) {
                len++; cs++;
            }
            set_digar(_digars+_n_digar, pos, BAM_CEQUAL, len, qi, 0, NULL);
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            pos += len; qi += len;
        } else if (*cs == '*') { // mismatch: 1bp
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            uint8_t *x_seq = (uint8_t*)malloc(sizeof(uint8_t));
            x_seq[0] = nst_nt4_table[(int)*(cs+2)];
            if (digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, 1, 1, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 0, x_seq);
            } else set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 1, x_seq);
            _n_digar++; pos++; qi++; cs += 3;
            n_total_cand_vars++;
        } else if (*cs == '+') { // insertion
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            cs++; int len = 0;
            while (isalpha(*cs)) {
                len++; cs++;
            }
            uint8_t *ins_seq = (uint8_t*)malloc(len * sizeof(uint8_t));
            for (int i = 0; i < len; ++i) ins_seq[i] = nst_nt4_table[(int)*(cs-len+i)];
            int is_low_qual = 1;
            for (int _i = 0; _i < len; ++_i) {
                if (digar->qual[qi+_i] >= opt->min_bq) {
                    is_low_qual = 0; break;
                }
            }
            if (!is_low_qual) push_xid_size_queue_win(q, pos, 0, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
            set_digar(_digars+_n_digar, pos, BAM_CINS, len, qi, is_low_qual, ins_seq);
            _n_digar++; qi += len;
            n_total_cand_vars++;
        } else if (*cs == '-') { // deletion
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            cs++; int len = 0;
            while (isalpha(*cs)) {
                len++; cs++;
            }
            if ((qi == 0 || digar->qual[qi-1] >= opt->min_bq) && digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, len, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 0, NULL);
            } else set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 1, NULL);
            _n_digar++; pos += len;
            n_total_cand_vars++;
        } else if (*cs == '~') { // intron/splice
            cs++;
            while (isalpha(*cs) || isdigit(*cs)) cs++;
        } else {
            _err_error_exit("Unexpected character in cs tag: %c", *cs);
        }
    }
    
    // right-end clipping
    if (bam_cigar_op(cigar[n_cigar-1]) == BAM_CSOFT_CLIP || bam_cigar_op(cigar[n_cigar-1]) == BAM_CHARD_CLIP) {
        int len = bam_cigar_oplen(cigar[n_cigar-1]);
        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
        if (right_clip_is_palindrome) set_digar(_digars+_n_digar, pos, BAM_CHARD_CLIP, len, qi, 0, NULL); // palindrome
        else set_digar(_digars+_n_digar, pos, bam_cigar_op(cigar[n_cigar-1]), len, qi, 0, NULL); // clipping
        _n_digar++; // push_xid_queue(q, pos, 0, 0);
        if (len > end_clip_reg && !right_clip_is_palindrome) {
            if (pos < tlen-10) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
            n_total_cand_vars++;
        }
        if (bam_cigar_op(cigar[n_cigar-1]) == BAM_CSOFT_CLIP) qi += len;
    }

    for (int i = 0; i < _n_digar; ++i) {
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
        var_size = var_size > (noisy_end - noisy_start + 1) ? var_size : (noisy_end - noisy_start + 1);
        // fprintf(stderr, "var_size: %d, noisy_start: %" PRIi64 ", noisy_end: %" PRIi64 "\n", var_size, noisy_start, noisy_end);
        cr_add(digar->noisy_regs, "cr", noisy_start-1, noisy_end, var_size);
    }
    cr_index(digar->noisy_regs);
    int skip = 0;
    int total_noisy_reg_len = collect_noisy_region_len(digar->noisy_regs);
    int mapped_len = digar->end - digar->beg + 1;
    if (total_noisy_reg_len > mapped_len * max_noisy_frac_per_read || n_total_cand_vars > mapped_len * max_var_ratio_per_read) {
        // fprintf(stderr, "Skip: %s %d-noisy %d-var %d-ref_span\n", bam_get_qname(read), total_noisy_reg_len, n_total_cand_vars, mapped_len);
        skip = 1;
    } else {
        for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
            if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "NoisyRead %s %s:%d-%d %d\n", bam_get_qname(read), chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
            if (is_overlap_reg(cr_start(digar->noisy_regs, i)+1, cr_end(digar->noisy_regs, i), reg_beg, reg_end))
                cr_add(chunk_noisy_regs, "cr", cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
        }
        if (LONGCALLD_VERBOSE >= 3) {
            fprintf(stderr, "DIGAR2: %s\n", bam_get_qname(read));
            print_digar(digar, stderr);
            if (digar->noisy_regs->n_r > 0) fprintf(stderr, "Noisy-%s\n", bam_get_qname(read));
            for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
                fprintf(stderr, "Noisy: %s:%d-%d\n", chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i));
            }
        }
    }
    free(_digars); free_xid_queue(q);
    return skip == 1 ? -1 : 0;
}

int collect_digar_from_MD_tag(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar) {
    bam1_t *read = chunk->reads[read_i];
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    uint8_t *s = bam_aux_get(read, "MD");
    if (s == NULL) { _err_error_exit("MD tag not found in the BAM file: %s", bam_get_qname(read)); }
    hts_pos_t pos = read->core.pos+1, qi = 0;
    digar->beg = pos; digar->end = bam_endpos(read); digar->is_rev = bam_is_rev(read);
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;

    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    uint32_t qlen = read->core.l_qseq;
    longcalld_copy_digar_read_buffers(chunk, read, digar);
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    char *md = bam_aux2Z(s); int md_i = 0;
    // printf("MD: %s\n", md);
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = chunk->whole_ref_len;
    xid_queue_t *q = init_xid_queue(rlen, max_s, win);
    hts_pos_t noisy_start = -1, noisy_end = -1; int cr_q_start = -1, cr_q_end = -1;
    int n_total_cand_vars = 0;

    int is_palindrome = is_ont_palindrome_clip(opt, read); int left_clip_is_palindrome = 0, right_clip_is_palindrome = 0;
    if (is_palindrome) {
        chunk->is_ont_palindrome[read_i] = 1;
        if (digar->is_rev) left_clip_is_palindrome = 1;
        else right_clip_is_palindrome = 1;
    }

    int last_eq_len = 0;
    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]), len = bam_cigar_oplen(cigar[i]);
        // fprintf(stderr, "%c %d, md_i: %d, MD: %s\n", BAM_CIGAR_STR[op], len, md_i, md);
        if (op == BAM_CMATCH) {
            // check if base is equal or mismatch based on MD tag
            int m_len = len, eq_len;
            while (1) { // MD
                if (last_eq_len > 0) {
                    if (last_eq_len >= m_len) {
                        eq_len = m_len;
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos, BAM_CEQUAL, eq_len, qi, 0, NULL);
                        _n_digar++; // push_xid_queue(q, pos, 0, 0);
                        pos += eq_len; qi += eq_len;
                        last_eq_len -= m_len; m_len = 0;
                    } else { // last_eq_len < m_len
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos, BAM_CEQUAL, last_eq_len, qi, 0, NULL);
                        _n_digar++; // push_xid_queue(q, pos, 0, 0);
                        pos += last_eq_len; qi += last_eq_len;
                        m_len -= last_eq_len; md_i = 0; // start of the next run
                        last_eq_len = 0;
                    }
                } else if (isdigit(md[md_i])) { // =
                    eq_len = strtol(&md[md_i], &md, 10);
                    if (eq_len > m_len) { //_err_error_exit("MD and CIGAR do not match: %s (%d : %d)", bam_get_qname(read), eq_len, m_len);
                        last_eq_len = eq_len - m_len;
                        eq_len = m_len;
                    } else if (eq_len == 0) {
                        md_i=0; continue;
                    }
                    _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                    set_digar(_digars+_n_digar, pos, BAM_CEQUAL, eq_len, qi, 0, NULL); // qual unset XXX
                    _n_digar++; // push_xid_queue(q, pos, 0, 0);
                    pos += eq_len; qi += eq_len;
                    m_len -= eq_len; md_i = 0; // start of the next run
                } else if (isalpha(md[md_i])) { // X
                    // if (qual[qi] >= min_bq) {
                    _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                    uint8_t *x_seq = (uint8_t*)malloc(sizeof(uint8_t));
                    x_seq[0] = seq_nt16_int[bam_seqi(digar->bseq, qi)];
                    if (digar->qual[qi] >= opt->min_bq) {
                        push_xid_size_queue_win(q, pos, 1, 1, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                        set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 0, x_seq);
                    } else set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 1, x_seq);
                    _n_digar++; // push_xid_queue(q, pos, 1, 1);
                    n_total_cand_vars++;
                    pos++; qi++; m_len -= 1; 
                    if (md[md_i+1] == '\0' || md[md_i+1] != '0') md_i++;
                    else md_i += 2; // skip the 0 after X
                } else _err_error_exit("MD and CIGAR do not match: %s (%s)", bam_get_qname(read), bam_aux2Z(s));
                if (m_len <= 0) break;
            }
        } else if (op == BAM_CDEL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((qi == 0 || digar->qual[qi-1] >= opt->min_bq) && digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, len, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 0, NULL);
            } else set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 1, NULL);
            _n_digar++; // push_xid_queue(q, pos, len, 1);
            n_total_cand_vars++;
            pos += len;
            // MD
            md_i++;
            while (md[md_i] && isalpha(md[md_i])) md_i++;
            if (md[md_i] == '0') md_i++; // skip 0 after D
        } else if (op == BAM_CINS) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            uint8_t *ins_seq = (uint8_t*)malloc(len * sizeof(uint8_t));
            for (int j = 0; j < len; ++j) ins_seq[j] = seq_nt16_int[bam_seqi(digar->bseq, qi+j)];
            int is_low_qual = 1;
            for (int _i = 0; _i < len; ++_i) {
                if (digar->qual[qi+_i] >= opt->min_bq) {
                    is_low_qual = 0; break;
                }
            }
            if (!is_low_qual) push_xid_size_queue_win(q, pos, 0, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
            set_digar(_digars+_n_digar, pos, BAM_CINS, len, qi, is_low_qual, ins_seq); // insertion
            _n_digar++; // push_xid_queue(q, pos, 0, 1);
            n_total_cand_vars++;
            qi += len;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((i == 0 && left_clip_is_palindrome) || (i != 0 && right_clip_is_palindrome))
                set_digar(_digars+_n_digar, pos, BAM_CHARD_CLIP, len, qi, 0, NULL); // palindromic: using HARD_CLIP
            else set_digar(_digars+_n_digar, pos, op, len, qi, 0, NULL); // normal clipping
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            if ((i == 0 && pos > 10) || (i != 0 && pos < tlen - 10)) {
                if (len > end_clip_reg) {
                    if (i == 0 && !left_clip_is_palindrome) {
                        if (pos > 1) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
                        n_total_cand_vars++;
                    } else if (i != 0 && !right_clip_is_palindrome) {
                        if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                        n_total_cand_vars++;
                    }
                }
            }
            if (op == BAM_CSOFT_CLIP) qi += len;
        } else if (op == BAM_CREF_SKIP) {
            pos += len;
        } else if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            _err_error_exit("CIGAR operation '=/X' is not expected: %s\n", bam_get_qname(read));
        }
    }
    for (int i = 0; i < _n_digar; ++i) {
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
        var_size = var_size > (noisy_end - noisy_start + 1) ? var_size : (noisy_end - noisy_start + 1);
        // fprintf(stderr, "var_size: %d, noisy_start: %" PRIi64 ", noisy_end: %" PRIi64 "\n", var_size, noisy_start, noisy_end);
        cr_add(digar->noisy_regs, "cr", noisy_start-1, noisy_end, var_size);
    }
    cr_index(digar->noisy_regs);
    int skip = 0;
    int total_noisy_reg_len = collect_noisy_region_len(digar->noisy_regs);
    int mapped_len = digar->end - digar->beg + 1;
    if (total_noisy_reg_len > mapped_len * max_noisy_frac_per_read || n_total_cand_vars > mapped_len * max_var_ratio_per_read) {
        // fprintf(stderr, "Skip: %s %d-noisy %d-var %d-ref_span\n", bam_get_qname(read), total_noisy_reg_len, n_total_cand_vars, mapped_len);
        skip = 1;
    } else {
        for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
            if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "NoisyRead %s %s:%d-%d %d\n", bam_get_qname(read), chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
            if (is_overlap_reg(cr_start(digar->noisy_regs, i)+1, cr_end(digar->noisy_regs, i), reg_beg, reg_end))
                cr_add(chunk_noisy_regs, "cr", cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
        }
        if (LONGCALLD_VERBOSE >= 3) {
            fprintf(stderr, "DIGAR2: %s\n", bam_get_qname(read));
            print_digar(digar, stderr);
            if (digar->noisy_regs->n_r > 0) fprintf(stderr, "Noisy-%s\n", bam_get_qname(read));
            for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
                fprintf(stderr, "Noisy: %s:%d-%d\n", chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i));
            }
        }
    }
    free(_digars); free_xid_queue(q);
    return skip == 1 ? -1 : 0;
}

int collect_digar_from_ref_seq(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar) {
    bam1_t *read = chunk->reads[read_i];
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    hts_pos_t pos = read->core.pos+1, qi = 0;
    digar->beg = pos; digar->end = bam_endpos(read); digar->is_rev = bam_is_rev(read);
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;

    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    uint32_t qlen = read->core.l_qseq;
    longcalld_copy_digar_read_buffers(chunk, read, digar);
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = chunk->whole_ref_len;
    xid_queue_t *q = init_xid_queue(rlen, max_s, win);
    hts_pos_t noisy_start = -1, noisy_end = -1; int cr_q_start = -1, cr_q_end = -1;
    int n_total_cand_vars = 0;

    int is_palindrome = is_ont_palindrome_clip(opt, read); int left_clip_is_palindrome = 0, right_clip_is_palindrome = 0;
    if (is_palindrome) {
        chunk->is_ont_palindrome[read_i] = 1;
        if (digar->is_rev) left_clip_is_palindrome = 1;
        else right_clip_is_palindrome = 1;
    }

    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) {
            int eq_len = 0;
            for (int j = 0; j < len; ++j) {
                // Get the reference base
                if (pos < ref_beg || pos > ref_end) {
                    if (LONGCALLD_VERBOSE >= 3) {
                        fprintf(stderr, "pos: %" PRIi64 " (%" PRIi64 "-%" PRIi64 ")\t", pos, ref_beg, ref_end);
                        fprintf(stderr, "Read exceed reference region sequence: %s", bam_get_qname(read));
                    }
                    pos++; qi++; continue;
                }
                int ref_base = nst_nt4_table[(int)ref_seq[pos-ref_beg]];
                int read_base = seq_nt16_int[bam_seqi(digar->bseq, qi)];
                if (ref_base != read_base) {
                    if (eq_len > 0) {
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos-eq_len, BAM_CEQUAL, eq_len, qi-eq_len, 0, NULL);
                        _n_digar++; // push_xid_queue(q, pos-eq_len, 0, 0);
                        eq_len = 0;
                    }
                    _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                    uint8_t *x_seq = (uint8_t*)malloc(sizeof(uint8_t));
                    x_seq[0] = read_base; //seq_nt16_int[bam_seqi(digar->bseq, qi)];
                    if (digar->qual[qi] >= opt->min_bq) {
                        push_xid_size_queue_win(q, pos, 1, 1, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                        set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 0, x_seq);
                    } else set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 1, x_seq);
                    _n_digar++;
                    n_total_cand_vars++;
                } else eq_len++;
                pos++; qi++;
            }
            if (eq_len > 0) {
                _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                set_digar(_digars+_n_digar, pos-eq_len, BAM_CEQUAL, eq_len, qi-eq_len, 0, NULL);
                _n_digar++; // push_xid_queue(q, pos-eq_len, 0, 0);
                eq_len = 0;
            }
        } else if (op == BAM_CDEL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((qi == 0 || digar->qual[qi-1] >= opt->min_bq) && digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, len, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 0, NULL);
            } else set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 1, NULL);
            _n_digar++; // push_xid_queue(q, pos, len, 1);
            n_total_cand_vars++;
            pos += len;
        } else if (op == BAM_CINS) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            uint8_t *ins_seq = (uint8_t*)malloc(len * sizeof(uint8_t));
            for (int j = 0; j < len; ++j) ins_seq[j] = seq_nt16_int[bam_seqi(digar->bseq, qi+j)];
            // set_digar(_digars+_n_digar, pos, BAM_CINS, len, qi); // insertion
            int is_low_qual = 1;
            for (int _i = 0; _i < len; ++_i) {
                if (digar->qual[qi+_i] >= opt->min_bq) {
                    is_low_qual = 0; break;
                }
            }
            if (!is_low_qual) push_xid_size_queue_win(q, pos, 0, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
            set_digar(_digars+_n_digar, pos, BAM_CINS, len, qi, is_low_qual, ins_seq);
            _n_digar++; //push_xid_queue(q, pos, 0, 1);
            n_total_cand_vars++;
            qi += len;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((i == 0 && left_clip_is_palindrome) || (i != 0 && right_clip_is_palindrome))
                set_digar(_digars+_n_digar, pos, BAM_CHARD_CLIP, len, qi, 0, NULL); // palindromic: using HARD_CLIP
            else set_digar(_digars+_n_digar, pos, op, len, qi, 0, NULL); // clipping 
            _n_digar++; //push_xid_queue(q, pos, 0, 0);
            if ((i==0 && pos >10) || (i!=0 && pos < tlen -10))  {
                if (len > end_clip_reg) {
                    if (i == 0 && !left_clip_is_palindrome) {
                        if (pos > 1) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
                        n_total_cand_vars++;
                    } else if (i != 0 && !right_clip_is_palindrome) {
                        if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                        n_total_cand_vars++;
                    }
                }
            }
            if (op == BAM_CSOFT_CLIP) qi += len;
        } else if (op == BAM_CREF_SKIP) {
            pos += len;
        } // else if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            // _err_error_exit("CIGAR operation '=/X' is not expected: %s\n", bam_get_qname(read));
        // }
    }
    for (int i = 0; i < _n_digar; ++i) {
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
        var_size = var_size > (noisy_end - noisy_start + 1) ? var_size : (noisy_end - noisy_start + 1);
        // fprintf(stderr, "var_size: %d, noisy_start: %" PRIi64 ", noisy_end: %" PRIi64 "\n", var_size, noisy_start, noisy_end);
        cr_add(digar->noisy_regs, "cr", noisy_start-1, noisy_end, var_size);
    }
    cr_index(digar->noisy_regs);
    int skip = 0;
    int total_noisy_reg_len = collect_noisy_region_len(digar->noisy_regs);
    int mapped_len = digar->end - digar->beg + 1;
    if (total_noisy_reg_len > mapped_len * max_noisy_frac_per_read || n_total_cand_vars > mapped_len * max_var_ratio_per_read) {
        skip = 1;
    } else {
        for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
            if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "NoisyRead %s %s:%d-%d %d\n", bam_get_qname(read), chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
            if (is_overlap_reg(cr_start(digar->noisy_regs, i)+1, cr_end(digar->noisy_regs, i), reg_beg, reg_end))
                cr_add(chunk_noisy_regs, "cr", cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
        }
        if (LONGCALLD_VERBOSE >= 3) {
            fprintf(stderr, "DIGAR3: %s\n", bam_get_qname(read));
            print_digar(digar, stderr);
            if (digar->noisy_regs->n_r > 0) fprintf(stderr, "Noisy-%s\n", bam_get_qname(read));
            for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
                fprintf(stderr, "Noisy: %s:%d-%d\n", chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i));
            }
        }
    }
    free(_digars); free_xid_queue(q);
    return skip == 1 ? -1 : 0;
}

int bam_chunk_init0(bam_chunk_t *chunk, const struct call_var_opt_t *opt, int n_reads, int n_bam) {
    // input
    chunk->n_reads = 0; chunk->m_reads = n_reads; chunk->ordered_read_ids = (int*)malloc(n_reads * sizeof(int));
    if (need_read_names(opt)) {
        chunk->read_names = (char**)malloc(n_reads * sizeof(char*));
        for (int i = 0; i < n_reads; i++) chunk->read_names[i] = NULL;
    }
    chunk->n_bam = n_bam;
    chunk->bam_read_off = (int*)calloc(n_bam, sizeof(int));
    chunk->n_up_ovlp_reads = (int*)calloc(n_bam, sizeof(int)); chunk->n_up_ovlp_skip_reads = (int*)calloc(n_bam, sizeof(int));
    chunk->up_ovlp_read_i = (int**)malloc(n_bam * sizeof(int*));
    chunk->n_down_ovlp_reads = (int*)calloc(n_bam, sizeof(int)); chunk->n_down_ovlp_skip_reads = (int*)calloc(n_bam, sizeof(int));
    chunk->down_ovlp_read_i = (int**)malloc(n_bam * sizeof(int*));
    for (int i = 0; i < n_bam; ++i) {
        chunk->up_ovlp_read_i[i] = (int*)malloc(n_reads * sizeof(int));
        chunk->down_ovlp_read_i[i] = (int*)malloc(n_reads * sizeof(int));
    }
    chunk->reads = (bam1_t**)malloc(n_reads * sizeof(bam1_t*));
    chunk->ref_seq = NULL;
    chunk->low_comp_cr = NULL;
    // intermediate
    chunk->qual_counts = (int*)calloc(256, sizeof(int));
    chunk->is_skipped = (uint8_t*)calloc(n_reads, sizeof(uint8_t));
    chunk->read_sigs = (uint64_t*)calloc(n_reads, sizeof(uint64_t));
    chunk->read_sample_idx = (uint8_t*)calloc(n_reads, sizeof(uint8_t));
    chunk->n_samples = opt->n_samples;
    chunk->n_sample_reads = (int*)calloc(opt->n_samples > 0 ? opt->n_samples : 1, sizeof(int));
    chunk->is_invalid = 0;
    chunk->n_clean_agree_snps = (int*)malloc(n_reads * sizeof(int));
    chunk->n_clean_conflict_snps = (int*)malloc(n_reads * sizeof(int));
    chunk->is_ont_palindrome = (uint8_t*)calloc(n_reads, sizeof(uint8_t));
    chunk->is_skipped_for_somatic = (uint8_t*)calloc(n_reads, sizeof(uint8_t));
    chunk->digars = (digar_t*)calloc(n_reads, sizeof(digar_t));
    for (int i = 0; i < n_reads; i++) {
        chunk->reads[i] = bam_init1();
        chunk->digars[i].n_digar = chunk->digars[i].m_digar = 0;
    }
    // noisy regions
    chunk->chunk_noisy_regs = NULL; chunk->noisy_reg_to_reads = NULL; chunk->noisy_reg_to_n_reads = NULL;
    // variant
    chunk->n_cand_vars = 0; chunk->cand_vars = NULL; chunk->var_i_to_cate = NULL;
    chunk->read_var_profile = NULL; chunk->read_var_cr = NULL;
    chunk->som_art_mask_hap = NULL; chunk->som_art_mask_pooled = NULL; chunk->n_som_art = 0;
    chunk->var_noisy_read_cov_cr = NULL; chunk->var_noisy_read_err_cr = NULL;
    chunk->var_noisy_read_marks = NULL; chunk->var_noisy_read_mark_id = 0;
    // output:
    chunk->flip_hap = 0;
    chunk->haps = (int*)calloc(n_reads, sizeof(int)); // read-wise
    chunk->phase_scores = (int*)calloc(n_reads, sizeof(int)); // read-wise
    chunk->phase_sets = (hts_pos_t*)malloc(n_reads * sizeof(hts_pos_t)); // read-wise
    for (int i = 0; i < n_reads; i++) chunk->phase_sets[i] = -1; // -1 means no phase set
    return 0;
}

int bam_chunk_realloc(bam_chunk_t *chunk, const struct call_var_opt_t *opt) {
    int m_reads = chunk->m_reads * 2;
    if (chunk->var_noisy_read_cov_cr != NULL) {
        cr_destroy(chunk->var_noisy_read_cov_cr);
        chunk->var_noisy_read_cov_cr = NULL;
    }
    if (chunk->var_noisy_read_err_cr != NULL) {
        cr_destroy(chunk->var_noisy_read_err_cr);
        chunk->var_noisy_read_err_cr = NULL;
    }
    if (chunk->var_noisy_read_marks != NULL) {
        free(chunk->var_noisy_read_marks);
        chunk->var_noisy_read_marks = NULL;
    }
    chunk->var_noisy_read_mark_id = 0;
    chunk->reads = (bam1_t**)realloc(chunk->reads, m_reads * sizeof(bam1_t*));
    if (need_read_names(opt)) 
        chunk->read_names = (char**)realloc(chunk->read_names, m_reads * sizeof(char*));
    chunk->ordered_read_ids = (int*)realloc(chunk->ordered_read_ids, m_reads * sizeof(int));
    for (int i = 0; i < chunk->n_bam; ++i) {
        chunk->up_ovlp_read_i[i] = (int*)realloc(chunk->up_ovlp_read_i[i], m_reads * sizeof(int));
        chunk->down_ovlp_read_i[i] = (int*)realloc(chunk->down_ovlp_read_i[i], m_reads * sizeof(int));
    }
    chunk->is_skipped = (uint8_t*)realloc(chunk->is_skipped, m_reads * sizeof(uint8_t));
    chunk->read_sigs = (uint64_t*)realloc(chunk->read_sigs, m_reads * sizeof(uint64_t));
    chunk->read_sample_idx = (uint8_t*)realloc(chunk->read_sample_idx, m_reads * sizeof(uint8_t));
    chunk->n_clean_agree_snps = (int*)realloc(chunk->n_clean_agree_snps, m_reads * sizeof(int));
    chunk->n_clean_conflict_snps = (int*)realloc(chunk->n_clean_conflict_snps, m_reads * sizeof(int));
    chunk->is_ont_palindrome = (uint8_t*)realloc(chunk->is_ont_palindrome, m_reads * sizeof(uint8_t));
    chunk->is_skipped_for_somatic = (uint8_t*)realloc(chunk->is_skipped_for_somatic, m_reads * sizeof(uint8_t));
    chunk->digars = (digar_t*)realloc(chunk->digars, m_reads * sizeof(digar_t));
    chunk->haps = (int*)realloc(chunk->haps, m_reads * sizeof(int));
    chunk->phase_scores = (int*)realloc(chunk->phase_scores, m_reads * sizeof(int));
    chunk->phase_sets = (hts_pos_t*)realloc(chunk->phase_sets, m_reads * sizeof(hts_pos_t));
    for (int i = chunk->m_reads; i < m_reads; i++) {
        chunk->reads[i] = bam_init1();
        chunk->digars[i].n_digar = chunk->digars[i].m_digar = 0;
        chunk->is_skipped[i] = 0; chunk->read_sigs[i] = 0; chunk->read_sample_idx[i] = 0;
        chunk->phase_scores[i] = 0; chunk->haps[i] = 0; chunk->phase_sets[i] = -1;
        chunk->is_skipped_for_somatic[i] = 0;
    }
    chunk->m_reads = m_reads;
    return 0;
}

void free_digar1(digar1_t *digar1, int n_digar) {
    for (int i = 0; i < n_digar; i++) {
        if (digar1[i].type == BAM_CINS || digar1[i].type == BAM_CDIFF) {
            free(digar1[i].alt_seq);
        }
    }
    free(digar1);
}

static void bam_chunk_clear_var_noisy_read_cache(bam_chunk_t *chunk) {
    if (chunk->var_noisy_read_cov_cr != NULL) {
        cr_destroy(chunk->var_noisy_read_cov_cr);
        chunk->var_noisy_read_cov_cr = NULL;
    }
    if (chunk->var_noisy_read_err_cr != NULL) {
        cr_destroy(chunk->var_noisy_read_err_cr);
        chunk->var_noisy_read_err_cr = NULL;
    }
    free(chunk->var_noisy_read_marks);
    chunk->var_noisy_read_marks = NULL;
    chunk->var_noisy_read_mark_id = 0;
}

// UNUSED. Superseded by bam_chunk_post_free(); nothing calls this, and it has
// already drifted out of sync with the fields the chunk actually owns. Left in
// place rather than deleted as part of an unrelated change, but do not call it:
// it would both miss newer allocations and double-free older ones.
void bam_chunk_free(bam_chunk_t *chunk) {
    free(chunk->qual_counts);
    if (chunk->ref_seq != NULL) free(chunk->ref_seq);
    if (chunk->low_comp_cr != NULL) cr_destroy(chunk->low_comp_cr);
    bam_chunk_clear_var_noisy_read_cache(chunk);
    for (int i = 0; i < chunk->m_reads; i++) {
        if (chunk->digars[i].m_digar > 0) {
            free_digar1(chunk->digars[i].digars, chunk->digars[i].n_digar);
            free(chunk->digars[i].bseq); free(chunk->digars[i].qual);
            cr_destroy(chunk->digars[i].noisy_regs);
        }
    }
    free(chunk->ordered_read_ids);
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < chunk->m_reads; i++) {
            bam_destroy1(chunk->reads[i]);
        }
        free(chunk->reads);
    }
    if (chunk->noisy_reg_to_reads != NULL) {
        for (int i = 0; i < chunk->chunk_noisy_regs->n_r; i++) {
            free(chunk->noisy_reg_to_reads[i]);
        } free(chunk->noisy_reg_to_reads);
    }
    if (chunk->noisy_reg_to_n_reads != NULL) free(chunk->noisy_reg_to_n_reads);
    if (chunk->chunk_noisy_regs != NULL) cr_destroy(chunk->chunk_noisy_regs);
    if (chunk->cand_vars != NULL) free_cand_vars(chunk->cand_vars, chunk->n_cand_vars);
    if (chunk->var_i_to_cate != NULL) free(chunk->var_i_to_cate);

    if (chunk->read_var_profile != NULL) free_read_var_profile(chunk->read_var_profile, chunk->n_reads);
    if (chunk->read_var_cr != NULL) cr_destroy(chunk->read_var_cr);
    free(chunk->digars); free(chunk->is_skipped); 
    free(chunk->n_clean_agree_snps); free(chunk->n_clean_conflict_snps);
    free(chunk->is_ont_palindrome); free(chunk->is_skipped_for_somatic);
    free(chunk->haps); free(chunk->phase_sets); free(chunk->phase_scores);
    for (int i = 0; i < chunk->n_bam; ++i) {
        if (chunk->up_ovlp_read_i[i] != NULL) free(chunk->up_ovlp_read_i[i]);
        if (chunk->down_ovlp_read_i[i] != NULL) free(chunk->down_ovlp_read_i[i]);
    }
    free(chunk->n_up_ovlp_reads); free(chunk->n_down_ovlp_reads); free(chunk->n_up_ovlp_skip_reads); free(chunk->n_down_ovlp_skip_reads);
}

void bam_chunk_free_digar(bam_chunk_t *chunk) {
    for (int i = 0; i < chunk->m_reads; i++) {
        if (chunk->digars[i].m_digar > 0) {
            free_digar1(chunk->digars[i].digars, chunk->digars[i].n_digar);
            free(chunk->digars[i].bseq); free(chunk->digars[i].qual);
            cr_destroy(chunk->digars[i].noisy_regs);
        }
    }
    free(chunk->digars);
}

void bam_chunk_post_free(bam_chunk_t *chunk, const struct call_var_opt_t *opt) {
    free(chunk->qual_counts);
    if (chunk->ref_seq != NULL) free(chunk->ref_seq);
    if (chunk->cand_vars != NULL) free_cand_vars(chunk->cand_vars, chunk->n_cand_vars);
    if (chunk->var_i_to_cate != NULL) free(chunk->var_i_to_cate);
    bam_chunk_clear_var_noisy_read_cache(chunk);
    free(chunk->is_skipped);
    free(chunk->read_sigs); // needed by write_read_to_bam in step 2, so not freed in mid_free
    free(chunk->read_sample_idx); free(chunk->n_sample_reads);
    free(chunk->n_clean_agree_snps); free(chunk->n_clean_conflict_snps);
    free(chunk->is_ont_palindrome); free(chunk->is_skipped_for_somatic);
    free(chunk->haps); free(chunk->phase_sets); free(chunk->phase_scores);
    for (int i = 0; i < chunk->n_bam; ++i) {
        if (chunk->up_ovlp_read_i[i] != NULL) free(chunk->up_ovlp_read_i[i]);
        if (chunk->down_ovlp_read_i[i] != NULL) free(chunk->down_ovlp_read_i[i]);
    }
    free(chunk->bam_read_off);
    free(chunk->var_slots);
    free(chunk->som_art_mask_hap); free(chunk->som_art_mask_pooled);
    chunk->som_art_mask_hap = NULL; chunk->som_art_mask_pooled = NULL; chunk->n_som_art = 0;
    free(chunk->n_up_ovlp_reads); free(chunk->n_down_ovlp_reads); free(chunk->n_up_ovlp_skip_reads); free(chunk->n_down_ovlp_skip_reads);
    free(chunk->up_ovlp_read_i); free(chunk->down_ovlp_read_i);
    if (opt->out_aln_fp != NULL && opt->refine_bam) bam_chunk_free_digar(chunk);
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < chunk->m_reads; i++) {
            bam_destroy1(chunk->reads[i]);
        }
        free(chunk->reads);
    }
    free(chunk->ordered_read_ids);
    // save read_var_profile for output_sv_rnames or output_somatic_sv_rnames, which will be used in make_variants
    if (chunk->read_var_profile != NULL) free_read_var_profile(chunk->read_var_profile, chunk->n_reads);
    if (need_read_names(opt)) {
        for (int i = 0; i < chunk->n_reads; i++) free(chunk->read_names[i]); 
        free(chunk->read_names);
    }
}

// free variables that will not be used in stitch & make variants
// Free what step 2 does not need. Every pointer freed here is set to NULL
// afterwards, which is not tidiness: bam_chunk_post_free() frees the same
// low_comp_cr and read_var_cr behind a "!= NULL" check, so leaving them
// dangling made that check always true and freed them twice. Anything in step 2
// that reaches for one of these now sees NULL and can say so, instead of
// reading freed memory -- which is how a cgranges query from step 2 came to
// crash inside cr_get_ctg().
void bam_chunk_mid_free(bam_chunk_t *chunk, const struct call_var_opt_t *opt) {
    if (chunk->low_comp_cr != NULL) { cr_destroy(chunk->low_comp_cr); chunk->low_comp_cr = NULL; }
    if (opt->out_aln_fp == NULL || opt->refine_bam == 0) bam_chunk_free_digar(chunk);
    bam_chunk_clear_var_noisy_read_cache(chunk);
    // The per-region read lists are sized by chunk_noisy_regs->n_r, so they have
    // to be freed while that object still exists.
    if (chunk->noisy_reg_to_reads != NULL) {
        int n_noisy_regs = chunk->chunk_noisy_regs != NULL ? chunk->chunk_noisy_regs->n_r : 0;
        for (int i = 0; i < n_noisy_regs; i++) free(chunk->noisy_reg_to_reads[i]);
        free(chunk->noisy_reg_to_reads); chunk->noisy_reg_to_reads = NULL;
    }
    if (chunk->chunk_noisy_regs != NULL) { cr_destroy(chunk->chunk_noisy_regs); chunk->chunk_noisy_regs = NULL; }
    if (chunk->noisy_reg_to_n_reads != NULL) { free(chunk->noisy_reg_to_n_reads); chunk->noisy_reg_to_n_reads = NULL; }
    if (chunk->read_var_cr != NULL) { cr_destroy(chunk->read_var_cr); chunk->read_var_cr = NULL; }
}

void bam_chunks_mid_free(bam_chunk_t *chunks, int n_chunks, const struct call_var_opt_t *opt) {
    for (int i = 0; i < n_chunks; i++) {
        bam_chunk_mid_free(chunks+i, opt);
    }
}

void bam_chunks_post_free(bam_chunk_t *chunks, int n_chunks, const struct call_var_opt_t *opt) {
    for (int i = 0; i < n_chunks; i++) {
        bam_chunk_post_free(chunks+i, opt);
    }
    free(chunks);
}

void get_bam_chunk_reg_ref_seq0(faidx_t *fai, bam_chunk_t *chunk, hts_pos_t beg, hts_pos_t end) {
    assert(fai != NULL); // faidx_t *fai = fai_load(ref_fasta);
    int ref_seq_len = faidx_seq_len(fai, chunk->tname);
    int flank_len = 50000; // [reg_beg-flank_len, reg_end+flank_len]
    hts_pos_t ref_beg = MAX_OF_TWO(flank_len, beg-1)-flank_len;
    hts_pos_t ref_end = MIN_OF_TWO(ref_seq_len-flank_len-1, end-1)+flank_len;
    chunk->whole_ref_len = ref_seq_len;

    int len;
    chunk->ref_seq = faidx_fetch_seq(fai, chunk->tname, ref_beg, ref_end, &len); // ref_beg & ref_end: 0-based
    if (chunk->ref_seq == NULL || len <= 0) {
        _err_error("Failed to fetch reference sequence for %s:%" PRId64 "-%" PRId64 " from the fasta file.\n", chunk->tname, beg, end);
        _err_error_exit("Please make sure the reference genome and the alignment file match.\n");
                         
    }
    chunk->ref_beg = ref_beg+1;   // 1-based
    chunk->ref_end = ref_beg+len; // 1-based

    // collect low-complexity regions
    chunk->low_comp_cr = cr_init();
    uint64_t *r; int n=0, T=LONGCALLD_SDUST_T, W=LONGCALLD_SDUST_W;
    r = sdust(0, (uint8_t*)chunk->ref_seq-chunk->ref_beg+chunk->reg_beg, chunk->reg_end - chunk->reg_beg+1, T, W, &n);
    for (int i = 0; i < n; ++i) {
        // fprintf(stderr, "low_comp_cr: %s %d-%d\n", chunk->tname, chunk->reg_beg+(int)(r[i]>>32)-1, chunk->reg_beg+(int)r[i]-1);
        cr_add(chunk->low_comp_cr, "cr", chunk->reg_beg+(int)(r[i]>>32)-1, chunk->reg_beg+(int)r[i]-1, 0);
    }
    cr_index(chunk->low_comp_cr); free(r);
}

static int is_ovlp_with_prev_region(const struct call_var_pl_t *pl, bam_chunk_t *chunk, bam1_t *read) {
    int tid = chunk->tid, reg_chunk_i = chunk->reg_chunk_i, reg_i = chunk->reg_i;
    if (reg_i <= 0) return 0;
    int pre_reg_chunk_i = reg_chunk_i, pre_reg_i = reg_i-1;
    int pre_tid = pl->reg_chunks[pre_reg_chunk_i].reg_tids[pre_reg_i];
    if (tid != pre_tid) return 0;
    hts_pos_t pre_reg_beg = pl->reg_chunks[pre_reg_chunk_i].reg_begs[pre_reg_i];
    hts_pos_t pre_reg_end = pl->reg_chunks[pre_reg_chunk_i].reg_ends[pre_reg_i];
    hts_pos_t read_beg = read->core.pos+1, read_end = bam_endpos(read);
    if (read_end < pre_reg_beg || read_beg > pre_reg_end) return 0;
    else {
        // fprintf(stderr, "ovlp_pre: %s %" PRIi64 "-%" PRIi64 " %" PRId64 "-%" PRId64 "\n", bam_get_qname(read), read_beg, read_end, chunk->reg_beg, chunk->reg_end);
        return 1;
    }
}

static int is_ovlp_with_next_region(const struct call_var_pl_t *pl, bam_chunk_t *chunk, bam1_t *read) {
    int tid = chunk->tid, reg_chunk_i = chunk->reg_chunk_i, reg_i = chunk->reg_i;
    if (reg_i >= pl->reg_chunks[reg_chunk_i].n_regions-1) return 0;
    int next_reg_chunk_i = reg_chunk_i, next_reg_i = reg_i+1;
    int next_tid = pl->reg_chunks[next_reg_chunk_i].reg_tids[next_reg_i];
    if (tid != next_tid) return 0;
    hts_pos_t next_reg_beg = pl->reg_chunks[next_reg_chunk_i].reg_begs[next_reg_i];
    hts_pos_t next_reg_end = pl->reg_chunks[next_reg_chunk_i].reg_ends[next_reg_i];
    hts_pos_t read_beg = read->core.pos+1, read_end = bam_endpos(read);
    if (read_end < next_reg_beg || read_beg > next_reg_end) return 0;
    else {
        // fprintf(stderr, "ovlp_next: %s %" PRIi64 "-%" PRIi64 " %" PRId64 "-%" PRId64 "\n", bam_get_qname(read), read_beg, read_end, chunk->reg_beg, chunk->reg_end);
        return 1;
    }
}

typedef struct {
    int read_i;
    hts_pos_t pos, end;
    int NM;
    char *qname;
} bam_read_sort_t;

int comp_bam_read_sort(const void *a, const void *b) {
    const bam_read_sort_t *ra = (const bam_read_sort_t *)a;
    const bam_read_sort_t *rb = (const bam_read_sort_t *)b;
    if (ra->pos != rb->pos) return (ra->pos < rb->pos) ? -1 : 1;
    if (ra->end != rb->end) return (ra->end < rb->end) ? 1 : -1;
    if (ra->NM != rb->NM) return (ra->NM < rb->NM) ? -1 : 1;
    return strcmp(ra->qname, rb->qname);
}

int bam_get_NM(bam1_t *read) {
    uint8_t *nm_ptr = bam_aux_get(read, "NM");
    if (nm_ptr != NULL) {
        return bam_aux2i(nm_ptr);
    } else {
        return 0;
    }
}

void sort_chunk_reads(bam_chunk_t *chunk) {
    bam_read_sort_t *read_sorts = (bam_read_sort_t*)malloc(chunk->n_reads * sizeof(bam_read_sort_t));
    for (int i = 0; i < chunk->n_reads; ++i) {
        read_sorts[i].read_i = i;
        read_sorts[i].pos = chunk->reads[i]->core.pos;
        read_sorts[i].end = bam_endpos(chunk->reads[i]);
        read_sorts[i].NM = bam_get_NM(chunk->reads[i]);
        read_sorts[i].qname = bam_get_qname(chunk->reads[i]);
    }
    // sort by pos, end, NM, qname
    qsort(read_sorts, chunk->n_reads, sizeof(bam_read_sort_t), comp_bam_read_sort);
    for (int i = 0; i < chunk->n_reads; ++i) {
        chunk->ordered_read_ids[i] = read_sorts[i].read_i;
    }
    free(read_sorts);
}

// load ref_seq/read in reg_chunks[reg_chunk_i]->tid/beg/end[reg_i] to chunks
int collect_ref_seq_bam_main(const struct call_var_pl_t *pl, struct call_var_io_aux_t *io_aux, int reg_chunk_i, int reg_i, bam_chunk_t *chunk) {
    assert(reg_chunk_i < pl->n_reg_chunks); assert(reg_i < pl->reg_chunks[reg_chunk_i].n_regions);
    int tid = pl->reg_chunks[reg_chunk_i].reg_tids[reg_i];
    hts_pos_t reg_beg = pl->reg_chunks[reg_chunk_i].reg_begs[reg_i];
    hts_pos_t reg_end = pl->reg_chunks[reg_chunk_i].reg_ends[reg_i];
    call_var_opt_t *opt = pl->opt;

    int min_mapq = opt->min_mq;
    bam_chunk_init0(chunk, opt, 4096, io_aux->n_bam);
    chunk->reg_chunk_i = reg_chunk_i; chunk->reg_i = reg_i;
    chunk->tid = tid; chunk->tname = io_aux->headers[0]->target_name[tid]; chunk->reg_beg = reg_beg; chunk->reg_end = reg_end;
    hts_pos_t min_read_beg = reg_beg, max_read_end = reg_end;
    for (int i = 0; i < io_aux->n_bam; ++i) {
        // Where this file's reads begin in chunk->reads[]. write_read_to_bam()
        // emits in coordinate order across files, so it needs this to map a
        // record back to the read_i the loader gave it.
        chunk->bam_read_off[i] = chunk->n_reads;
        // create iterator for the region
        hts_itr_t *iter = sam_itr_queryi(io_aux->idxs[i], tid, reg_beg-1, reg_end); // (reg_beg-1, reg_end] ==> [reg_beg, reg_end]
        if (iter == NULL) {
            _err_warning("Failed to create iterator for region: %s:%" PRId64 "-%" PRId64 "\n", io_aux->headers[i]->target_name[tid], reg_beg, reg_end);
            // The chunk is only half-loaded; mark it so neither variant calling
            // nor either output path touches it. Previously this return value was
            // ignored and the partial chunk was processed anyway.
            chunk->is_invalid = 1;
            return -1;
        }
        // load bam records
        samFile *in_bam = io_aux->bams[i];
        // fprintf(stderr, "Loading reads from region %s:%" PRId64 "-%" PRId64 "\n", io_aux->headers[i]->target_name[tid], reg_beg, reg_end);
        while (sam_itr_next(in_bam, iter, chunk->reads[chunk->n_reads]) >= 0) {
            // fprintf(stderr, "load read %s\n", bam_get_qname(chunk->reads[chunk->n_reads]));
            if (!bam_record_is_callable(chunk->reads[chunk->n_reads], min_mapq)) {
                if (is_ovlp_with_prev_region(pl, chunk, chunk->reads[chunk->n_reads])) chunk->n_up_ovlp_skip_reads[i]++;
                if (is_ovlp_with_next_region(pl, chunk, chunk->reads[chunk->n_reads])) chunk->n_down_ovlp_skip_reads[i]++;
                continue; // BAM_FSUPPLEMENTARY
            }
            // fprintf(stderr, "CHUNK-READ: %s %d-%d\n", bam_get_qname(chunk->reads[chunk->n_reads]), chunk->reg_beg, chunk->reg_end);
            // check if read is overlapping with previous region
            if (is_ovlp_with_prev_region(pl, chunk, chunk->reads[chunk->n_reads])) {
                chunk->up_ovlp_read_i[i][chunk->n_up_ovlp_reads[i]++] = chunk->n_reads;
            }
            if (is_ovlp_with_next_region(pl, chunk, chunk->reads[chunk->n_reads])) {
                chunk->down_ovlp_read_i[i][chunk->n_down_ovlp_reads[i]++] = chunk->n_reads;
            }
            // update min_read_beg, max_read_end
            if (chunk->reads[chunk->n_reads]->core.pos+1 < min_read_beg) min_read_beg = chunk->reads[chunk->n_reads]->core.pos+1;
            if (bam_endpos(chunk->reads[chunk->n_reads]) > max_read_end) max_read_end = bam_endpos(chunk->reads[chunk->n_reads]);
            // check if read is overlapping with next region
            if (need_read_names(opt)) {
                chunk->read_names[chunk->n_reads] = strdup(bam_get_qname(chunk->reads[chunk->n_reads]));
            }
            chunk->read_sigs[chunk->n_reads] = bam_record_sig(chunk->reads[chunk->n_reads]);
            // Each input file belongs to exactly one sample, and the per-BAM loop
            // fills chunk->reads with one monotonic counter, so every read_i
            // carries its sample origin without any extra bookkeeping.
            {
                int smp = opt->bam_i_to_sample != NULL ? opt->bam_i_to_sample[i] : 0;
                chunk->read_sample_idx[chunk->n_reads] = (uint8_t)smp;
                if (smp >= 0 && smp < chunk->n_samples) chunk->n_sample_reads[smp]++;
            }
            chunk->n_reads++;
            if (chunk->n_reads == chunk->m_reads) bam_chunk_realloc(chunk, opt);
        }
        bam_itr_destroy(iter);
    }
    if (chunk->n_reads <= 0) return 0;
    // load ref seq
    get_bam_chunk_reg_ref_seq0(io_aux->fai, chunk, min_read_beg, max_read_end);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "CHUNK: tname: %s, tid: %d, beg: %" PRId64 ", end: %" PRId64 ", n_reads: %d\n", chunk->tname, chunk->tid, chunk->reg_beg, chunk->reg_end, chunk->n_reads);
    }
    sort_chunk_reads(chunk);
    return chunk->n_reads;
}

int cigar_is_idential(const uint32_t *cigar, int n_cigar, const uint32_t *new_cigar, int new_n_cigar) {
    if (n_cigar != new_n_cigar) return 0; // Different lengths
    for (int i = 0; i < n_cigar; ++i) {
        if (cigar[i] != new_cigar[i]) return 0; // Different CIGAR operations
    }
    return 1; // Identical CIGAR strings
}

static inline int get_nm_from_digar(digar1_t *d, int n_digar) {
    // Calculate NM tag from CIGAR string
    int NM = 0;
    for (int i = 0; i < n_digar; ++i) {
        int op = d[i].type;
        int len = d[i].len;
        if (op == BAM_CDIFF || op == BAM_CINS || op == BAM_CDEL) {
            NM += len;
        }
    }
    return NM;
}

static inline kstring_t *get_md_from_digar(digar1_t *d, int n_digar, const char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end) {
    // Calculate MD tag from CIGAR string
    kstring_t *MD_str = (kstring_t*)calloc(1, sizeof(kstring_t));
    hts_pos_t pos; // Start position in reference sequence
    int eq_len = 0, prev_was_del = 0;
    for (int i = 0; i < n_digar; ++i) {
        int op = d[i].type, len = d[i].len;
        if (op == BAM_CEQUAL) {
            eq_len += len;
        } else if (op == BAM_CDIFF) {
            kputw(eq_len, MD_str);
            for (int j = 0; j < len; ++j) {
                pos = d[i].pos + j;
                char ref_base = 'N';
                if (pos >= ref_beg && pos < ref_end) ref_base = ref_seq[pos - ref_beg];
                kputc(ref_base, MD_str);
            }
            eq_len = 0;
        } else if (op == BAM_CDEL) {
            // Consecutive D ops describe one contiguous deletion, so extend the
            // run rather than starting a second "0^..." group. samtools calmd
            // writes "^CC" where emitting each op separately gives "^C0^C".
            if (!prev_was_del) {
                kputw(eq_len, MD_str);
                kputc('^', MD_str);
            }
            for (int j = 0; j < len; ++j) {
                pos = d[i].pos + j;
                char ref_base = 'N';
                if (pos >= ref_beg && pos < ref_end) ref_base = ref_seq[pos - ref_beg];
                kputc(ref_base, MD_str);
            }
            eq_len = 0; prev_was_del = 1; continue;
        }
        prev_was_del = 0;
    }
    // Unconditional: the MD grammar is [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*, so the
    // string must END with a number as well as start with one. Emitting this
    // only when eq_len > 0 left a trailing "^A" whenever the alignment ended on
    // a deletion or mismatch, which is not valid MD.
    kputw(eq_len, MD_str);
    return MD_str;
}

static inline kstring_t *get_cs_from_digar(digar1_t *d, int n_digar, const char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end) {
    // Calculate cs tag from CIGAR string
    kstring_t *cs_str = (kstring_t*)calloc(1, sizeof(kstring_t));
    hts_pos_t pos; // Start position in reference sequence
    int prev_was_del = 0;
    for (int i = 0; i < n_digar; ++i) {
        int op = d[i].type, len = d[i].len;
        if (op == BAM_CEQUAL) {
            kputc(':', cs_str); kputw(len, cs_str);
        } else if (op == BAM_CDIFF) {
            for (int j = 0; j < len; ++j) {
                kputc('*', cs_str); 
                char ref_base = 'n';
                if (d[i].pos + j >= ref_beg && d[i].pos + j < ref_end) ref_base = tolower(ref_seq[d[i].pos - ref_beg + j]);
                kputc(ref_base, cs_str); kputc("acgtn"[d[i].alt_seq[j]], cs_str);
            }
        } else if (op == BAM_CINS) {
            kputc('+', cs_str);
            for (int j = 0; j < len; ++j) {
                char ins_base = "acgtn"[d[i].alt_seq[j]];
                kputc(ins_base, cs_str);
            }
        } else if (op == BAM_CDEL) {
            // Same reasoning as get_md_from_digar(): consecutive D ops are one
            // deletion, so "-cc" rather than "-c-c". NOTE: unverified against a
            // reference implementation, because neither bundled test file
            // carries a cs tag and update_bam1_tags() only rewrites cs when the
            // input already had one.
            if (!prev_was_del) kputc('-', cs_str);
            for (int j = 0; j < len; ++j) {
                pos = d[i].pos + j;
                char ref_base = 'n';
                if (pos >= ref_beg && pos < ref_end) ref_base = tolower(ref_seq[pos - ref_beg]);
                kputc(ref_base, cs_str);
            }
            prev_was_del = 1; continue;
        }
        prev_was_del = 0;
    }
    return cs_str;
}

// int strcicmp(char const *a, char const *b) {
//     for (;; a++, b++) {
//         int d = tolower((unsigned char)*a) - tolower((unsigned char)*b);
//         if (d != 0 || !*a)
//             return d;
//     }
// }

// collect size difference for MD and cs tags
int dry_update_bam1_tags(bam1_t *b, digar1_t *digar, int n_digar, const char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end) {
    int size_diff = 0;
    // check if MD tag exists
    uint8_t *md_tag = bam_aux_get(b, "MD");
    if (md_tag != NULL) {
        kstring_t *MD = get_md_from_digar(digar, n_digar, ref_seq, ref_beg, ref_end);
        if (strcmp(MD->s, (char*)bam_aux2Z(md_tag)) != 0) {
            int old_md_len = strlen((char*)bam_aux2Z(md_tag));
            size_diff += (MD->l + 1) - (old_md_len + 1);
        }
        free(MD->s); free(MD);
    }
    uint8_t *cs_tag = bam_aux_get(b, "cs");
    if (cs_tag != NULL) {
        kstring_t *cs = get_cs_from_digar(digar, n_digar, ref_seq, ref_beg, ref_end);
        if (strcmp(cs->s, (char*)bam_aux2Z(cs_tag)) != 0) {
            int old_cs_len = strlen((char*)bam_aux2Z(cs_tag));
            size_diff += (cs->l + 1) - (old_cs_len + 1);
        }
        free(cs->s); free(cs);
    }
    return size_diff;
}

int update_bam1_tags(bam1_t *b, digar1_t *digar, int n_digar, const char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end) {
    // Update NM/MD/cs tags based on the new CIGAR string
    // check if NM tag exists
    uint8_t *nm_tag = bam_aux_get(b, "NM");
    if (nm_tag != NULL) {
        int NM = get_nm_from_digar(digar, n_digar);
        if (NM != bam_aux2i(nm_tag)) {
            bam_aux_del(b, nm_tag);
            bam_aux_append(b, "NM", 'i', 4, (uint8_t *)&(NM));
        }
    }
    // check if MD tag exists
    uint8_t *md_tag = bam_aux_get(b, "MD");
    if (md_tag != NULL) {
        kstring_t *MD = get_md_from_digar(digar, n_digar, ref_seq, ref_beg, ref_end);
        if (strcmp(MD->s, (char*)bam_aux2Z(md_tag)) != 0) {
            // fprintf(stderr, "MD tag is not identical %s: %s != %s\n", bam_get_qname(b), MD->s, (char*)bam_aux2Z(md_tag));
            bam_aux_del(b, md_tag);
            bam_aux_append(b, "MD", 'Z', MD->l+1, (uint8_t*)MD->s);
        }
        free(MD->s); free(MD);
    }
    uint8_t *cs_tag = bam_aux_get(b, "cs");
    if (cs_tag != NULL) {
        kstring_t *cs = get_cs_from_digar(digar, n_digar, ref_seq, ref_beg, ref_end);
        if (strcmp(cs->s, (char*)bam_aux2Z(cs_tag)) != 0) {
            // fprintf(stderr, "CS tag is not identical %s: %s != %s\n", bam_get_qname(b), cs->s, (char*)bam_aux2Z(cs_tag));
            bam_aux_del(b, cs_tag);
            bam_aux_append(b, "cs", 'Z', cs->l+1, (uint8_t*)cs->s);
        }
        free(cs->s); free(cs);
    }
    return 0;
}

int refine_bam1(bam_chunk_t *chunk, int read_i, bam1_t *b) {
    // update b->core.pos based on chunk->digars[read_i]
    hts_pos_t new_pos = chunk->digars[read_i].digars[0].pos; // assume the first digar is the start position
    if (new_pos < 1) {
        fprintf(stderr, "Error: Invalid position %" PRId64 " for read %s in chunk %s:%" PRId64 "-%" PRId64 "\n",
                new_pos, bam_get_qname(b), chunk->tname, chunk->reg_beg, chunk->reg_end);
        return -1; // Invalid position
    }
    // fprintf(stderr, "refining read: %s %" PRId64 "\n", bam_get_qname(b), b->core.pos);
    if (new_pos != b->core.pos + 1) { // bam1_t->core.pos is 0-based, chunk->digars[read_i].digars[0].pos is 1-based
        b->core.pos = new_pos - 1; // Update to 0-based
    }
    // update b->cigar based on chunk->digars[read_i]
    uint32_t *new_cigar = (uint32_t*)malloc(chunk->digars[read_i].n_digar * sizeof(uint32_t));
    int new_n_cigar = 0, new_m_cigar = chunk->digars[read_i].n_digar;
    for (int i = 0; i < chunk->digars[read_i].n_digar; ++i) {
        new_cigar = longcalld_push_cigar(&new_n_cigar, &new_m_cigar, new_cigar,
                             chunk->digars[read_i].digars[i].type,
                             chunk->digars[read_i].digars[i].len);
    }
    // XXX currently, no supp alignment will be processed/generated
    // fix CLIP in CIGAR, change HARD_CLIP in primary CIGAR to SOFT_CLIP
    if (!(b->core.flag & BAM_FSUPPLEMENTARY)) {
        if ((new_cigar[0] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP) new_cigar[0] = (new_cigar[0] & ~BAM_CIGAR_MASK) | BAM_CSOFT_CLIP;
        if ((new_cigar[new_n_cigar-1] & BAM_CIGAR_MASK) == BAM_CHARD_CLIP) new_cigar[new_n_cigar-1] = (new_cigar[new_n_cigar-1] & ~BAM_CIGAR_MASK) | BAM_CSOFT_CLIP;
    } else {
        // Deliberately not refined. refine_bam1() rebuilds the CIGAR from the
        // read's digar, which describes the primary alignment: the H->S
        // conversion above is only valid for a primary, and a supplementary's
        // clipping must stay hard so the record remains a consistent member of
        // the read's alignment set. Supplementary records are therefore written
        // with their original CIGAR; the count is reported at end of run.
    }
    uint32_t *old_cigar = bam_get_cigar(b); int old_n_cigar = b->core.n_cigar;
    if (cigar_is_idential(old_cigar, old_n_cigar, new_cigar, new_n_cigar)) {
        free(new_cigar);
        return 0; // No need to update
    }
    // fprintf(stderr, "Refining BAM record: %s, old CIGAR: %d %d\n", bam_get_qname(b), new_n_cigar, old_n_cigar);
    // for (int i = 0; i < old_n_cigar; ++i) {
        // fprintf(stderr, "%d%c", old_cigar[i] >> BAM_CIGAR_SHIFT , "MIDNSHP=X"[bam_cigar_op(old_cigar[i])]);
    // } 
    // fprintf(stderr, "\nNew CIGAR: ");
    // for (int i = 0; i < new_n_cigar; ++i) {
        // fprintf(stderr, "%d%c", new_cigar[i] >> BAM_CIGAR_SHIFT, "MIDNSHP=X"[bam_cigar_op(new_cigar[i])]);
    // } fprintf(stderr, "\n");

    int old_data_len = b->l_data, new_data_len = b->l_data + (new_n_cigar-old_n_cigar) * sizeof(uint32_t);
    if (new_data_len > b->m_data) sam_realloc_bam_data(b, new_data_len);
    // fprintf(stderr, "old_l_data: %d, new_l_data: %d, %d\n", old_data_len, new_data_len, b->m_data);
    uint8_t *new_data = (uint8_t*)malloc(b->m_data);

    // Copy: QNAME
    int qname_len = b->core.l_qname;
    memcpy(new_data, b->data, qname_len);

    // Copy: new CIGAR
    memcpy(new_data + qname_len, new_cigar, new_n_cigar * sizeof(uint32_t));

    // Copy: SEQ + QUAL + AUX
    int rest_len = old_data_len - (qname_len + old_n_cigar * sizeof(uint32_t));
    // fprintf(stderr, "memcpy: new_size: %d, %d, new_offset: %d, old_size: %d, rest_len: %d\n", b->m_data, new_data_len, qname_len + new_n_cigar * sizeof(uint32_t), qname_len + old_n_cigar * sizeof(uint32_t), rest_len);
    memcpy(new_data + qname_len + new_n_cigar * sizeof(uint32_t),
           b->data + qname_len + old_n_cigar * sizeof(uint32_t), rest_len);

    free(b->data); free(new_cigar);
    b->data = new_data;
    b->core.n_cigar = new_n_cigar;
    b->l_data = new_data_len;
    // calculate NM/MD/cs tags
    // fprintf(stderr, "%s\n", bam_get_qname(b));
    update_bam1_tags(b, chunk->digars[read_i].digars, chunk->digars[read_i].n_digar, chunk->ref_seq, chunk->ref_beg, chunk->ref_end);
    return 0;
}

// ---------------------------------------------------------------------------
// aln_out_t: bounded reorder buffer for the phased alignment output
// ---------------------------------------------------------------------------

aln_out_t *aln_out_init(htsFile *fp, bam_hdr_t *hdr, int is_cram, int no_reorder, hts_pos_t window) {
    aln_out_t *o = (aln_out_t*)calloc(1, sizeof(aln_out_t));
    o->fp = fp; o->hdr = hdr; o->is_cram = is_cram != 0;
    o->no_reorder = no_reorder != 0;
    o->window = no_reorder ? 0 : window;
    o->m_buf = 256; o->buf = (bam1_t**)malloc(o->m_buf * sizeof(bam1_t*));
    o->n_buf = 0; o->last_tid = -1; o->last_pos = -1; o->n_written = 0;
    return o;
}

// Emit one record, checking that the file stays coordinate-sorted. "Silently
// invalid BAM" is not an acceptable outcome, so a violation is fatal.
static void aln_out_emit(aln_out_t *o, bam1_t *b, const char *tname) {
    // Unplaced records (tid < 0) sort after every placed one, so they are
    // exempt from the position check but must all come last.
    if (b->core.tid < 0) {
        o->seen_unplaced = 1;
    } else {
        if (!o->no_reorder && o->seen_unplaced) {
            _err_error_exit("Placed record %s at %s:%" PRId64 " written after an unmapped record; "
                            "unmapped records must come last.\n",
                            bam_get_qname(b), tname ? tname : "?", b->core.pos+1);
        }
        if (!o->no_reorder && o->last_tid >= 0 &&
            (b->core.tid < o->last_tid || (b->core.tid == o->last_tid && b->core.pos < o->last_pos))) {
            _err_error_exit("Output alignment went backwards at %s:%" PRId64 " %s (previous record was tid %d pos %" PRId64 "). "
                            "The reorder window (%" PRId64 " bp) was too small.\n",
                            tname ? tname : "?", b->core.pos+1, bam_get_qname(b), o->last_tid, o->last_pos+1, o->window);
        }
    }
    if (sam_write1(o->fp, o->hdr, b) < 0)
        _err_error_exit("Failed to write BAM record. %s:%" PRId64 " %s\n", tname ? tname : "?", b->core.pos+1, bam_get_qname(b));
    if (b->core.tid >= 0) { o->last_tid = b->core.tid; o->last_pos = b->core.pos; }
    o->n_written++;
}

// Drop the first n buffered records (already emitted).
static void aln_out_pop_front(aln_out_t *o, int n) {
    for (int i = 0; i < n; ++i) bam_destroy1(o->buf[i]);
    if (n < o->n_buf) memmove(o->buf, o->buf + n, (o->n_buf - n) * sizeof(bam1_t*));
    o->n_buf -= n;
}

void aln_out_flush(aln_out_t *o) {
    if (o == NULL) return;
    for (int i = 0; i < o->n_buf; ++i) aln_out_emit(o, o->buf[i], NULL);
    aln_out_pop_front(o, o->n_buf);
}

void aln_out_push(aln_out_t *o, bam1_t *b, const char *tname) {
    // Nothing can move, so there is nothing to reorder: write straight through
    // and avoid a copy per record. The sortedness check still applies.
    // Unplaced records never need reordering either.
    if (o->window == 0 || b->core.tid < 0) {
        if (b->core.tid < 0) aln_out_flush(o); // everything buffered is placed, so it precedes this
        aln_out_emit(o, b, tname); return;
    }

    // A new reference means everything buffered is already behind us.
    if (o->n_buf > 0 && o->buf[0]->core.tid != b->core.tid) aln_out_flush(o);

    if (o->n_buf == o->m_buf) {
        o->m_buf <<= 1;
        o->buf = (bam1_t**)realloc(o->buf, o->m_buf * sizeof(bam1_t*));
    }
    // Insertion sort from the back: the stream is nearly sorted, so this
    // almost always lands at the end after a single comparison.
    bam1_t *copy = bam_dup1(b);
    int i = o->n_buf;
    while (i > 0 && (o->buf[i-1]->core.tid > copy->core.tid ||
                     (o->buf[i-1]->core.tid == copy->core.tid && o->buf[i-1]->core.pos > copy->core.pos))) {
        o->buf[i] = o->buf[i-1]; --i;
    }
    o->buf[i] = copy; o->n_buf++;

    // Everything more than `window` behind the current input position is safe
    // to emit: no record still to come can start earlier than that.
    hts_pos_t safe = b->core.pos - o->window;
    int n_ready = 0;
    while (n_ready < o->n_buf && o->buf[n_ready]->core.tid == b->core.tid && o->buf[n_ready]->core.pos < safe) ++n_ready;
    for (int k = 0; k < n_ready; ++k) aln_out_emit(o, o->buf[k], tname);
    aln_out_pop_front(o, n_ready);
}

void aln_out_free(aln_out_t *o) {
    if (o == NULL) return;
    aln_out_flush(o);
    free(o->buf);
    free(o);
}

// Records the caller did not phase itself: secondary, unmapped, low-MQ, and
// supplementary alignments.
//
// A supplementary alignment is not an independent observation -- it is part of a
// read whose primary alignment we did phase -- so when its read's primary is in
// this same chunk it inherits that primary's HP/PS. Everything else has both
// tags stripped, since a stale HP/PS from the input would be a lie about our
// phasing. A supplementary whose primary lies in another chunk is not in the map
// and so gets nothing, which is the safe answer rather than a guess.
int write_untagged_read_to_bam(bam_chunk_t *chunk, aln_out_t *out, bam1_t *read, void *phase_by_qname) {
    int hap = 0; hts_pos_t ps = 0;
    if (phase_by_qname != NULL && (read->core.flag & BAM_FSUPPLEMENTARY) &&
        !(read->core.flag & (BAM_FSECONDARY | BAM_FUNMAP))) {
        khash_t(rphase) *h = (khash_t(rphase)*)phase_by_qname;
        khiter_t k = kh_get(rphase, h, bam_get_qname(read));
        if (k != kh_end(h)) { hap = kh_value(h, k).hap; ps = kh_value(h, k).ps; }
        out->n_unrefined_supp++; // reported at end of run when --refine-aln is on
    }

    uint8_t *hp = bam_aux_get(read, "HP");
    if (hp != NULL) bam_aux_del(read, hp);
    uint8_t *ps_tag = bam_aux_get(read, "PS");
    if (ps_tag != NULL) bam_aux_del(read, ps_tag);

    if (hap != 0 && ps > 0) { // both or neither, as everywhere else
        bam_aux_append(read, "HP", 'i', 4, (uint8_t *)&hap);
        bam_aux_append(read, "PS", 'i', 4, (uint8_t *)&ps);
    }
    aln_out_push(out, read, chunk->tname);
    return 0;
}

int write_processed_read_to_bam(const call_var_opt_t *opt, bam_chunk_t *chunk, aln_out_t *out, bam1_t *read, int read_i) {
    // The caller rebuilt read_i by re-reading the region and replaying
    // bam_record_is_callable() plus the n_up_ovlp_skip_reads[] counters. If any
    // of that drifts, we would silently stamp this record with another read's
    // haplotype, phase set and refined CIGAR. Verify identity instead.
    if (read_i < 0 || read_i >= chunk->n_reads) {
        _err_error_exit("read_i %d out of range [0, %d) when writing BAM record: %s:%" PRId64 " %s\n",
                        read_i, chunk->n_reads, chunk->tname, read->core.pos+1, bam_get_qname(read));
    }
    if (chunk->read_sigs[read_i] != bam_record_sig(read)) {
        _err_error_exit("BAM record/read_i mismatch at %s:%" PRId64 " %s (read_i %d): the phased-BAM "
                        "writer and the chunk loader disagree about which records are callable.\n",
                        chunk->tname, read->core.pos+1, bam_get_qname(read), read_i);
    }
    // refine bam cigars
    if (opt->refine_bam) {
        // check if bam and chunk->digars[read_i] match
        if (digar2qlen(chunk->digars+read_i) != read->core.l_qseq) {
            _err_error_exit("Read length mismatch when writing BAM record: %s:%" PRId64 " %s, digar_qlen: %d, bam_qlen: %d\n",
                            chunk->tname, read->core.pos+1, bam_get_qname(read),
                            digar2qlen(chunk->digars+read_i), read->core.l_qseq);
        }
        if (refine_bam1(chunk, read_i, read) != 0)
            _err_error_exit("Failed to refine BAM record: %s:%" PRId64 " %s\n", chunk->tname, read->core.pos+1, bam_get_qname(read));
    }
    // Add HP tag
    int hap = chunk->haps[read_i];
    hts_pos_t ps = chunk->phase_sets[read_i];

    // HP and PS are meaningless alone: a haplotype with no phase block to
    // anchor it says nothing, and a phase set with no haplotype assignment says
    // nothing either. Emit both or neither, decided once. Guarding them
    // independently is what let 6 of 587 ONT records out with exactly one.
    int emit_phase = (hap != 0 && ps > 0);
    if (emit_phase) {
        // check if HP tag exists
        uint8_t *hp = bam_aux_get(read, "HP");
        if (hp != NULL) {
            if (hap != bam_aux2i(hp)) {
                bam_aux_del(read, hp);
                bam_aux_append(read, "HP", 'i', 4, (uint8_t *)&(hap));
            }
        } else {
            bam_aux_append(read, "HP", 'i', 4, (uint8_t *)&(hap));
        }
    } else {
        // remove HP tag if we are not emitting a phase for this read
        uint8_t *hp = bam_aux_get(read, "HP");
        if (hp != NULL) bam_aux_del(read, hp);
    }
    // Add PS tag
    if (emit_phase) {
        // check if PS tag exists
        uint8_t *ps_tag = bam_aux_get(read, "PS");
        if (ps_tag != NULL) {
            if (ps != bam_aux2i(ps_tag)) {
                bam_aux_del(read, ps_tag);
                bam_aux_append(read, "PS", 'i', 4, (uint8_t *)&(ps));
            }
        } else {
            bam_aux_append(read, "PS", 'i', 4, (uint8_t *)&(ps));
        }
    } else {
        // remove PS tag if we are not emitting a phase for this read
        uint8_t *ps_tag = bam_aux_get(read, "PS");
        if (ps_tag != NULL) bam_aux_del(read, ps_tag);
    }
    // write to bam
    aln_out_push(out, read, chunk->tname);
    return 0;
}

// open & read bam, then output phased bam
int write_read_to_bam(bam_chunk_t *chunk, const struct call_var_opt_t *opt, const struct call_var_io_aux_t *io_aux) {
    int n_out_reads = 0;
    int min_mapq = opt->min_mq;
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end; int tid = chunk->tid;
    int n_bam = chunk->n_bam;

    // qname -> phase, over this chunk's phased primary reads, so a
    // supplementary alignment can inherit its read's HP/PS.
    khash_t(rphase) *phase_by_qname = NULL;
    if (chunk->read_names != NULL) {
        phase_by_qname = kh_init(rphase);
        for (int ri = 0; ri < chunk->n_reads; ++ri) {
            if (chunk->read_names[ri] == NULL) continue;
            if (chunk->haps[ri] == 0 || chunk->phase_sets[ri] <= 0) continue;
            int absent; khiter_t k = kh_put(rphase, phase_by_qname, chunk->read_names[ri], &absent);
            if (absent) {
                read_phase_t v; v.hap = chunk->haps[ri]; v.ps = chunk->phase_sets[ri];
                kh_value(phase_by_qname, k) = v;
            }
        }
    }

    // Emit records in coordinate order across all input files.
    //
    // Iterating one file to exhaustion before starting the next -- which is what
    // this used to do -- groups the output by input file, so position jumps
    // backwards at every file boundary. With two bundled halves as -X inputs that
    // was a measured 185 kb back-jump, producing a BAM that claimed
    // SO:coordinate and that samtools index rejected.
    //
    // chunk->reads[] is still filled in BAM-major order, so read_i is recovered
    // as bam_read_off[i] + (this file's callable ordinal). The per-read signature
    // check in write_processed_read_to_bam() verifies that arithmetic on every
    // record, which is what makes this safe to restructure.
    hts_itr_t **iters = (hts_itr_t**)calloc(n_bam, sizeof(hts_itr_t*));
    bam1_t **cur = (bam1_t**)calloc(n_bam, sizeof(bam1_t*));
    int *have = (int*)calloc(n_bam, sizeof(int));
    int *read_i = (int*)calloc(n_bam, sizeof(int));        // callable reads seen in this file
    int *skip_read_i = (int*)calloc(n_bam, sizeof(int));   // filtered-out reads seen in this file

    for (int i = 0; i < n_bam; ++i) {
        bam_hdr_t *header = opt->headers[i];
        iters[i] = sam_itr_queryi(opt->idxs[i], tid, reg_beg-1, reg_end);
        if (iters[i] == NULL) {
            _err_warning("Failed to create iterator for region: %s:%" PRId64 "-%" PRId64 ". Skipping.\n",
                         header->target_name[tid], reg_beg, reg_end);
            continue;
        }
        cur[i] = bam_init1();
        have[i] = (sam_itr_next(opt->bams[i], iters[i], cur[i]) >= 0);
    }

    for (;;) {
        // Pick the file whose pending record sorts first; ties go to the lower
        // file index, so the output is deterministic.
        int pick = -1;
        for (int i = 0; i < n_bam; ++i) {
            if (!have[i]) continue;
            if (pick < 0) { pick = i; continue; }
            if (cur[i]->core.tid < cur[pick]->core.tid ||
                (cur[i]->core.tid == cur[pick]->core.tid && cur[i]->core.pos < cur[pick]->core.pos)) pick = i;
        }
        if (pick < 0) break;

        bam1_t *read = cur[pick];
        // Which file this input's reads belong in. NULL means this sample has
        // no alignment output: the record is not written, but read_i still has
        // to advance over it, because it indexes chunk->haps for every sample.
        int smp = opt->bam_i_to_sample != NULL ? opt->bam_i_to_sample[pick] : 0;
        aln_out_t *out = (smp >= 0 && smp < LONGCALLD_MAX_SAMPLES) ? opt->aln_outs[smp] : NULL;

        if (!bam_record_is_callable(read, min_mapq)) {
            if (out != NULL && skip_read_i[pick] >= chunk->n_up_ovlp_skip_reads[pick]) {
                write_untagged_read_to_bam(chunk, out, read, phase_by_qname);
                n_out_reads++;
            }
            skip_read_i[pick]++;
        } else {
            if (out != NULL && read_i[pick] >= chunk->n_up_ovlp_reads[pick]) {
                write_processed_read_to_bam(opt, chunk, out, read,
                                            chunk->bam_read_off[pick] + read_i[pick]);
                n_out_reads++;
            }
            read_i[pick]++;
        }
        have[pick] = (sam_itr_next(opt->bams[pick], iters[pick], read) >= 0);
    }

    for (int i = 0; i < n_bam; ++i) {
        if (iters[i] != NULL) bam_itr_destroy(iters[i]);
        if (cur[i] != NULL) bam_destroy1(cur[i]);
    }
    free(iters); free(cur); free(have); free(read_i); free(skip_read_i);
    if (phase_by_qname != NULL) kh_destroy(rphase, phase_by_qname);
    return n_out_reads;
}

// check if multiple RG/SM tag are the same, if not same, output warning message
char *extract_sample_name_from_bam_header(bam_hdr_t *header) {
    int n_rg = sam_hdr_count_lines(header, "RG");
    char *sample_name = NULL;
    for (int i = 0; i < n_rg; ++i) {
        kstring_t ks = KS_INITIALIZE;
        if (sam_hdr_find_tag_pos(header, "RG", i, "SM", &ks) == 0) {
            if (i == 0) {
                sample_name = ks.s;
            } else {
                if (strcmp(sample_name, ks.s) != 0) {
                    free(ks.s);
                    _err_warning("Multiple RG/SM found in alignment file header, using the first one: %s\n", sample_name);
                    return sample_name;
                } free(ks.s);
            }
        }
    }
    if (sample_name != NULL) _err_info("Sample name from alignment file header: %s\n", sample_name);
    return sample_name;
}
