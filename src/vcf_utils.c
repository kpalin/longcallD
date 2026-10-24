#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "main.h"
#include "bam_utils.h"
#include "utils.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "collect_var.h"
#include "call_var_main.h"
#include "htslib/vcf.h"
#include "htslib/sam.h"

extern int LONGCALLD_VERBOSE;

void write_vcf_header(bam_hdr_t *hdr, struct call_var_opt_t *opt) {
    htsFile *out_vcf = opt->out_vcf;
    char *sample_name = opt->sample_name;
    bcf_hdr_t *vcf_hdr = bcf_hdr_init("w");
    if (!vcf_hdr) _err_error_exit("Could not allocate VCF header.\n");
    // File format
    // bcf_hdr_append(vcf_hdr, "##fileformat=VCFv4.3");

    // Get current date
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    char date[11];
    strftime(date, sizeof(date), "%Y%m%d", tm);
    char date_str[50];
    snprintf(date_str, sizeof(date_str), "##fileDate=%s", date);
    bcf_hdr_append(vcf_hdr, date_str);

    // Source information
    char source_str[100];
    snprintf(source_str, sizeof(source_str), "##source=%s version=%s", PROG, LONGCALLD_VERSION);
    bcf_hdr_append(vcf_hdr, source_str);

    // Command line
    char *cmd_str = (char*)_err_malloc(strlen(CMD) + 6);
    snprintf(cmd_str, strlen(CMD) + 6, "##CL=%s", CMD);
    bcf_hdr_append(vcf_hdr, cmd_str);
    free(cmd_str);

    // Reference sequence information
    for (int i = 0; i < hdr->n_targets; i++) {
        char contig_str[256];
        snprintf(contig_str, sizeof(contig_str), "##contig=<ID=%s,length=%d>", hdr->target_name[i], hdr->target_len[i]);
        bcf_hdr_append(vcf_hdr, contig_str);
    }

    // FILTER fields
    bcf_hdr_append(vcf_hdr, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
    bcf_hdr_append(vcf_hdr, "##FILTER=<ID=LowQual,Description=\"Low quality variant\">");
    bcf_hdr_append(vcf_hdr, "##FILTER=<ID=RefCall,Description=\"Reference call\">");
    bcf_hdr_append(vcf_hdr, "##FILTER=<ID=NoCall,Description=\"Site has depth=0 resulting in no call\">");
    if (opt->n_samples > 1) {
        bcf_hdr_append(vcf_hdr, "##FILTER=<ID=NormalOnly,Description=\"Present in the normal but not the tumour; treated as false\">");
        bcf_hdr_append(vcf_hdr, "##FILTER=<ID=LowDP,Description=\"Too little depth in one sample to classify\">");
        bcf_hdr_append(vcf_hdr, "##FILTER=<ID=Ambiguous,Description=\"Evidence fits no tumour/normal class cleanly\">");
        bcf_hdr_append(vcf_hdr, "##INFO=<ID=UNPHASED,Number=0,Type=Flag,Description=\"Classified on pooled per-sample counts; no usable haplotype grid at this site\">");
    }

    // INFO fields
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of the variant described in this record\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic/mosaic variant\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=CLEAN,Number=0,Type=Flag,Description=\"Clean-region variant (SNP or simple indel in non-repetitive region)\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variant\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=SVLEN,Number=A,Type=Integer,Description=\"Difference in length between REF and ALT alleles\">");
    // TSD info
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=TSD,Number=A,Type=String,Description=\"Target site duplication sequence\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=TSDLEN,Number=A,Type=Integer,Description=\"Length of target site duplication\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=POLYALEN,Number=A,Type=Integer,Description=\"Length of polyA/T sequence\">");
    // there will be 2 TSDs for DEL
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=MEI,Number=0,Type=Flag,Description=\"Mobile element insertion\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=TSDPOS1,Number=A,Type=Integer,Description=\"Start position of first target site duplication on CHROM\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=TSDPOS2,Number=A,Type=Integer,Description=\"Start position of second target site duplication on CHROM\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=REPNAME,Number=A,Type=String,Description=\"Repeat name\">");
    if (opt->n_samples > 1) {
        bcf_hdr_append(vcf_hdr, "##INFO=<ID=NORMALONLY,Number=0,Type=Flag,Description=\"Seen in the normal and not the tumour\">");
        bcf_hdr_append(vcf_hdr, "##INFO=<ID=LOH,Number=0,Type=Flag,Description=\"Normal-only because the tumour has lost that chromosome copy\">");
        bcf_hdr_append(vcf_hdr, "##INFO=<ID=SOMHAP,Number=1,Type=Integer,Description=\"Tumour haplotype carrying the somatic alt (1 or 2)\">");
        bcf_hdr_append(vcf_hdr, "##INFO=<ID=HAPVAF,Number=1,Type=Float,Description=\"Alt fraction within that haplotype; estimates purity x cancer-cell fraction, unlike the pooled VAF which reads purity/2 for a clonal het\">");
        bcf_hdr_append(vcf_hdr, "##INFO=<ID=SPV,Number=1,Type=Float,Description=\"Fisher exact p, the somatic slot against the other three\">");
    }
    // bcf_hdr_append(vcf_hdr, "##INFO=<ID=REPCLSFAM,Number=1,Type=String,Description=\"Repeat class/family\">");

    // FORMAT fields
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for each allele\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant allele frequency\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods rounded to the closest integer\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">");
    // per-sample variant supporting read IDs
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=ALTREADS,Number=.,Type=String,Description=\"IDs of reads supporting the variant\">");
    if (opt->n_samples > 1) {
        bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=HDP,Number=2,Type=Integer,Description=\"Depth on haplotype 1 and 2\">");
        bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=HAD,Number=2,Type=Integer,Description=\"Alt depth on haplotype 1 and 2\">");
    }

    // Add samples. Across both columns, HDP/HAD are the four-slot grid the
    // paired verdict was decided on, so a disputed call can be re-litigated
    // from the VCF alone rather than by re-piling the BAM.
    if (opt->n_samples > 1) {
        // NORMAL first, then TUMOR: the GATK/Strelka convention.
        bcf_hdr_add_sample(vcf_hdr, opt->sample_names != NULL && opt->sample_names[LONGCALLD_SAMPLE_NORMAL] != NULL
                                    ? opt->sample_names[LONGCALLD_SAMPLE_NORMAL] : "NORMAL");
        bcf_hdr_add_sample(vcf_hdr, opt->sample_names != NULL && opt->sample_names[LONGCALLD_SAMPLE_TUMOR] != NULL
                                    ? opt->sample_names[LONGCALLD_SAMPLE_TUMOR] : "TUMOR");
    } else {
        bcf_hdr_add_sample(vcf_hdr, sample_name);
    }
    bcf_hdr_add_sample(vcf_hdr, NULL); // Finalize sample list

    // Write header to file
    if (bcf_hdr_write(out_vcf, vcf_hdr) < 0) {
        fprintf(stderr, "Error: Could not write VCF header.\n");
    }
    opt->vcf_hdr = vcf_hdr;
}

// FILTER column for a paired call. More than one can apply, so they are
// semicolon-joined the way the VCF spec intends.
static void append_paired_filter(kstring_t *ks, const var1_t *v) {
    int n = 0;
    switch (v->paired_class) {
        case LONGCALLD_PAIRED_NORMAL_ONLY:
        case LONGCALLD_PAIRED_LOH:   kputs("NormalOnly", ks); n++; break;
        case LONGCALLD_PAIRED_LOWDP: kputs("LowDP", ks); n++; break;
        case LONGCALLD_PAIRED_AMBIG: kputs("Ambiguous", ks); n++; break;
        default: break; // SOMATIC and GERMLINE pass
    }
    if (n == 0) kputs("PASS", ks);
    // Note: a call made without a usable haplotype grid is reported through
    // INFO/UNPHASED, not here. As a FILTER it made every such call non-PASS,
    // which means som.py -- and any `bcftools view -f PASS` consumer -- drops it
    // entirely, so the two evidence regimes could not be compared at all. In
    // INFO it is scored and still stratifiable:
    //   bcftools view -f PASS -i 'SOMATIC=1 && UNPHASED=1'
}

// One sample's FORMAT fields, from the slot grid.
//
// Genotype is read off the evidence rather than inherited from the pooled
// consensus: in a tumour/normal pair the two samples genuinely can differ, and
// reporting the pooled genotype twice would hide exactly what the caller is
// for.
static void append_paired_sample(kstring_t *ks, const struct call_var_opt_t *opt,
                                 const var1_t *v, int smp, int want_ps) {
    const var_slot_counts_t *sc = &v->slots;
    int dp = smp_dp(sc, smp), alt = smp_alt(sc, smp);
    int d1 = slot_dp(sc, smp, 1), d2 = slot_dp(sc, smp, 2);
    int min_hap_dp = (smp == LONGCALLD_SAMPLE_NORMAL) ? opt->min_normal_hap_dp : opt->min_tumor_hap_dp;
    int phased = (d1 >= min_hap_dp && d2 >= min_hap_dp);

    if (dp == 0) {
        kputs("./.", ks);
    } else if (phased) {
        kputw(sc->alt[smp][1] > opt->max_hap_alt ? 1 : 0, ks);
        kputc('|', ks);
        kputw(sc->alt[smp][2] > opt->max_hap_alt ? 1 : 0, ks);
    } else {
        double af = (double)alt / dp;
        if (af > opt->max_af)      kputs("1/1", ks);
        else if (af >= opt->min_af) kputs("0/1", ks);
        else                        kputs("0/0", ks);
    }
    ksprintf(ks, ":%d:%d,%d:%.3f:%d", dp, dp - alt, alt, dp > 0 ? (double)alt / dp : 0.0, v->GQ);
    if (want_ps) ksprintf(ks, ":%" PRId64, v->PS);
    ksprintf(ks, ":%d,%d:%d,%d", d1, d2, sc->alt[smp][1], sc->alt[smp][2]);
}

int write_var_to_vcf(var_t *vars, const struct call_var_opt_t *opt, bam_chunk_t *chunk) {
    char *chrom = chunk->tname;
    htsFile *out_vcf = opt->out_vcf;  // htsFile pointer
    int n_vars = vars->n;
    int n_output_vars = 0;
    int buf_m = 50000; char *buffer = (char*)_err_malloc(buf_m * sizeof(char));  // Buffer for output string
    int len = 0;        // Length of the output string

    for (int var_i = 0; var_i < n_vars; var_i++) {
        var1_t var = vars->vars[var_i];
        if (var.n_alt_allele == 0) continue;
        if (var.DP < opt->min_dp) continue;
        if (opt->out_somatic && var.is_somatic) {
            if (var.AD[1] < opt->min_somatic_te_dp) continue;
            else if (var.AD[1] < opt->min_alt_dp && var.tsd_len <= 0) continue; 
        } else {
            if (var.AD[1] < opt->min_alt_dp) continue;
        }

        // Validate bases
        if (opt->out_amb_base == 0) {
            uint8_t skip = 0;
            for (int j = 0; j < var.ref_len; j++) {
                if (var.ref_bases[j] >= 4) {
                    if (LONGCALLD_VERBOSE >= 2) 
                        fprintf(stderr, "Invalid ref base: %s %" PRId64 " %d\n", chrom, var.pos + j, var.ref_bases[j]);
                    skip = 1; break;
                }
            }
            if (skip) continue;

            for (int j = 0; j < var.n_alt_allele; j++) {
                for (int k = 0; k < var.alt_len[j]; k++) {
                    if (var.alt_bases[j][k] >= 4) {
                        if (LONGCALLD_VERBOSE >= 2) 
                            fprintf(stderr, "Invalid alt base: %s %" PRId64 " %d\n", chrom, var.pos, var.alt_bases[j][k]);
                        skip = 1; break;
                    }
                }
                if (skip) break;
            }
            if (skip) continue;
        }

        // Write to buffer
        len = snprintf(buffer, buf_m, "%s\t%" PRId64 "\t.\t", chrom, var.pos);

        // Write REF
        int base_len = var.ref_len;
        for (int j = 0; j < var.n_alt_allele; j++) base_len += var.alt_len[j];
        if (base_len > buf_m) {
            buf_m = base_len + 10000;
            buffer = (char*)_err_realloc(buffer, buf_m * sizeof(char));
        }
        for (int j = 0; j < var.ref_len; j++) 
            len += snprintf(buffer + len, buf_m - len, "%c", "ACGTN"[var.ref_bases[j]]);

        // Write ALT
        len += snprintf(buffer + len, buf_m - len, "\t");
        for (int j = 0; j < var.n_alt_allele; j++) {
            for (int k = 0; k < var.alt_len[j]; k++) 
                len += snprintf(buffer + len, buf_m - len, "%c", "ACGTN"[var.alt_bases[j][k]]);
            if (j < var.n_alt_allele - 1) len += snprintf(buffer + len, buf_m - len, ",");
        }

        // Structural Variant (SV) annotation
        int k = 0;
        char SVLEN[1024] = "SVLEN=", tmp[1024];
        char SVTYPE[1024] = "SVTYPE=";
        for (int i = 0; i < var.n_alt_allele; i++) { // assert(var.n_alt_allele == 1)
            if (var.is_sv) { 
                if (k > 0) {
                    strcat(SVLEN, ",");
                    strcat(SVTYPE, ",");
                }
                sprintf(tmp, "%d", var.alt_len[i] - var.ref_len);
                strcat(SVLEN, tmp);
                sprintf(tmp, "%s", var.alt_len[i] > var.ref_len ? "INS" : "DEL");
                strcat(SVTYPE, tmp);
                k++;
            }
        }
        
        // Write QUAL, FILTER, INFO
        if (opt->n_samples > 1) {
            kstring_t fks = KS_INITIALIZE;
            append_paired_filter(&fks, &var);
            len += snprintf(buffer + len, buf_m - len, "\t%d\t%s\t", var.QUAL, fks.s ? fks.s : "PASS");
            ks_free(&fks);
        } else {
            len += snprintf(buffer + len, buf_m - len, "\t%d\tPASS\t", var.QUAL);
        }
        if (var.is_clean) len += snprintf(buffer + len, buf_m - len, "CLEAN;");
        if (opt->n_samples > 1) {
            if (var.paired_class == LONGCALLD_PAIRED_SOMATIC) len += snprintf(buffer + len, buf_m - len, "SOMATIC;");
            if (var.paired_class == LONGCALLD_PAIRED_NORMAL_ONLY || var.paired_class == LONGCALLD_PAIRED_LOH)
                len += snprintf(buffer + len, buf_m - len, "NORMALONLY;");
            if (var.paired_class == LONGCALLD_PAIRED_LOH) len += snprintf(buffer + len, buf_m - len, "LOH;");
            if (var.is_unphased_call) len += snprintf(buffer + len, buf_m - len, "UNPHASED;");
        } else if (var.is_somatic) len += snprintf(buffer + len, buf_m - len, "SOMATIC;");
        if (var.te_seq_i >= 0) len += snprintf(buffer + len, buf_m - len, "MEI;");
        len += snprintf(buffer + len, buf_m - len, "END=%" PRId64 "", var.pos + var.ref_len - 1);
        if (opt->n_samples > 1 && var.paired_class == LONGCALLD_PAIRED_SOMATIC) {
            if (var.som_hap > 0) len += snprintf(buffer + len, buf_m - len, ";SOMHAP=%d", var.som_hap);
            len += snprintf(buffer + len, buf_m - len, ";HAPVAF=%.4f;SPV=%.3g", var.hap_vaf, var.paired_pval);
        }
        if (var.is_sv) { 
            len += snprintf(buffer + len, buf_m - len, ";%s;%s", SVTYPE, SVLEN);
            if (var.tsd_len > 0) {
                len += snprintf(buffer + len, buf_m - len, ";TSD=");
                for (int i = 0; i < var.tsd_len; ++i) len += snprintf(buffer + len, buf_m - len, "%c", "ACGTN"[var.tsd_seq[i]]);
                len += snprintf(buffer + len, buf_m - len, ";TSDLEN=%d;POLYALEN=%d;TSDPOS1=%" PRId64 "", var.tsd_len, var.polya_len, var.tsd_pos1);
                if (var.tsd_pos2 > 0) len += snprintf(buffer + len, buf_m - len, ";TSDPOS2=%" PRId64 "", var.tsd_pos2);
            }
            if (var.te_seq_i >= 0) len += snprintf(buffer + len, buf_m - len, ";REPNAME=%c%s", "+-"[var.te_is_rev], opt->te_seq_names[var.te_seq_i]);
        }
        len += snprintf(buffer + len, buf_m - len, "\t");

        // Write FORMAT and Genotype Data
        if (opt->n_samples > 1) {
            // Two sample columns, NORMAL then TUMOR. Built in a kstring so the
            // doubled per-sample block cannot overrun the record buffer: the
            // old hand-rolled growth check only accounted for REF+ALT.
            int want_ps = (var.PS != 0);
            kstring_t sks = KS_INITIALIZE;
            kputs("GT:DP:AD:VAF:GQ", &sks);
            if (want_ps) kputs(":PS", &sks);
            kputs(":HDP:HAD", &sks);
            kputc('\t', &sks);
            append_paired_sample(&sks, opt, &var, LONGCALLD_SAMPLE_NORMAL, want_ps);
            kputc('\t', &sks);
            append_paired_sample(&sks, opt, &var, LONGCALLD_SAMPLE_TUMOR, want_ps);
            if (len + (int)sks.l + 2 >= buf_m) {
                buf_m = len + (int)sks.l + 1024;
                buffer = (char*)_err_realloc(buffer, buf_m);
            }
            len += snprintf(buffer + len, buf_m - len, "%s", sks.s);
            ks_free(&sks);
        } else {
            int gt1 = var.GT[0], gt2 = var.GT[1];
            int is_hom = gt1 == gt2; int gt_seperator = '|';
            if (var.PS == 0) {
                gt_seperator = '/';
                if (gt1 > gt2) { int tmp = gt1; gt1 = gt2; gt2 = tmp;  }
            }
            // add SVREADS to FORMAT if outputting supporting read names for SVs
            len += snprintf(buffer + len, buf_m - len, "GT:DP:AD:VAF:GQ");
            if (is_hom == 0 && var.PS != 0) len += snprintf(buffer + len, buf_m - len, ":PS");
            if (opt->output_var_rnames || (var.is_sv && (opt->output_sv_rnames)) || (var.is_somatic && opt->output_somatic_var_rnames)) {
                len += snprintf(buffer + len, buf_m - len, ":ALTREADS");
            }
            len += snprintf(buffer + len, buf_m - len, "\t");

            // GT:DP
            len += snprintf(buffer + len, buf_m - len, "%d%c%d:%d:", gt1, gt_seperator, gt2, var.DP);
            // AD
            for (int j = 0; j < 1 + var.n_alt_allele; j++) {
                if (j > 0) len += snprintf(buffer + len, buf_m - len, ",");
                len += snprintf(buffer + len, buf_m - len, "%d", var.AD[j]);
            }
            // VAF
            for (int j = 0; j < var.n_alt_allele; j++) {
                if (j == 0) len += snprintf(buffer + len, buf_m - len, ":");
                if (j > 0) len += snprintf(buffer + len, buf_m - len, ",");
                float vaf = (float)var.AD[j+1] / var.DP;
                len += snprintf(buffer + len, buf_m - len, "%.3f", vaf);
            }
            // GQ
            len += snprintf(buffer + len, buf_m - len, ":%d", var.GQ);
            // PS
            if (is_hom == 0 && var.PS != 0) 
                len += snprintf(buffer + len, buf_m - len, ":%" PRId64 "", var.PS);
            // ALTREADS: output_var_rnames for all variants, output_sv_rnames for SVs, output_somatic_var_rnames for somatic SVs
            // if (var.is_sv && (opt->output_sv_rnames || (opt->output_somatic_var_rnames && var.is_somatic))) {
            if (opt->output_var_rnames || (var.is_sv && (opt->output_sv_rnames)) || (var.is_somatic && opt->output_somatic_var_rnames)) {
                // check if the buffer is large enough for read names
                int n_reads_needed = 0; // 
                for (int i = 0; i < var.AD[1]; ++i)
                    n_reads_needed += strlen(chunk->read_names[var.alt_read_i[i]]) + 1;
                if (len + n_reads_needed >= buf_m) {
                    buf_m = len + n_reads_needed + 1024;
                    buffer = (char*)_err_realloc(buffer, buf_m);
                }
                len += snprintf(buffer + len, buf_m - len, ":");
                if (var.AD[1] > 0) {
                    for (int i = 0; i < var.AD[1]; ++i) {
                        if (i > 0) len += snprintf(buffer + len, buf_m - len, ",");
                        int alt_read_i = var.alt_read_i[i];
                        len += snprintf(buffer + len, buf_m - len, "%s", chunk->read_names[alt_read_i]);
                    }
                } else { // for SVs without read support, output "." for SVREADS
                    len += snprintf(buffer + len, buf_m - len, ".");
                }
            }
        }
        len += snprintf(buffer + len, buf_m - len, "\n");

        // Write to htsFile
        if (out_vcf->format.compression!=no_compression) {
            if (bgzf_write(out_vcf->fp.bgzf, buffer, len) < 0) {
                _err_error_exit("Could not write to VCF file.\n");
            }
        } else {
            if (hwrite(out_vcf->fp.hfile, buffer, len) < 0) {
                _err_error_exit("Could not write to VCF file.\n");
            }
        }
        n_output_vars++;
    }
    free(buffer);
    return n_output_vars;
}
