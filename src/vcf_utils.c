#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
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
    {
        char lq[256];
        snprintf(lq, sizeof(lq), "##FILTER=<ID=LowQual,Description=\"QUAL below --min-qual %d\">", opt->min_qual);
        bcf_hdr_append(vcf_hdr, lq);
    }
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
        bcf_hdr_append(vcf_hdr, "##INFO=<ID=TVAF,Number=1,Type=Float,Description=\"Alt fraction over the pooled tumour depth, reported in place of HAPVAF on UNPHASED calls where no haplotype was resolved; reads purity/2 for a clonal het\">");
        bcf_hdr_append(vcf_hdr, "##INFO=<ID=SPV,Number=1,Type=Float,Description=\"Fisher exact p, the somatic slot against the other three\">");
    }
    // bcf_hdr_append(vcf_hdr, "##INFO=<ID=REPCLSFAM,Number=1,Type=String,Description=\"Repeat class/family\">");
    // QUAL has no header line of its own, so what it means is stated here.
    bcf_hdr_append(vcf_hdr, "##longcallD_QUAL=\"-10 log10 P(false), capped at 999. For a SOMATIC record, P(not somatic: absent, or germline); otherwise P(no variant). P(false) is a per-technology logistic model over the read likelihood (INFO/QLIK: per-read base quality, mapping quality and indel repeat context) and features read likelihoods cannot see: strand bias, ref/alt mapping- and base-quality ranks, repeat context, unphased and mis-phased reads, and depth. Where no model is set, QUAL is QLIK\"");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=QLIK,Number=1,Type=Float,Description=\"QUAL from the read likelihoods alone, before the artifact model\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=SOR,Number=1,Type=Float,Description=\"Strand odds ratio (GATK definition) over the alt-carrying sample\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=MQRS,Number=1,Type=Float,Description=\"Rank-sum z of alt against ref read mapping quality; negative: alt lower\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=BQRS,Number=1,Type=Float,Description=\"Rank-sum z of alt against ref base quality; negative: alt lower\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=RU,Number=1,Type=String,Description=\"Repeat unit of the homopolymer or tandem repeat at the variant\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=RPL,Number=1,Type=Integer,Description=\"Length in reference bases of that repeat tract\">");
    if (getenv("LONGCALLD_QUALCHK") != NULL)
        bcf_hdr_append(vcf_hdr, "##INFO=<ID=QF,Number=.,Type=Float,Description=\"Debug: quality class, model P(false) (0 where none is set), then the model features in var_qual.h order\">");

    // FORMAT fields
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality: this sample's PL of the best other genotype minus that of the one printed, capped at 99; 0 when the printed genotype is not the most likely\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for each allele\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=VAF,Number=A,Type=Float,Description=\"Variant allele frequency\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods rounded to the closest integer, from per-read base and mapping quality; a tumour's (or a mosaic call's) 0/1 is the best alt fraction up to 1/2, not exactly 1/2\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=UPR,Number=2,Type=Integer,Description=\"Unphased reads: of all reads, and of alt reads, those with no haplotype in the site's phase block\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=UPF,Number=2,Type=Float,Description=\"Unphased proportion: UPR over this sample's DP, and over its alt reads\">");
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
static void append_paired_filter(kstring_t *ks, const struct call_var_opt_t *opt, const var1_t *v) {
    int n = 0;
    switch (v->paired_class) {
        case LONGCALLD_PAIRED_NORMAL_ONLY:
        case LONGCALLD_PAIRED_LOH:   kputs("NormalOnly", ks); n++; break;
        case LONGCALLD_PAIRED_LOWDP: kputs("LowDP", ks); n++; break;
        case LONGCALLD_PAIRED_AMBIG: kputs("Ambiguous", ks); n++; break;
        default: break; // SOMATIC and GERMLINE pass
    }
    if (opt->min_qual > 0 && v->QUAL < opt->min_qual) { if (n++) kputc(';', ks); kputs("LowQual", ks); }
    if (n == 0) kputs("PASS", ks);
    // Note: a call made without a usable haplotype grid is reported through
    // INFO/UNPHASED, not here. As a FILTER it made every such call non-PASS,
    // which means som.py -- and any `bcftools view -f PASS` consumer -- drops it
    // entirely, so the two evidence regimes could not be compared at all. In
    // INFO it is scored and still stratifiable:
    //   bcftools view -f PASS -i 'SOMATIC=1 && UNPHASED=1'
}

// INFO/UNPHASED says a VERDICT was reached on pooled per-sample counts because
// the haplotype grid was not usable, so it belongs on the records that carry a
// verdict. A germline or ambiguous record has no verdict for it to qualify, and
// that is where almost all of it landed: 129732 of the 131189 occurrences on the
// HiFi HCC1395 chr1 pair, 98.9%. Those tags also put the flag at odds with the
// genotypes beside it, because the GT separator is decided per sample from that
// sample's own slot depths while the flag is decided once for the site.
static int paired_call_is_unphased(const var1_t *v) {
    if (!v->is_unphased_call) return 0;
    return v->paired_class == LONGCALLD_PAIRED_SOMATIC ||
           v->paired_class == LONGCALLD_PAIRED_NORMAL_ONLY ||
           v->paired_class == LONGCALLD_PAIRED_LOH;
}

// FORMAT UPR and UPF: how much of one sample's evidence the phasing could not
// place on a haplotype of the site's block. Read from the quality model's own
// tallies, which use the slot grid's rule, so UPR[0] == DP - HDP[0] - HDP[1].
static void append_unphased_fields(kstring_t *ks, const var1_t *v, int smp) {
    int dp = v->vq.n_dp[smp], u = v->vq.n_unph[smp], ua = v->vq.n_unph_alt[smp];
    int alt = smp_alt(&v->slots, smp);
    ksprintf(ks, ":%d,%d:%.3f,%.3f", u, ua, dp > 0 ? (double)u / dp : 0.0, alt > 0 ? (double)ua / alt : 0.0);
}

// INFO fields behind QUAL.
static void append_qual_info(kstring_t *ks, const var1_t *v) {
    const var_qual_t *q = &v->vq;
    if (!q->is_set) return;
    ksprintf(ks, ";QLIK=%.1f;SOR=%.3f;MQRS=%.2f;BQRS=%.2f", q->qlik, q->sor, q->mq_rs, q->bq_rs);
    if (q->rep_period > 0 && q->rep_unit[0] != '\0') ksprintf(ks, ";RU=%s;RPL=%d", q->rep_unit, q->rep_len);
    if (getenv("LONGCALLD_QUALCHK") != NULL) {
        ksprintf(ks, ";QF=%d,%.4g", q->vq_class, q->p_art);
        for (int f = 0; f < VQF_N; ++f) ksprintf(ks, ",%.4g", q->feat[f]);
    }
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
    // A phased genotype has to say which phase block it is phased in: "|" is a
    // claim about haplotype 1 and haplotype 2 of some block, and without a PS
    // to name that block the claim cannot be read. So want_ps gates the
    // separator, which is also what keeps INFO/UNPHASED and the separator
    // answering the same question -- see paired_call_is_unphased().
    int phased = (want_ps && d1 >= min_hap_dp && d2 >= min_hap_dp);

    int dosage = -1; // alt alleles in the printed genotype; -1 for "./."
    if (dp == 0) {
        kputs("./.", ks);
    } else if (phased) {
        int g1 = paired_slot_is_alt(opt, sc, smp, 1) ? 1 : 0, g2 = paired_slot_is_alt(opt, sc, smp, 2) ? 1 : 0;
        kputw(g1, ks); kputc('|', ks); kputw(g2, ks);
        dosage = g1 + g2;
    } else {
        double af = (double)alt / dp;
        if (af > opt->max_af)      { kputs("1/1", ks); dosage = 2; }
        else if (af >= opt->min_af) { kputs("0/1", ks); dosage = 1; }
        else                        { kputs("0/0", ks); dosage = 0; }
    }
    // GQ is this sample's own: the genotype printed here is decided per sample,
    // so the confidence in it has to be too.
    const int16_t *pl = v->vq.pl[smp];
    ksprintf(ks, ":%d:%d,%d:%.3f:%d:%d,%d,%d", dp, dp - alt, alt, dp > 0 ? (double)alt / dp : 0.0,
             dosage >= 0 ? var_qual_gq(pl, dosage) : 0, pl[0], pl[1], pl[2]);
    if (want_ps) ksprintf(ks, ":%" PRId64, v->PS);
    ksprintf(ks, ":%d,%d:%d,%d", d1, d2, sc->alt[smp][1], sc->alt[smp][2]);
    append_unphased_fields(ks, v, smp);
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
            append_paired_filter(&fks, opt, &var);
            len += snprintf(buffer + len, buf_m - len, "\t%d\t%s\t", var.QUAL, fks.s ? fks.s : "PASS");
            ks_free(&fks);
        } else {
            len += snprintf(buffer + len, buf_m - len, "\t%d\t%s\t", var.QUAL,
                            (opt->min_qual > 0 && var.QUAL < opt->min_qual) ? "LowQual" : "PASS");
        }
        if (var.is_clean) len += snprintf(buffer + len, buf_m - len, "CLEAN;");
        if (opt->n_samples > 1) {
            if (var.paired_class == LONGCALLD_PAIRED_SOMATIC) len += snprintf(buffer + len, buf_m - len, "SOMATIC;");
            if (var.paired_class == LONGCALLD_PAIRED_NORMAL_ONLY || var.paired_class == LONGCALLD_PAIRED_LOH)
                len += snprintf(buffer + len, buf_m - len, "NORMALONLY;");
            if (var.paired_class == LONGCALLD_PAIRED_LOH) len += snprintf(buffer + len, buf_m - len, "LOH;");
            if (paired_call_is_unphased(&var)) len += snprintf(buffer + len, buf_m - len, "UNPHASED;");
        } else if (var.is_somatic) len += snprintf(buffer + len, buf_m - len, "SOMATIC;");
        if (var.te_seq_i >= 0) len += snprintf(buffer + len, buf_m - len, "MEI;");
        len += snprintf(buffer + len, buf_m - len, "END=%" PRId64 "", var.pos + var.ref_len - 1);
        if (opt->n_samples > 1 && var.paired_class == LONGCALLD_PAIRED_SOMATIC) {
            if (var.som_hap > 0) {
                len += snprintf(buffer + len, buf_m - len, ";SOMHAP=%d;HAPVAF=%.4f", var.som_hap, var.hap_vaf);
            } else {
                // An unphased call has no haplotype, so it has no within-haplotype
                // fraction either: classify_paired_vars() fills hap_vaf from the
                // pooled tumour depth there. That is half the HAPVAF of the same
                // clonal het variant, and emitting it under HAPVAF contradicted
                // that key's own header and silently mixed two scales in whatever
                // read it -- notably the purity curve, where HAPVAF is supposed to
                // track purity directly while the pooled VAF tracks purity/2.
                len += snprintf(buffer + len, buf_m - len, ";TVAF=%.4f", var.hap_vaf);
            }
            len += snprintf(buffer + len, buf_m - len, ";SPV=%.3g", var.paired_pval);
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
        {
            kstring_t qks = KS_INITIALIZE;
            append_qual_info(&qks, &var);
            if (qks.l > 0) {
                if (len + (int)qks.l + 2 >= buf_m) { buf_m = len + (int)qks.l + 1024; buffer = (char*)_err_realloc(buffer, buf_m); }
                len += snprintf(buffer + len, buf_m - len, "%s", qks.s);
            }
            ks_free(&qks);
        }
        len += snprintf(buffer + len, buf_m - len, "\t");

        // Write FORMAT and Genotype Data
        if (opt->n_samples > 1) {
            // Two sample columns, NORMAL then TUMOR. Built in a kstring so the
            // doubled per-sample block cannot overrun the record buffer: the
            // old hand-rolled growth check only accounted for REF+ALT.
            // Absent, not zero: a somatic variant drives no phasing of its own
            // and carries assign_hap.c's -1 sentinel until the block it was
            // decided on is filled in. Testing != 0 printed that sentinel.
            // A call the grid did not decide reports no phase block, as
            // FORMAT/PS is documented to do. The somatic backfill below already
            // honoured that, but a variant that inherited a phase set from the
            // germline chain kept it regardless, which is how 1069 of the 1425
            // pooled-decided somatic calls on the HiFi HCC1395 chr1 pair named a
            // block they were not decided on.
            int want_ps = (var.PS > 0 && !paired_call_is_unphased(&var));
            kstring_t sks = KS_INITIALIZE;
            kputs("GT:DP:AD:VAF:GQ:PL", &sks);
            if (want_ps) kputs(":PS", &sks);
            kputs(":HDP:HAD:UPR:UPF", &sks);
            kputc('\t', &sks);
            append_paired_sample(&sks, opt, &var, LONGCALLD_SAMPLE_NORMAL, want_ps);
            kputc('\t', &sks);
            append_paired_sample(&sks, opt, &var, LONGCALLD_SAMPLE_TUMOR, want_ps);
            // The two ways this record can describe its own phasing must not
            // contradict each other: if INFO/UNPHASED says no haplotype grid
            // decided the call, no sample may then print a phased genotype, and
            // the record may not name a phase block. Decided independently,
            // they disagreed on 120157 of the 131189 UNPHASED records of the
            // HiFi HCC1395 chr1 pair. Holds by construction: want_ps is 0 on
            // every UNPHASED record, and append_paired_sample() prints "|"
            // only under want_ps.
            if (len + (int)sks.l + 2 >= buf_m) {
                buf_m = len + (int)sks.l + 1024;
                buffer = (char*)_err_realloc(buffer, buf_m);
            }
            len += snprintf(buffer + len, buf_m - len, "%s", sks.s);
            ks_free(&sks);
        } else {
            int gt1 = var.GT[0], gt2 = var.GT[1];
            int is_hom = gt1 == gt2; int gt_seperator = '|';
            // <= 0, not == 0: assign_hap.c's somatic path uses -1 for "no phase
            // set" where the germline path uses 0, and a record with PS = -1
            // printed "|" while the ":PS" field it refers to was suppressed.
            if (var.PS <= 0) {
                gt_seperator = '/';
                if (gt1 > gt2) { int tmp = gt1; gt1 = gt2; gt2 = tmp;  }
            }
            // add SVREADS to FORMAT if outputting supporting read names for SVs
            len += snprintf(buffer + len, buf_m - len, "GT:DP:AD:VAF:GQ:PL");
            if (is_hom == 0 && var.PS > 0) len += snprintf(buffer + len, buf_m - len, ":PS");
            len += snprintf(buffer + len, buf_m - len, ":UPR:UPF");
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
            // GQ:PL
            len += snprintf(buffer + len, buf_m - len, ":%d:%d,%d,%d", var.GQ, var.vq.pl[0][0], var.vq.pl[0][1], var.vq.pl[0][2]);
            // PS
            if (is_hom == 0 && var.PS > 0)
                len += snprintf(buffer + len, buf_m - len, ":%" PRId64 "", var.PS);
            {
                kstring_t uks = KS_INITIALIZE;
                append_unphased_fields(&uks, &var, 0);
                len += snprintf(buffer + len, buf_m - len, "%s", uks.s);
                ks_free(&uks);
            }
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
