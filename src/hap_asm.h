#ifndef LONGCALLD_HAP_ASM_H
#define LONGCALLD_HAP_ASM_H

// Read-backed haplotype assembly and haplotype-informed genotyping, after
// Longshot (Edge & Bansal 2019) and the HapCUT2 objective it uses.
//
// Each read is a fragment of allele calls, every call weighed by its own error
// model. Haplotypes are found by local search on the fragment likelihood
// sum_f log(P(f|H) + P(f|~H)), and every variant is then re-genotyped on the
// reads' haplotype likelihoods; the two steps alternate until the likelihood
// settles. Phase blocks are the connected components of the reads: two hets a
// read informatively covers are in the same block, always. A het the reads do
// not agree on is pruned from the phasing (HapInconsistent) instead of being
// allowed to switch or split the block.

#include <stdint.h>
#include "htslib/hts.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    int v;          // variant index
    uint8_t a;      // allele the read shows: 0 ref, 1 alt
    uint8_t ok;     // trusted enough to link variants into blocks
    uint8_t tag_only; // from a read calling skipped: used for its haplotype tag and nothing else
    float lm, lx;   // ln P(call | haplotype allele == a), ln P(call | haplotype allele != a)
} ha_call_t;

typedef struct {
    int off, n;     // calls call[off .. off+n-1], sorted by v
    uint8_t smp;    // sample the read came from
} ha_frag_t;

typedef struct {
    hts_pos_t pos;          // record position; a block is named after its first variant
    uint8_t phasable;       // may join the phasing graph
    uint8_t geno;           // genotype here: 2 full, 1 orientation only (zygosity fixed), 0 not at all
    uint8_t is_indel;
    int8_t g[2];            // ordered genotype: allele on haplotype 1 and 2; in: seed, out: result
    int8_t ug;              // alt alleles of the genotype before phasing (0..2)
    uint8_t pruned;         // removed from the phasing graph: its reads contradict each other
    uint8_t inconsistent;   // HapInconsistent: pruned, or alt on no haplotype once phased
    uint8_t linked;         // in the final phasing graph, tied by reads to another het of its block
    uint8_t discordant;     // HapDiscordant: too many reads of a haplotype show the other allele
    uint16_t hc[2][2];      // per haplotype: reads showing that haplotype's allele, reads in all
    hts_pos_t ps;           // phase set, 0 if none
    float pq, sq, gq;       // phred: phase quality, switch quality to the previous het, genotype quality; -1 unset
} ha_var_t;

typedef struct {
    double max_miscall;     // a call links variants only with error below this (Longshot --min_allele_qual)
    double min_pq;          // a het whose single flip costs less (phred) is pruned, given >=2 linking reads
    double hap_assign_p;    // posterior a read needs for HP
    double prior[2][3];     // [snv, indel][0/0, het (both orders), 1/1]
    int geno_smp;           // genotype from this sample's reads only; -1 all
    int max_rounds;         // assembly / genotyping rounds
    int solo_needs_phasable; // an unlinked het may be a block of its own only if phasable (paired anchoring)
    double min_hap_conc;    // share of a haplotype's reads that must show its allele; 0 off
    int min_hap_conc_dp;    // reads a haplotype needs before it is held to that
    double conc_assign_p;   // posterior a read needs, on its other sites, to count toward a haplotype
    double conc_max_miscall; // a call counts toward concordance only with error below this
    double conc_err[2];     // [snv, indel] read error the off-allele reads are tested against
    double conc_pval;       // ...and the binomial tail they must fall below to be called discordant
} ha_opt_t;

typedef struct {
    int n_var; ha_var_t *var;
    int n_frag; ha_frag_t *frag;
    int n_call; ha_call_t *call;
    // out, per fragment: phase set and haplotype (0 untagged, 1, 2)
    hts_pos_t *frag_ps; int8_t *frag_hap;
} ha_t;

void ha_opt_init(ha_opt_t *o);
void ha_run(const ha_opt_t *o, ha_t *h);
void ha_free(ha_t *h);

struct call_var_opt_t;
struct bam_chunk_t;
// Assemble the haplotypes of one stitch batch, chromosome by chromosome, and
// write genotypes, phase sets and read haplotypes back into the chunks.
void hap_asm_chunks(const struct call_var_opt_t *opt, struct bam_chunk_t *chunks, int n_chunks);
// End of step 0, while digars exist: alleles of skipped reads at germline SNV
// candidates, so the reads can still be tagged.
void hap_asm_collect_skipped(const struct call_var_opt_t *opt, struct bam_chunk_t *chunk);

#ifdef __cplusplus
}
#endif

#endif
