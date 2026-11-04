#!/usr/bin/env python3
"""Phasing QC for a longcallD VCF and its HP/PS-tagged alignment output.

Reports phase-block statistics and counts violations of the invariants a
read-backed phasing must satisfy:

  I1  a read carries informative alleles at phased hets of two different PS
      (the reads link the blocks, so they should have been one block)
  I2  a PASS germline het in the phasing sample is printed unphased
  I3  a PASS somatic SNV has no PS although tagged reads carry its alt
  I4  an untagged read has >=2 informative alleles at phased hets of one PS,
      all on the same haplotype (it should have been tagged)
  I5  a tagged read disagrees with its own HP at >20% of >=3 sites of its PS,
      or carries a PS that none of its informative sites belongs to

It also recounts, for every pair of consecutive phased hets, the reads that
link them cis or trans relative to the VCF, so a block held together by a weak
or contradicted link and a block boundary that the reads clearly bridge are
both visible.

With --truth (a phased VCF, e.g. GIAB HG002 trio-phased) it counts switch and
flip errors over hets phased in both, and genotype concordance of the phasing
sample, optionally restricted to --bed.

Read-level checks use biallelic SNVs only. Requires pysam.

  phase_qc.py --vcf out.vcf.gz --bam tumor.bam --bam normal.bam [--region chr1]
"""
import argparse
import bisect
import collections
import json
import re
import sys

import pysam

LINK_MIN = 2         # reads for a link to count as supported
LINK_CONFLICT = 0.2  # minority fraction above which a link is contradicted


def n50(lengths):
    if not lengths:
        return 0
    s = sorted(lengths, reverse=True)
    half, acc = sum(s) / 2.0, 0
    for x in s:
        acc += x
        if acc >= half:
            return x
    return s[-1]


def load_bed(path):
    if not path:
        return None
    iv = collections.defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(('#', 'track', 'browser')):
                continue
            f = line.split()
            iv[f[0]].append((int(f[1]), int(f[2])))
    for c in iv:
        iv[c].sort()
    return iv


def in_bed(bed, chrom, pos0):
    if bed is None:
        return True
    lst = bed.get(chrom)
    if not lst:
        return False
    i = bisect.bisect_right(lst, (pos0, float('inf'))) - 1
    return i >= 0 and lst[i][0] <= pos0 < lst[i][1]


def iter_region(vf, region):
    """Fetch a region, scanning the file when it has no index."""
    if not region:
        return vf.fetch()
    try:
        return vf.fetch(region=region)
    except ValueError:
        m = re.match(r'^([^:]+)(?::([\d,]+)-([\d,]+))?$', region)
        chrom = m.group(1)
        beg = int(m.group(2).replace(',', '')) if m.group(2) else 1
        end = int(m.group(3).replace(',', '')) if m.group(3) else float('inf')
        return (r for r in vf.fetch() if r.chrom == chrom and beg <= r.pos <= end)


def is_pass(rec):
    f = list(rec.filter.keys())
    return not f or f == ['PASS']


def gt_of(sample):
    gt = sample.get('GT')
    if gt is None or None in gt:
        return None, False
    return tuple(gt), bool(sample.phased)


def load_calls(vcf, region, sidx, tidx, min_qual=0):
    """Het records of the phasing sample, somatic records, and all PASS germline genotypes."""
    hets, soms, gts = [], [], {}
    vf = pysam.VariantFile(vcf)
    names = list(vf.header.samples)
    for rec in iter_region(vf, region):
        if min_qual > 0 and (rec.qual is None or rec.qual < min_qual):
            continue
        s = rec.samples[names[sidx]]
        gt, phased = gt_of(s)
        som = 'SOMATIC' in rec.info
        alts = rec.alts or ()
        snv = len(rec.ref) == 1 and len(alts) == 1 and len(alts[0]) == 1 and alts[0] in 'ACGT'
        if som:
            t = rec.samples[names[tidx]] if tidx is not None else s
            tgt, tph = gt_of(t)
            tps = t.get('PS') if 'PS' in rec.format else None
            soms.append(dict(chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alt=alts[0] if alts else '',
                             snv=snv, passed=is_pass(rec), ps=tps if tph else None))
            continue
        if is_pass(rec) and gt is not None and any(gt):
            gts[(rec.chrom, rec.pos, rec.ref, alts)] = gt
        if gt is None or len(gt) != 2 or gt[0] == gt[1]:
            continue
        ps = s.get('PS') if 'PS' in rec.format else None
        ok = phased and ps is not None and ps > 0
        hets.append(dict(chrom=rec.chrom, pos=rec.pos, ref=rec.ref, alts=alts, gt=gt,
                         phased=ok, ps=ps if ok else None, passed=is_pass(rec), snv=snv))
    return names, hets, soms, gts


def read_bases_at(read, positions0):
    """Map sorted 0-based ref positions to the read's base there (absent if deleted)."""
    out = {}
    seq = read.query_sequence
    if not positions0 or seq is None:
        return out
    rpos, qpos, k, n = read.reference_start, 0, 0, len(positions0)
    for op, ln in read.cigartuples:
        if k >= n:
            break
        if op in (0, 7, 8):  # M = X
            end = rpos + ln
            while k < n and positions0[k] < end:
                p = positions0[k]
                if p >= rpos:
                    out[p] = seq[qpos + p - rpos]
                k += 1
            rpos, qpos = end, qpos + ln
        elif op in (2, 3):  # D N
            end = rpos + ln
            while k < n and positions0[k] < end:
                k += 1
            rpos = end
        elif op in (1, 4):  # I S
            qpos += ln
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--vcf', required=True)
    ap.add_argument('--bam', action='append', default=[], help='HP/PS-tagged BAM/CRAM; repeatable')
    ap.add_argument('--ref', help='reference for CRAM input')
    ap.add_argument('--region')
    ap.add_argument('--sample', type=int, default=0, help='VCF column used for phasing (normal; default 0)')
    ap.add_argument('--tumor-sample', type=int, default=None, help='VCF column of the tumour (default 1 if paired)')
    ap.add_argument('--truth', help='phased truth VCF')
    ap.add_argument('--bed', help='restrict the truth comparison to these regions')
    ap.add_argument('--min-mapq', type=int, default=1)
    ap.add_argument('--min-qual', type=float, default=0, help='ignore records below this QUAL (default: all)')
    ap.add_argument('--json', help='also write the summary as JSON here')
    ap.add_argument('--list', type=int, default=0, help='print up to N examples of each violation')
    a = ap.parse_args()

    tidx = a.tumor_sample
    if tidx is None and len(pysam.VariantFile(a.vcf).header.samples) > 1:
        tidx = 1
    names, hets, soms, gts = load_calls(a.vcf, a.region, a.sample, tidx, a.min_qual)
    R = collections.OrderedDict()
    ex = collections.defaultdict(list)

    def example(key, s):
        if len(ex[key]) < a.list:
            ex[key].append(s)

    # ---------------------------------------------------------------- blocks
    pass_hets = [h for h in hets if h['passed']]
    phased = [h for h in pass_hets if h['phased']]
    blocks = collections.OrderedDict()
    for h in phased:
        blocks.setdefault((h['chrom'], h['ps']), []).append(h)
    spans = [b[-1]['pos'] - b[0]['pos'] + 1 for b in blocks.values()]
    nvar = [len(b) for b in blocks.values()]
    R['pass_hets'] = len(pass_hets)
    R['phased_hets'] = len(phased)
    R['I2_unphased_pass_hets'] = len(pass_hets) - len(phased)
    for h in pass_hets:
        if not h['phased']:
            example('I2', '%s:%d' % (h['chrom'], h['pos']))
    R['blocks'] = len(blocks)
    R['block_N50_bp'] = n50(spans)
    R['block_N50_vars'] = n50(nvar)
    R['block_max_bp'] = max(spans) if spans else 0
    R['singleton_blocks'] = sum(1 for n in nvar if n == 1)
    R['blocks_lt_1kb_multi'] = sum(1 for s, n in zip(spans, nvar) if n > 1 and s < 1000)
    R['blocks_lt_10kb'] = sum(1 for s in spans if s < 10000)
    ivs = collections.defaultdict(list)
    for (c, ps), b in blocks.items():
        ivs[c].append((b[0]['pos'], b[-1]['pos']))
    nest = 0
    for c in ivs:
        mx = -1
        for s, e in sorted(ivs[c]):
            if s <= mx:
                nest += 1
            mx = max(mx, e)
    R['interleaved_blocks'] = nest

    # ----------------------------------------------------------- read checks
    by_chrom = collections.defaultdict(list)
    for h in pass_hets:
        if h['snv']:
            by_chrom[h['chrom']].append(h)
    som_snv = collections.defaultdict(dict)
    for s in soms:
        if s['snv'] and s['passed']:
            som_snv[s['chrom']][s['pos'] - 1] = s
    link = collections.defaultdict(lambda: [0, 0])  # (chrom, i, j) -> [cis, trans]
    som_tagged_alt = collections.Counter()
    v = collections.Counter()
    for bam in a.bam:
        af = pysam.AlignmentFile(bam, reference_filename=a.ref)
        for chrom in sorted(set(by_chrom) | set(som_snv)):
            hl = by_chrom[chrom]
            pos0 = [h['pos'] - 1 for h in hl]
            spos = sorted(som_snv[chrom])
            fetch = af.fetch(region=a.region) if a.region else af.fetch(chrom)
            for r in fetch:
                if r.is_unmapped or r.is_secondary or r.is_supplementary or r.mapping_quality < a.min_mapq:
                    continue
                v['reads'] += 1
                lo = bisect.bisect_left(pos0, r.reference_start)
                hi = bisect.bisect_left(pos0, r.reference_end)
                slo = bisect.bisect_left(spos, r.reference_start)
                shi = bisect.bisect_left(spos, r.reference_end)
                hp = r.get_tag('HP') if r.has_tag('HP') else 0
                rps = r.get_tag('PS') if r.has_tag('PS') else None
                if hp:
                    v['tagged'] += 1
                bases = read_bases_at(r, sorted(pos0[lo:hi] + spos[slo:shi]))
                if hp:
                    for p in spos[slo:shi]:
                        if bases.get(p) == som_snv[chrom][p]['alt']:
                            som_tagged_alt[(chrom, p)] += 1
                ph = []  # (het index, ps, haplotype 1/2 relative to the VCF)
                for i in range(lo, hi):
                    h = hl[i]
                    b = bases.get(pos0[i])
                    alle = [h['ref']] + list(h['alts'])
                    if b is None or b not in alle:
                        continue
                    ai = alle.index(b)
                    if ai not in h['gt'] or not h['phased']:
                        continue
                    ph.append((i, h['ps'], 1 if h['gt'][0] == ai else 2))
                if not ph:
                    continue
                v['reads_informative'] += 1
                for x, y in zip(ph, ph[1:]):
                    link[(chrom, x[0], y[0])][0 if x[2] == y[2] else 1] += 1
                pss = set(x[1] for x in ph)
                if len(pss) > 1:
                    v['I1_reads_spanning_blocks'] += 1
                    example('I1', '%s %s:%d-%d PS=%s' % (r.query_name, chrom, r.reference_start, r.reference_end,
                                                         ','.join(map(str, sorted(pss)))))
                if not hp:
                    per = collections.defaultdict(list)
                    for x in ph:
                        per[x[1]].append(x[2])
                    if any(len(l) >= 2 for l in per.values()):
                        v['untagged_reads_ge2_phased'] += 1
                    if any(len(l) >= 2 and len(set(l)) == 1 for l in per.values()):
                        v['I4_untagged_consistent_reads'] += 1
                        example('I4', '%s %s:%d' % (r.query_name, chrom, r.reference_start))
                    continue
                mine = [x for x in ph if x[1] == rps]
                if not mine:
                    v['I5_tagged_ps_not_covered'] += 1
                    example('I5ps', '%s %s:%d PS=%s sites_PS=%s' % (r.query_name, chrom, r.reference_start, rps,
                                                                    ','.join(map(str, sorted(pss)))))
                    continue
                if len(mine) >= 3:
                    dis = sum(1 for x in mine if x[2] != hp)
                    v['tagged_reads_ge3'] += 1
                    if dis > LINK_CONFLICT * len(mine):
                        v['I5_tagged_discordant'] += 1
                        example('I5', '%s %s:%d HP=%d %d/%d' % (r.query_name, chrom, r.reference_start, hp, dis, len(mine)))
    for k in ('reads', 'tagged', 'reads_informative', 'I1_reads_spanning_blocks', 'untagged_reads_ge2_phased',
              'I4_untagged_consistent_reads', 'tagged_reads_ge3', 'I5_tagged_discordant', 'I5_tagged_ps_not_covered'):
        if a.bam:
            R[k] = v[k]
    # link support between consecutive phased SNV hets
    if a.bam:
        cin, cut = collections.Counter(), collections.Counter()
        for chrom, hl in by_chrom.items():
            idx = [i for i, h in enumerate(hl) if h['phased']]
            for i, j in zip(idx, idx[1:]):
                cis, trans = link.get((chrom, i, j), (0, 0))
                tot, mn = cis + trans, min(cis, trans)
                tag = '%s:%d-%d cis=%d trans=%d' % (chrom, hl[i]['pos'], hl[j]['pos'], cis, trans)
                if hl[i]['ps'] == hl[j]['ps']:
                    cin['pairs'] += 1
                    if tot < LINK_MIN:
                        cin['unsupported'] += 1
                    elif mn > LINK_CONFLICT * tot:
                        cin['contradicted'] += 1
                        example('link_contradicted', tag)
                    elif trans > cis:
                        cin['phase_disagrees_with_reads'] += 1
                        example('link_flipped', tag)
                else:
                    cut['boundaries'] += 1
                    if tot >= LINK_MIN and mn <= LINK_CONFLICT * tot:
                        cut['decisively_bridged'] += 1
                        example('cut_bridged', tag)
        for k in ('pairs', 'unsupported', 'contradicted', 'phase_disagrees_with_reads'):
            R['inblock_link_' + k] = cin[k]
        for k in ('boundaries', 'decisively_bridged'):
            R['cut_' + k] = cut[k]

    # --------------------------------------------------------------- somatic
    if soms:
        sp = [s for s in soms if s['passed']]
        R['somatic_pass'] = len(sp)
        R['somatic_phased'] = sum(1 for s in sp if s['ps'])
        if a.bam:
            reach = [s for s in sp if s['snv'] and som_tagged_alt[(s['chrom'], s['pos'] - 1)] > 0]
            R['somatic_snv_reachable'] = len(reach)
            R['I3_somatic_reachable_unphased'] = sum(1 for s in reach if not s['ps'])
            for s in reach:
                if not s['ps']:
                    example('I3', '%s:%d' % (s['chrom'], s['pos']))

    # ----------------------------------------------------------------- truth
    if a.truth:
        bed = load_bed(a.bed)
        tf = pysam.VariantFile(a.truth)
        tsn = list(tf.header.samples)[0]
        truth = {}
        for rec in iter_region(tf, a.region):
            gt, ph = gt_of(rec.samples[tsn])
            if gt is None or not any(gt) or not in_bed(bed, rec.chrom, rec.pos - 1):
                continue
            truth[(rec.chrom, rec.pos, rec.ref, rec.alts)] = (gt, ph)
        q = {k: g for k, g in gts.items() if in_bed(bed, k[0], k[1] - 1)}
        tp = [k for k in q if k in truth]
        R['truth_variants'] = len(truth)
        R['call_variants'] = len(q)
        R['allele_TP'] = len(tp)
        R['allele_precision'] = round(len(tp) / len(q), 5) if q else 0
        R['allele_recall'] = round(len(tp) / len(truth), 5) if truth else 0
        R['genotype_concordance_of_TP'] = round(sum(1 for k in tp if sorted(q[k]) == sorted(truth[k][0])) / len(tp), 5) if tp else 0
        R['het_called_hom_in_truth'] = sum(1 for k in tp if q[k][0] != q[k][1] and truth[k][0][0] == truth[k][0][1])
        R['hom_called_het_in_truth'] = sum(1 for k in tp if q[k][0] == q[k][1] and truth[k][0][0] != truth[k][0][1])
        sw = fl = pairs = cover = 0
        for (c, ps), b in blocks.items():
            seq, at = [], []
            for h in b:
                t = truth.get((h['chrom'], h['pos'], h['ref'], h['alts']))
                if not t or not t[1] or t[0][0] == t[0][1] or sorted(t[0]) != sorted(h['gt']):
                    continue
                seq.append(0 if tuple(h['gt']) == tuple(t[0]) else 1)
                at.append(h['pos'])
            cover += len(seq)
            d = [x != y for x, y in zip(seq, seq[1:])]
            for i, x in enumerate(d):
                if x:
                    example('switch', '%s:%d-%d PS=%s' % (c, at[i], at[i + 1], ps))
            pairs += len(d)
            sw += sum(d)
            i = 0
            while i < len(d) - 1:  # two switches in a row around one variant = one flip
                if d[i] and d[i + 1]:
                    fl += 1
                    i += 2
                else:
                    i += 1
        R['truth_phased_hets_compared'] = cover
        R['switch_pairs'] = pairs
        R['switch_errors'] = sw
        R['switch_rate'] = round(sw / pairs, 5) if pairs else 0
        R['flip_errors'] = fl
        R['switch_minus_flip_rate'] = round((sw - 2 * fl) / pairs, 5) if pairs else 0

    w = max(len(k) for k in R)
    for k, x in R.items():
        print('%-*s %s' % (w, k, x))
    for k, l in ex.items():
        for s in l:
            print('example\t%s\t%s' % (k, s))
    if a.json:
        with open(a.json, 'w') as fh:
            json.dump(R, fh, indent=1)


if __name__ == '__main__':
    sys.exit(main())
