#!/usr/bin/env python3
"""Synthetic tumour/normal pair from one real long-read BAM, with exact truth.

Reads are split in two by a hash of their name. HG002's hets are phased by
trio (HG003/HG004 GIAB truth), so phase is independent of the caller under
test; each read is assigned the parent whose alleles it carries. Synthetic
somatic SNVs and indels are then written into the tumour half's reads of one
haplotype only -- every read, or a hash-chosen half of them (subclonal) -- so
the injected alleles inherit the real per-read error profile and the real
haplotype structure the paired caller has to work with.

usage: make_pair.py BAM REF.fa HG002.vcf HG002.bed HG003.vcf HG003.bed HG004.vcf HG004.bed
                    REGION(chr:beg-end) OUTDIR N_SNV N_INDEL SEED
"""
import sys, random, hashlib, bisect, os
import pysam

(bam_fn, ref_fn, kid_vcf, kid_bed, dad_vcf, dad_bed, mom_vcf, mom_bed,
 region, outdir, n_snv, n_indel, seed) = sys.argv[1:14]
n_snv, n_indel = int(n_snv), int(n_indel)
rng = random.Random(int(seed))
ctg, rr = region.split(':'); rbeg, rend = (int(x.replace(',', '')) for x in rr.split('-'))
os.makedirs(outdir, exist_ok=True)
ref = pysam.FastaFile(ref_fn)
REF = ref.fetch(ctg, 0, rend + 100000).upper()  # 0-based string

def load_bed(fn):
    iv = []
    for l in open(fn):
        t = l.split()
        if t[0] == ctg and int(t[2]) > rbeg and int(t[1]) < rend: iv.append((int(t[1]), int(t[2])))
    iv.sort(); return iv

def in_bed(iv, starts, p0):  # p0 0-based
    i = bisect.bisect_right(starts, p0) - 1
    return i >= 0 and iv[i][0] <= p0 < iv[i][1]

def load_vcf(fn):
    d = {}
    for l in open(fn):
        if l[0] == '#': continue
        t = l.rstrip('\n').split('\t')
        if t[0] != ctg: continue
        gt = t[9].split(':')[0].replace('|', '/')
        d[int(t[1])] = (t[3], t[4], gt)
    return d

kid = load_vcf(kid_vcf); dad = load_vcf(dad_vcf); mom = load_vcf(mom_vcf)
kb, db, mb = load_bed(kid_bed), load_bed(dad_bed), load_bed(mom_bed)
kbs, dbs, mbs = [x[0] for x in kb], [x[0] for x in db], [x[0] for x in mb]

def parent_alt(d, iv, ivs, pos, ref_, alt):
    """alt copies in a parent at this site, None if unknown"""
    if pos in d:
        r, a, gt = d[pos]
        if r != ref_ or a != alt: return None
        return {'0/0': 0, '0/1': 1, '1/0': 1, '1/1': 2}.get(gt)
    return 0 if in_bed(iv, ivs, pos - 1) else None

# trio-phased het SNVs: pos -> (ref, alt, hap carrying alt: 1 paternal, 2 maternal)
phased = {}
for pos, (r, a, gt) in kid.items():
    if gt not in ('0/1', '1/0') or len(r) != 1 or len(a) != 1 or ',' in a: continue
    da = parent_alt(dad, db, dbs, pos, r, a); ma = parent_alt(mom, mb, mbs, pos, r, a)
    if da is None or ma is None: continue
    if da == 0 and ma >= 1: phased[pos] = (r, a, 2)
    elif ma == 0 and da >= 1: phased[pos] = (r, a, 1)
    elif da == 2 and ma == 1: phased[pos] = (r, a, 1)
    elif ma == 2 and da == 1: phased[pos] = (r, a, 2)
ppos = sorted(phased)
print(f"het SNVs {sum(1 for v in kid.values() if v[2] in ('0/1','1/0'))}, trio-phased {len(ppos)}", file=sys.stderr)

def h01(name, salt):  # deterministic [0,1)
    return int(hashlib.md5((salt + name).encode()).hexdigest()[:8], 16) / 2**32

def read_hap(rd):
    if rd.is_unmapped or rd.is_secondary or rd.is_supplementary or rd.mapping_quality < 20: return 0
    i = bisect.bisect_left(ppos, rd.reference_start + 1)
    want = set()
    while i < len(ppos) and ppos[i] <= rd.reference_end: want.add(ppos[i] - 1); i += 1
    if not want: return 0
    sup = [0, 0, 0]
    seq = rd.query_sequence
    for qp, rp in rd.get_aligned_pairs(matches_only=True):
        if rp in want:
            r, a, h = phased[rp + 1]
            b = seq[qp]
            if b == a: sup[h] += 1
            elif b == r: sup[3 - h] += 1
    if sup[1] >= 2 and sup[1] >= 4 * sup[2]: return 1
    if sup[2] >= 2 and sup[2] >= 4 * sup[1]: return 2
    return 0

# somatic sites: inside HG002's confident regions, >= 50 bp from any HG002
# variant, >= 5 kb apart, no N in the window
kvpos = sorted(kid)
def near_germline(p, d=50):
    i = bisect.bisect_left(kvpos, p - d)
    return i < len(kvpos) and kvpos[i] <= p + d
sites = []; taken = []
tries = 0
types = ['snv'] * n_snv + ['del'] * (n_indel // 2) + ['ins'] * (n_indel - n_indel // 2)
rng.shuffle(types)
for ty in types:
    while True:
        tries += 1
        if tries > 10**6: sys.exit("could not place sites")
        p = rng.randint(rbeg + 20000, rend - 20000)  # 1-based
        if not in_bed(kb, kbs, p - 1) or near_germline(p): continue
        if any(abs(p - q) < 5000 for q in taken): continue
        win = REF[p - 60:p + 60]
        if 'N' in win: continue
        break
    taken.append(p)
    hap = rng.choice((1, 2)); frac = rng.choice((1.0, 0.5))
    if ty == 'snv':
        alt = rng.choice([b for b in 'ACGT' if b != REF[p - 1]]); sites.append((p, ty, alt, hap, frac))
    elif ty == 'del':
        L = rng.choice((1, 1, 1, 2, 2, 3, 4, 5, 8, 12)); sites.append((p, ty, L, hap, frac))
    else:
        L = rng.choice((1, 1, 1, 2, 2, 3, 4, 5, 8, 12))
        sites.append((p, ty, ''.join(rng.choice('ACGT') for _ in range(L)), hap, frac))
sites.sort()
spos = [s[0] for s in sites]

def expand(rd):
    """per-op list [(op, len)] -- pysam codes: 0 M, 1 I, 2 D, 4 S, 5 H, 7 =, 8 X"""
    return [list(c) for c in rd.cigartuples]

def find_qpos(cig, rstart, target0, need_len=1):
    """index of the M/=/X op that aligns ref [target0, target0+need_len) and query offset of target0"""
    rp, qp = rstart, 0
    for k, (op, n) in enumerate(cig):
        if op in (0, 7, 8):
            if rp <= target0 and target0 + need_len <= rp + n: return k, qp + (target0 - rp), target0 - rp
            rp += n; qp += n
        elif op in (2, 3):
            if rp <= target0 < rp + n: return None
            rp += n
        elif op in (1, 4): qp += n
    return None

def inject(rd, site):
    p, ty, x, hap, frac = site
    t0 = p - 1
    seq = rd.query_sequence; qual = rd.query_qualities
    cig = expand(rd)
    if ty == 'snv':
        f = find_qpos(cig, rd.reference_start, t0)
        if f is None: return False
        k, q, off = f
        seq = seq[:q] + x + seq[q + 1:]
    elif ty == 'del':
        f = find_qpos(cig, rd.reference_start, t0, x)
        if f is None: return False
        k, q, off = f
        op, n = cig[k]
        seq = seq[:q] + seq[q + x:]; qual = qual[:q] + qual[q + x:]
        new = [[op, off]] if off else []
        new.append([2, x])
        if n - off - x: new.append([op, n - off - x])
        cig = cig[:k] + new + cig[k + 1:]
    else:  # insertion before ref base p
        f = find_qpos(cig, rd.reference_start, t0)
        if f is None: return False
        k, q, off = f
        op, n = cig[k]
        if off == 0: return False  # keep it simple: insertion strictly inside an aligned run
        qq = qual[q - 1] if qual is not None else None
        seq = seq[:q] + x + seq[q:]
        if qual is not None: qual = qual[:q] + pysam.qualitystring_to_array(chr(qq + 33) * len(x)) + qual[q:]
        cig = cig[:k] + [[op, off], [1, len(x)], [op, n - off]] + cig[k + 1:]
    rd.query_sequence = seq
    rd.cigartuples = [tuple(c) for c in cig]
    if qual is not None: rd.query_qualities = qual
    return True

inb = pysam.AlignmentFile(bam_fn)
hdr = inb.header
outN = pysam.AlignmentFile(f"{outdir}/normal.unsorted.bam", 'wb', header=hdr)
outT = pysam.AlignmentFile(f"{outdir}/tumor.unsorted.bam", 'wb', header=hdr)
n_inj = [0] * len(sites); n_cov = [0] * len(sites)
nr = [0, 0]
for rd in inb.fetch(ctg, rbeg, rend):
    tumour = h01(rd.query_name, 'split') < 0.5
    if not tumour:
        outN.write(rd); nr[0] += 1; continue
    nr[1] += 1
    if not rd.is_unmapped and not rd.is_secondary:
        i = bisect.bisect_left(spos, rd.reference_start + 30)
        hits = []
        while i < len(spos) and spos[i] <= rd.reference_end - 30: hits.append(i); i += 1
        if hits:
            h = read_hap(rd)
            # inject right to left so earlier coordinates stay valid
            for i in reversed(hits):
                s = sites[i]
                n_cov[i] += 1
                if h != s[3] or h01(rd.query_name, f"cl{s[0]}") >= s[4]: continue
                for tag in ('MD', 'NM', 'cs', 'de', 'ms', 'AS'):
                    if rd.has_tag(tag): rd.set_tag(tag, None)
                if inject(rd, s): n_inj[i] += 1
    outT.write(rd)
outN.close(); outT.close()
for nm in ('normal', 'tumor'):
    pysam.sort('-o', f"{outdir}/{nm}.bam", f"{outdir}/{nm}.unsorted.bam"); pysam.index(f"{outdir}/{nm}.bam")
    os.remove(f"{outdir}/{nm}.unsorted.bam")
print(f"reads normal {nr[0]} tumour {nr[1]}", file=sys.stderr)

with open(f"{outdir}/truth.raw.vcf", 'w') as o:
    o.write("##fileformat=VCFv4.2\n")
    o.write(f"##contig=<ID={ctg},length={ref.get_reference_length(ctg)}>\n")
    o.write('##INFO=<ID=HAP,Number=1,Type=Integer,Description="1 paternal, 2 maternal">\n')
    o.write('##INFO=<ID=CCF,Number=1,Type=Float,Description="share of the haplotype\'s tumour reads injected">\n')
    o.write('##INFO=<ID=NINJ,Number=1,Type=Integer,Description="tumour reads carrying it">\n')
    o.write('##INFO=<ID=TDP,Number=1,Type=Integer,Description="tumour reads spanning it">\n')
    o.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
    for i, (p, ty, x, hap, frac) in enumerate(sites):
        if n_inj[i] == 0: continue
        if ty == 'snv': pos, r, a = p, REF[p - 1], x
        elif ty == 'del': pos, r, a = p - 1, REF[p - 2:p - 1 + x], REF[p - 2]
        else: pos, r, a = p - 1, REF[p - 2], REF[p - 2] + x
        o.write(f"{ctg}\t{pos}\t.\t{r}\t{a}\t.\tPASS\tHAP={hap};CCF={frac};NINJ={n_inj[i]};TDP={n_cov[i]}\n")
inj = [n for n in n_inj if n]
print(f"sites {len(sites)}, injected {len(inj)}, reads/site median {sorted(inj)[len(inj)//2] if inj else 0}", file=sys.stderr)
