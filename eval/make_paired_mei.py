#!/usr/bin/env python3
"""Build a synthetic tumour/normal pair carrying mobile element insertions.

usage: make_paired_mei.py ref.fa te.fa outdir

Writes outdir/{ref.fa,tumor.sam,normal.sam,events.tsv}. The reference is a 60 kb
slice of chr11 renamed "t"; two haplotypes differ from it by a het SNP every
kb, which is what phasing needs. Reads are error-free 8 kb HiFi-like reads
with exact =/X/I CIGARs, so no aligner is involved and the alt read counts are
exactly the ones requested. Each event is "insert SEQ before reference
position POS on haplotype 1" in the reads listed for it:

  A  Alu (+), TSD 12, polyA 30   2 tumour reads, no normal read   -> somatic MEI
  B  Alu (-), TSD 12, polyT 30   1 tumour read,  no normal read   -> somatic MEI
  C  Alu (+), TSD 12, polyA 30   2 tumour reads, 1 normal read    -> not somatic
  D  300 bp random, no TSD/polyA 2 tumour reads, no normal read   -> not somatic
"""
import random, sys, os

ref_fn, te_fn, out = sys.argv[1:4]
os.makedirs(out, exist_ok=True)
rng = random.Random(11)
COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]


def read_fa(fn):
    seqs, name = {}, None
    for l in open(fn):
        l = l.strip()
        if l.startswith(">"): name = l[1:].split()[0]; seqs[name] = []
        elif name: seqs[name].append(l.upper())
    return {k: "".join(v) for k, v in seqs.items()}


ref = read_fa(ref_fn)["chr11"][400000:460000]
assert "N" not in ref
alu = read_fa(te_fn)["AluY"]
L, RLEN, DEPTH = len(ref), 8000, 15

# het SNPs, alternating haplotypes, kept 1.5 kb clear of the events
EV = {"A": 15000, "B": 30000, "C": 45000, "D": 52000}
snp = {}
for k, x in enumerate(range(700, L - 700, 1000)):
    if min(abs(x - p) for p in EV.values()) < 1500: continue
    snp[x] = (1 + k % 2, rng.choice([c for c in "ACGT" if c != ref[x]]))


def ins_seq(e):
    p = EV[e]
    tsd = ref[p:p + 12]
    if e in "AC": return tsd + alu + "A" * 30
    if e == "B": return tsd + rc(alu + "A" * 30)
    return "".join(rng.choice("ACGT") for _ in range(300))


INS = {e: ins_seq(e) for e in EV}


def hap_bases(h, events):
    """Haplotype h as (base, ref_pos or None) with the given insertions."""
    out = []
    for x in range(L):
        for e in events:
            if EV[e] == x: out += [(c, None) for c in INS[e]]
        b = snp[x][1] if x in snp and snp[x][0] == h else ref[x]
        out.append((b, x))
    return out


def sam_line(name, bases, beg, end):
    seg = bases[beg:end]
    while seg and seg[0][1] is None: seg = seg[1:]   # a read must start and end aligned
    while seg and seg[-1][1] is None: seg = seg[:-1]
    ops = []
    for b, x in seg:
        op = "I" if x is None else ("=" if b == ref[x] else "X")
        if ops and ops[-1][0] == op: ops[-1][1] += 1
        else: ops.append([op, 1])
    cig = "".join(f"{n}{o}" for o, n in ops)
    seq = "".join(b for b, _ in seg)
    flag = rng.choice((0, 16))
    return f"{name}\t{flag}\tt\t{seg[0][1] + 1}\t60\t{cig}\t*\t0\t0\t{seq}\t{'?' * len(seq)}\n"


def tiled(h, prefix, f):
    bases = hap_bases(h, [])
    n = 0
    while n * RLEN < DEPTH * L:
        st = rng.randint(0, L - RLEN)
        f.write(sam_line(f"{prefix}{h}_{n}", bases, st, st + RLEN)); n += 1


def carriers(e, k, prefix, f):
    bases = hap_bases(1, [e])
    pos = next(i for i, (_, x) in enumerate(bases) if x is None)
    for j in range(k):
        st = pos - rng.randint(2000, RLEN - len(INS[e]) - 2000)
        f.write(sam_line(f"{prefix}{e}_{j}", bases, st, st + RLEN))


hdr = f"@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:t\tLN:{L}\n"
with open(f"{out}/ref.fa", "w") as f:
    f.write(">t\n" + "\n".join(ref[i:i + 80] for i in range(0, L, 80)) + "\n")
with open(f"{out}/tumor.sam", "w") as f:
    f.write(hdr + "@RG\tID:T\tSM:tumor\n")
    for h in (1, 2): tiled(h, "t", f)
    for e, k in (("A", 2), ("B", 1), ("C", 2), ("D", 2)): carriers(e, k, "tsom", f)
with open(f"{out}/normal.sam", "w") as f:
    f.write(hdr + "@RG\tID:N\tSM:normal\n")
    for h in (1, 2): tiled(h, "n", f)
    carriers("C", 1, "nsom", f)
with open(f"{out}/events.tsv", "w") as f:
    for e, p in EV.items(): f.write(f"{e}\t{p}\t{len(INS[e])}\n")
