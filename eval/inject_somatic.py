#!/usr/bin/env python3
"""Inject a synthetic somatic SNV or deletion into one haplotype of a SAM stream.

The bundled test files are a single sample, so splitting them in half gives a
tumour/normal pair with no somatic variants at all -- a negative control. To get
a positive control we have to manufacture one.

A real clonal somatic SNV sits on one chromosome copy of the tumour and nowhere
in the normal. To reproduce that, reads are partitioned by their allele at a
nearby heterozygous site (which is what distinguishes the two haplotypes), and
the new base is written only into the reads on one side.

CIGARs here are explicit =/X, so a substitution is not just a SEQ edit: the '='
run covering the site has to be split into '=', '1X', '='. Changing the base
without that leaves a record whose CIGAR contradicts its sequence.

Usage: inject_somatic.py HET_POS HET_ALT VAR_POS NEW_BASE [FRACTION] < in.sam > out.sam
NEW_BASE is a single base for a substitution, or "del:N" to delete N reference
bases starting at VAR_POS -- which is what an artifact check keyed on repeat
context needs, since a substitution can never exercise it.
All positions are 1-based reference coordinates. FRACTION (default 1.0) injects
into only that share of the haplotype's reads, simulating a subclonal variant;
the choice is a hash of the read name, so it is deterministic.
"""
import hashlib, re, sys

CIG = re.compile(r'(\d+)([MIDNSHP=X])')
CONSUMES_REF = set("MDN=X")
CONSUMES_QRY = set("MIS=X")


def query_offset(pos, cigar_ops, target):
    """Query index of reference position `target`, or None if not aligned there."""
    rp, qp = pos, 0
    for n, op in cigar_ops:
        if op in CONSUMES_REF and op in CONSUMES_QRY:      # =, X, M
            if rp <= target < rp + n:
                return qp + (target - rp)
            rp += n; qp += n
        elif op in CONSUMES_REF:                            # D, N
            if rp <= target < rp + n:
                return None                                 # deleted in this read
            rp += n
        elif op in CONSUMES_QRY:                            # I, S
            qp += n
    return None


def split_match_at(cigar_ops, pos, target):
    """Rewrite the =/M op covering `target` as (left=, 1X, right=)."""
    out, rp = [], pos
    for n, op in cigar_ops:
        if op in ("=", "M") and rp <= target < rp + n:
            left, right = target - rp, n - (target - rp) - 1
            if left:  out.append((left, op))
            out.append((1, "X"))
            if right: out.append((right, op))
            rp += n
        else:
            out.append((n, op))
            if op in CONSUMES_REF: rp += n
    merged = []
    for n, op in out:                       # keep the CIGAR canonical
        if merged and merged[-1][1] == op: merged[-1] = (merged[-1][0] + n, op)
        else: merged.append((n, op))
    return merged


def delete_at(cigar_ops, pos, target, dlen):
    """Rewrite the =/M op covering `target` as (left=, dlen D, right=).

    Returns None when the deletion would not fit inside a single match run,
    which is the honest answer: splicing it across an existing indel would
    manufacture a read the aligner would never have produced.
    """
    out, rp, done = [], pos, False
    for n, op in cigar_ops:
        if not done and op in ("=", "M") and rp <= target < rp + n:
            left = target - rp
            right = n - left - dlen
            if right < 0: return None                    # runs off the end of the run
            if left:  out.append((left, op))
            out.append((dlen, "D"))
            if right: out.append((right, op))
            done = True
            rp += n
        else:
            out.append((n, op))
            if op in CONSUMES_REF: rp += n
    if not done: return None
    merged = []
    for n, op in out:
        if merged and merged[-1][1] == op: merged[-1] = (merged[-1][0] + n, op)
        else: merged.append((n, op))
    return merged


def keep(name, fraction):
    if fraction >= 1.0: return True
    h = int(hashlib.sha1(name.encode()).hexdigest()[:8], 16) / 0xffffffff
    return h < fraction


def main():
    het_pos, het_alt, var_pos, new_base = int(sys.argv[1]), sys.argv[2], int(sys.argv[3]), sys.argv[4]
    fraction = float(sys.argv[5]) if len(sys.argv) > 5 else 1.0
    del_len = int(new_base.split(':')[1]) if new_base.startswith('del:') else 0
    injected = skipped_wrong_hap = skipped_no_cover = skipped_subclonal = 0
    for line in sys.stdin:
        if line.startswith('@'):
            sys.stdout.write(line); continue
        f = line.rstrip('\n').split('\t')
        pos, cigar, seq = int(f[3]), f[5], f[9]
        if cigar == '*' or int(f[1]) & 0x904:               # unmapped/secondary/supplementary
            sys.stdout.write(line); continue
        ops = [(int(n), op) for n, op in CIG.findall(cigar)]

        qh = query_offset(pos, ops, het_pos)
        qv = query_offset(pos, ops, var_pos)
        if qh is None or qv is None:
            skipped_no_cover += 1; sys.stdout.write(line); continue
        if seq[qh].upper() != het_alt.upper():              # the other haplotype
            skipped_wrong_hap += 1; sys.stdout.write(line); continue
        if not del_len and seq[qv].upper() == new_base.upper():
            sys.stdout.write(line); continue                # already that base
        if not keep(f[0], fraction):
            skipped_subclonal += 1; sys.stdout.write(line); continue

        if del_len:
            new_ops = delete_at(ops, pos, var_pos, del_len)
            if new_ops is None:
                skipped_no_cover += 1; sys.stdout.write(line); continue
            f[9] = seq[:qv] + seq[qv+del_len:]
            if f[10] != '*': f[10] = f[10][:qv] + f[10][qv+del_len:]
            f[5] = ''.join(f"{n}{op}" for n, op in new_ops)
        else:
            f[9] = seq[:qv] + new_base + seq[qv+1:]
            f[5] = ''.join(f"{n}{op}" for n, op in split_match_at(ops, pos, var_pos))
        # NM/MD would now be stale; drop them rather than leave them wrong.
        f = [x for x in f if not x.startswith(('NM:i:', 'MD:Z:', 'cs:Z:'))]
        sys.stdout.write('\t'.join(f) + '\n')
        injected += 1
    print(f"injected into {injected} reads (wrong haplotype {skipped_wrong_hap}, "
          f"not covering {skipped_no_cover}, held back as subclonal {skipped_subclonal})", file=sys.stderr)


if __name__ == '__main__':
    main()
