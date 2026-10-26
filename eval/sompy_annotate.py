#!/usr/bin/env python3
"""Annotate a query VCF with som.py's verdicts.

som.py (hap.py v0.3.9, --feature-table generic) writes <prefix>.features.csv
with one row per scored record and a tag: TP and FP for query records, FN for
truth records the query missed. This writes them back onto VCFs:

  <out>.vcf.gz      every record of the query VCF; scored ones carry
                    INFO/SOMPY=TP|FP, unscored ones (not PASS+SOMATIC, outside
                    the BED, longer than the length cap) carry nothing
  <out>.fn.vcf.gz   the truth records som.py counted as FN, INFO/SOMPY=FN
  <out>.tsv         one row per scored record: tag, QUAL, QLIK, class, ... for
                    curves (FN rows have no QUAL)

Matching is on (CHROM, POS, REF, ALT), then on (CHROM, POS) for records whose
alleles som.py re-represented. Unmatched counts are reported, not hidden.

    sompy_annotate.py query.vcf(.gz) sompy.features.csv truth.vcf.gz out_prefix
"""
import csv
import gzip
import subprocess
import sys


def vopen(fn):
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)


def bgzip_index(path_txt, out_gz):
    with open(out_gz, "wb") as fo:
        subprocess.run(["bgzip", "-c", path_txt], stdout=fo, check=True)
    subprocess.run(["tabix", "-f", "-p", "vcf", out_gz], check=True)


def main(query, feat, truth, out):
    tags, tags_pos, fn_keys = {}, {}, set()
    with open(feat) as fh:
        for r in csv.DictReader(fh):
            if r["tag"] == "FN":
                fn_keys.add((r["CHROM"], int(r["POS"]), r["REF.truth"] or r["REF"], r["ALT.truth"] or r["ALT"]))
                continue
            k = (r["CHROM"], int(r["POS"]), r["REF"], r["ALT"])
            tags[k] = r["tag"]
            tags_pos.setdefault((r["CHROM"], int(r["POS"])), []).append(r["tag"])
    n_scored = len(tags)
    used, used_pos, by_pos = set(), set(), 0
    # Positions where some query record carries exactly the alleles som.py scored.
    exact_pos = set()
    with vopen(query) as fi:
        for line in fi:
            if line.startswith("#"):
                continue
            f = line.split("\t", 5)
            if (f[0], int(f[1]), f[3], f[4]) in tags:
                exact_pos.add((f[0], int(f[1])))
    tmp_q, tmp_fn = out + ".vcf", out + ".fn.vcf"
    tsv = open(out + ".tsv", "w")
    tsv.write("chrom\tpos\tref\talt\ttag\tqual\tqlik\tqual_old\tsomatic\tunphased\tfilter\n")
    with vopen(query) as fi, open(tmp_q, "w") as fo:
        for line in fi:
            if line.startswith("##"):
                fo.write(line); continue
            if line.startswith("#CHROM"):
                fo.write('##INFO=<ID=SOMPY,Number=1,Type=String,Description="som.py verdict: TP or FP (FN records are in the .fn.vcf.gz companion)">\n')
                fo.write(line); continue
            f = line.rstrip("\n").split("\t")
            k = (f[0], int(f[1]), f[3], f[4])
            tag = tags.get(k)
            if tag is None:
                # Position-only fallback, for alleles som.py re-represented: only
                # where som.py scored nothing at this position under the record's
                # own alleles, and only once -- an indel and a substitution can
                # share a POS, and only one of them was scored.
                pos = (f[0], int(f[1]))
                pt = tags_pos.get(pos)
                if pt and len(set(pt)) == 1 and pos not in exact_pos and pos not in used_pos:
                    tag = pt[0]; by_pos += 1; used_pos.add(pos)
            if tag is not None:
                used.add(k)
                f[7] = f[7] + ";SOMPY=" + tag if f[7] != "." else "SOMPY=" + tag
                info = f[7].split(";")
                ql = next((x[5:] for x in info if x.startswith("QLIK=")), "")
                # The pre-model QUAL, from the summed sample AD (fixed 0.001 error, capped at 60).
                fmt = f[8].split(":"); rn = an = 0
                for s in f[9:]:
                    ad = dict(zip(fmt, s.split(":"))).get("AD", "0,0").split(",")
                    rn += int(ad[0]); an += int(ad[1]) if len(ad) > 1 else 0
                old = min(60, int(-10 * (rn * -0.000434512 + an * -3.0)))
                tsv.write("\t".join([f[0], f[1], f[3], f[4], tag, f[5], ql, str(old), str(int("SOMATIC" in info)),
                                     str(int("UNPHASED" in info)), f[6]]) + "\n")
            fo.write("\t".join(f) + "\n")
    n_fn = 0
    with vopen(truth) as fi, open(tmp_fn, "w") as fo:
        for line in fi:
            if line.startswith("##"):
                fo.write(line); continue
            if line.startswith("#CHROM"):
                fo.write('##INFO=<ID=SOMPY,Number=1,Type=String,Description="som.py verdict: FN">\n')
                fo.write(line); continue
            f = line.rstrip("\n").split("\t")
            hit = any((f[0], int(f[1]), f[3], a) in fn_keys for a in f[4].split(","))
            if hit:
                f[7] = f[7] + ";SOMPY=FN" if f[7] != "." else "SOMPY=FN"
                fo.write("\t".join(f) + "\n"); n_fn += 1
                tsv.write("\t".join([f[0], f[1], f[3], f[4], "FN", "", "", "", "", "", ""]) + "\n")
    tsv.close()
    bgzip_index(tmp_q, out + ".vcf.gz")
    bgzip_index(tmp_fn, out + ".fn.vcf.gz")
    import os
    os.remove(tmp_q); os.remove(tmp_fn)
    n_tagged = len(used)
    sys.stderr.write("som.py scored %d query records: %d annotated (%d by position only), %d unmatched; "
                     "%d of %d FN written\n" % (n_scored, n_tagged, by_pos, n_scored - n_tagged, n_fn, len(fn_keys)))


if __name__ == "__main__":
    if len(sys.argv) != 5:
        sys.exit(__doc__)
    main(*sys.argv[1:])
