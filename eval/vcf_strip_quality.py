#!/usr/bin/env python3
"""Print a VCF's records with every quality-model field removed.

The quality model (src/var_qual.c) only scores records: it must not change
which records are written, their alleles, FILTER, classification INFO or
genotypes. Diffing two runs through this script checks exactly that.

Removed: QUAL; INFO QLIK SOR MQRS BQRS RU RPL QF; FORMAT GQ PL UPR UPF.

    vcf_strip_quality.py a.vcf > a.txt; vcf_strip_quality.py b.vcf > b.txt; diff a.txt b.txt
"""
import gzip
import sys

INFO_DROP = {"QLIK", "SOR", "MQRS", "BQRS", "RU", "RPL", "QF"}
FMT_DROP = {"GQ", "PL", "UPR", "UPF"}


def main(fn):
    op = gzip.open if fn.endswith(".gz") else open
    with op(fn, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            f[5] = "."
            f[7] = ";".join(x for x in f[7].split(";") if x.split("=", 1)[0] not in INFO_DROP)
            if len(f) > 9:
                keys = f[8].split(":")
                keep = [i for i, k in enumerate(keys) if k not in FMT_DROP]
                f[8] = ":".join(keys[i] for i in keep)
                for s in range(9, len(f)):
                    vals = f[s].split(":")
                    f[s] = ":".join(vals[i] for i in keep if i < len(vals))
            print("\t".join(f))


if __name__ == "__main__":
    main(sys.argv[1])
