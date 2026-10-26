#!/usr/bin/env python3
"""Label, fit and report the longcallD quality model. Standard library only.

  label   VCF (run with LONGCALLD_QUALCHK=1, normalised with bcftools norm)
          against a truth VCF and BED, one TSV row per scored record.
  fit     the per-technology model from labelled rows; writes weight lines for
          src/var_qual_weights.inc or --qual-weights.
  report  ranking (AUPRC), the F1-optimal QUAL cut and calibration, optionally
          against a second QUAL column (the previous binary's) and against QUAL
          recomputed offline from a weight file.

Labels
  somatic   SOMATIC records only; true if the allele is in the somatic truth.
  germline  every record; true if the allele is in the truth (genotype ignored),
            and gt_ok says whether the genotype matches too.

The model, as src/var_qual.c evaluates it, per technology and class:
    logit P(false) = w . x
where x is the QF features plus lik = min(QLIK, 60)/10 and lik_low =
min(QLIK, 20)/10, which are recomputed here from QLIK. The QF field is the
quality class, the model's P(false) and the features in the order of
src/var_qual.h (LONGCALLD_VQF_NAMES).
"""
import argparse
import bisect
import gzip
import math
import random
import sys

FEATURES = ["bias", "sor", "one_strand", "mq_rs", "mq_diff", "bq_rs", "bq_low",
            "hp_indel", "str_indel", "hp_snv", "unph_frac", "unph_alt_excess", "log_unph",
            "alt_minor_hap", "alt_conflict", "unphased_call", "dp_high", "dp_low",
            "log_alt", "vaf", "lik", "lik_low"]
I_LIK, I_LIK_LOW = FEATURES.index("lik"), FEATURES.index("lik_low")
CLASSES = ["germ_snv", "germ_indel", "som_snv", "som_indel"]


def vopen(fn):
    return gzip.open(fn, "rt") if fn.endswith(".gz") else open(fn)


def read_bed(fn):
    regs = {}
    with vopen(fn) as fh:
        for line in fh:
            if line.startswith(("#", "track", "browser")) or not line.strip():
                continue
            c, b, e = line.split("\t")[:3]
            regs.setdefault(c, []).append((int(b), int(e)))
    out = {}
    for c, iv in regs.items():
        iv.sort()
        merged = []
        for b, e in iv:
            if merged and b <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], e)
            else:
                merged.append([b, e])
        out[c] = ([b for b, _ in merged], [e for _, e in merged])
    return out


def in_bed(bed, c, pos1):
    if bed is None:
        return True
    if c not in bed:
        return False
    starts, ends = bed[c]
    i = bisect.bisect_right(starts, pos1 - 1) - 1
    return i >= 0 and ends[i] >= pos1


def dosage(gt):
    a = gt.replace("|", "/").split("/")
    if "." in a:
        return -1
    return sum(1 for x in a if x != "0")


def parse_region(s):
    if s is None:
        return None
    c, _, r = s.partition(":")
    if not r:
        return (c, 0, 1 << 62)
    b, _, e = r.replace(",", "").partition("-")
    return (c, int(b), int(e))


def in_region(c, p, reg):
    return c == reg[0] and reg[1] <= p <= reg[2]


def read_truth(fn, region=None):
    t = {}
    with vopen(fn) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if region and not in_region(f[0], int(f[1]), region):
                continue
            gt = f[9].split(":")[0] if len(f) > 9 else "0/1"
            for alt in f[4].split(","):
                t[(f[0], int(f[1]), f[3], alt)] = dosage(gt)
    return t


def cmd_label(a):
    bed = read_bed(a.bed) if a.bed else None
    reg = parse_region(a.region)
    truth = read_truth(a.truth, reg)
    old = {}
    if a.old:  # the same calls from another binary, for its QUAL
        with vopen(a.old) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                f = line.split("\t", 6)
                old[(f[0], int(f[1]), f[3], f[4])] = f[5]
    out = open(a.out, "w") if a.out else sys.stdout
    out.write("\t".join(["chrom", "pos", "ref", "alt", "filter", "cls", "y", "gt_ok", "qual", "qlik", "p_art",
                         "dataset", "gq", "gq_raw"] + FEATURES + (["qual_old"] if a.old or a.old_formula else [])) + "\n")
    n = {"tp": 0, "fp": 0}
    with vopen(a.vcf) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            c, p, ref, alt = f[0], int(f[1]), f[3], f[4]
            if reg and not in_region(c, p, reg):
                continue
            if not in_bed(bed, c, p):
                continue
            info = dict(kv.split("=", 1) if "=" in kv else (kv, "1") for kv in f[7].split(";"))
            if "QF" not in info:
                sys.exit("label: record without QF; run the caller with LONGCALLD_QUALCHK=1")
            qf = [float(x) for x in info["QF"].split(",")]
            qf += [0.0] * (2 + len(FEATURES) - len(qf))  # a binary from before lik/lik_low; read_rows derives them
            cls = CLASSES[int(qf[0])]
            if a.mode == "somatic" and "SOMATIC" not in info:
                continue
            key = (c, p, ref, alt)
            y_true = key in truth
            # Germline: a record is right if the allele is; gt_ok asks the genotype too.
            gt_ok = 0
            smp = dict(zip(f[8].split(":"), f[-1].split(":")))  # single sample, or the tumour column
            gq = smp.get("GQ", "0")
            # GQ before the 99 cap, from PL: the best other genotype minus the printed one.
            gq_raw = "0"
            if "PL" in smp and dosage(smp.get("GT", "./.")) >= 0:
                pl = [int(x) for x in smp["PL"].split(",")]
                d = dosage(smp["GT"])
                gq_raw = str(max(0, min(pl[g] for g in range(3) if g != d) - pl[d]))
            if a.mode == "germline" and y_true:
                gt_ok = int(dosage(smp.get("GT", "./.")) == truth[key])
            n["tp" if y_true else "fp"] += 1
            row = [c, str(p), ref, alt, f[6], cls, str(0 if y_true else 1), str(gt_ok), f[5],
                   info.get("QLIK", "0"), "%g" % qf[1], a.dataset, gq, gq_raw] + ["%g" % x for x in qf[2:]]
            if a.old:
                row.append(old.get(key, "0"))
            elif a.old_formula:
                # The pre-model QUAL: fixed 0.001 error on the record's ref and alt
                # depth, truncated, capped at 60. Paired AD is per sample here, so
                # it is summed over the columns.
                fmt = f[8].split(":")
                ref_n = alt_n = 0
                for s in f[9:]:
                    ad = dict(zip(fmt, s.split(":"))).get("AD", "0,0").split(",")
                    ref_n += int(ad[0]); alt_n += int(ad[1]) if len(ad) > 1 else 0
                row.append(str(min(60, int(-10 * (ref_n * math.log10(0.999) + alt_n * -3.0)))))
            out.write("\t".join(row) + "\n")
    sys.stderr.write("label: %d true, %d false records\n" % (n["tp"], n["fp"]))


def read_rows(fns, cls=None):
    rows = []
    for fn in fns:
        with open(fn) as fh:
            hdr = fh.readline().rstrip("\n").split("\t")
            ix = {h: i for i, h in enumerate(hdr)}
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if cls and f[ix["cls"]] not in cls:
                    continue
                r = {"chrom": f[ix["chrom"]], "pos": int(f[ix["pos"]]), "cls": f[ix["cls"]],
                     "y": int(f[ix["y"]]), "gt_ok": int(f[ix["gt_ok"]]), "qual": float(f[ix["qual"]]),
                     "qlik": float(f[ix["qlik"]]), "dataset": f[ix["dataset"]], "gq": int(f[ix["gq"]]),
                     "gq_raw": int(f[ix["gq_raw"]]) if "gq_raw" in ix else int(f[ix["gq"]]),
                     "filter": f[ix["filter"]],
                     "x": [float(f[ix[n]]) if n in ix else 0.0 for n in FEATURES]}
                r["x"][I_LIK] = min(r["qlik"], 60.0) / 10.0
                r["x"][I_LIK_LOW] = min(r["qlik"], 20.0) / 10.0
                for k in ix:
                    if k.startswith("qual_"):
                        r[k] = float(f[ix[k]])
                rows.append(r)
    return rows


# ------------------------------------------------------------------ model --

def sigmoid(z):
    if z >= 0:
        return 1.0 / (1.0 + math.exp(-z))
    e = math.exp(z)
    return e / (1.0 + e)


def clampp(p):
    return min(max(p, 1e-12), 1 - 1e-12)


def predict(w, x):
    return sigmoid(sum(wi * xi for wi, xi in zip(w, x)))


def objective(w, rows, idx, l2):
    ll = 0.0
    for r in rows:
        p = clampp(predict(w, [r["x"][i] for i in idx]))
        ll += r["w"] * (math.log(p) if r["y"] else math.log(1 - p))
    return ll - 0.5 * l2 * sum(w[j] ** 2 for j in range(len(idx)) if idx[j] != 0)


def solve(A, b):
    n = len(b)
    M = [row[:] + [b[i]] for i, row in enumerate(A)]
    for c in range(n):
        piv = max(range(c, n), key=lambda r: abs(M[r][c]))
        M[c], M[piv] = M[piv], M[c]
        if abs(M[c][c]) < 1e-12:
            continue
        for r in range(n):
            if r != c:
                fct = M[r][c] / M[c][c]
                if fct:
                    for k in range(c, n + 1):
                        M[r][k] -= fct * M[c][k]
    return [M[i][n] / M[i][i] if abs(M[i][i]) > 1e-12 else 0.0 for i in range(n)]


def fit(rows, idx, l2=1.0, iters=50, start=None):
    """Weighted, ridge-penalised logistic regression by Newton's method (the
    intercept is not penalised). Convex, so a warm start changes speed only."""
    d = len(idx)
    if start is not None:
        w = start[:]
    else:
        w = [0.0] * d
        rate = (sum(r["w"] * r["y"] for r in rows) + 0.5) / (sum(r["w"] for r in rows) + 1.0)
        if 0 in idx:
            w[idx.index(0)] = math.log(rate / (1 - rate))
    obj = objective(w, rows, idx, l2)
    for _ in range(iters):
        g = [0.0] * d
        H = [[0.0] * d for _ in range(d)]
        for r in rows:
            x = [r["x"][i] for i in idx]
            p = predict(w, x)
            gi = r["w"] * (r["y"] - p)
            hi = r["w"] * p * (1 - p)
            nz = [j for j in range(d) if x[j] != 0.0]
            for j in nz:
                g[j] += gi * x[j]
                hx = hi * x[j]
                Hj = H[j]
                for k in nz:
                    if k >= j:
                        Hj[k] += hx * x[k]
        for j in range(d):
            for k in range(j):
                H[j][k] = H[k][j]
            if idx[j] != 0:
                g[j] -= l2 * w[j]
                H[j][j] += l2
            H[j][j] += 1e-9
        step = solve(H, g)
        t = 1.0
        while t > 1e-4:
            w2 = [a + t * s for a, s in zip(w, step)]
            obj2 = objective(w2, rows, idx, l2)
            if obj2 >= obj - 1e-9:
                break
            t /= 2
        if t <= 1e-4:
            break
        gain = obj2 - obj
        w, obj = w2, obj2
        if gain < 1e-7:
            break
    return w


def per_row_loss(w, rows, idx):
    out = []
    for r in rows:
        p = clampp(predict(w, [r["x"][i] for i in idx]))
        out.append(-(math.log(p) if r["y"] else math.log(1 - p)))
    return out


def lik_only_loss(rows):
    out = []
    for r in rows:
        p = clampp(10 ** (-r["qlik"] / 10))
        out.append(-(math.log(p) if r["y"] else math.log(1 - p)))
    return out


def split(rows, how, win):
    if how == "dataset":
        return rows, rows
    tr = [r for r in rows if (r["pos"] // win) % 2 == 1]
    te = [r for r in rows if (r["pos"] // win) % 2 == 0]
    return tr, te


def subsample(rows, max_rows, seed=1):
    """Every false record, and true ones at rate r with weight 1/r, so the fit
    stays unbiased while the pure-Python Newton step stays affordable."""
    fl = [r for r in rows if r["y"]]
    tr = [r for r in rows if not r["y"]]
    for r in fl:
        r["w"] = 1.0
    keep = max(max_rows - len(fl), min(len(tr), 5 * len(fl)))
    if keep >= len(tr):
        for r in tr:
            r["w"] = 1.0
        return fl + tr
    wt = len(tr) / float(keep)
    tr = random.Random(seed).sample(tr, keep)
    for r in tr:
        r["w"] = wt
    return fl + tr


def cmd_fit(a):
    out = open(a.out, "w") if a.out else sys.stdout
    always = {0, I_LIK, I_LIK_LOW}  # intercept and the likelihood are never dropped
    for cls in a.classes.split(","):
        train = read_rows(a.train, {cls})
        test = read_rows(a.test, {cls}) if a.test else None
        if test is None:
            train, test = split(train, "window", a.win)
        for r in test:
            r["w"] = 1.0
        n_false = sum(r["y"] for r in train)
        sys.stderr.write("%s: train %d (%d false), test %d (%d false)\n" %
                         (cls, len(train), n_false, len(test), sum(r["y"] for r in test)))
        if n_false < a.min_false or len(train) - n_false < a.min_false or not test:
            sys.stderr.write("%s: too few of one label to fit; left unset\n" % cls)
            continue
        tfit = subsample(train, a.max_rows)
        idx = [i for i, nm in enumerate(FEATURES) if nm not in set(a.exclude.split(","))
               and (i in always or any(r["x"][i] != 0.0 for r in tfit))]
        # Drop, in rounds, every feature whose removal does not cost at least one
        # standard error of held-out log-loss. Each must earn its place.
        for rnd in range(5):
            w = fit(tfit, idx, a.l2)
            base = per_row_loss(w, test, idx)
            drop = []
            for j, fi in enumerate(idx):
                if fi in always:
                    continue
                idx2 = idx[:j] + idx[j + 1:]
                w2 = fit(tfit, idx2, a.l2, start=w[:j] + w[j + 1:])
                diff = [b2 - b1 for b1, b2 in zip(base, per_row_loss(w2, test, idx2))]  # > 0: removing hurts
                m = sum(diff) / len(diff)
                sd = math.sqrt(sum((x - m) ** 2 for x in diff) / max(1, len(diff) - 1))
                se = sd / math.sqrt(len(diff))
                # A feature with no effect differs only by rounding; its ratio
                # to an equally tiny SE means nothing.
                gain = 0.0 if abs(m) < 1e-9 else m / se
                if gain < 1.0:
                    drop.append((FEATURES[fi], gain))
            if not drop:
                break
            sys.stderr.write("%s: round %d drops %s\n" % (cls, rnd + 1,
                             ", ".join("%s (%.2f SE)" % t for t in drop)))
            gone = {nm for nm, _ in drop}
            idx = [i for i in idx if FEATURES[i] not in gone]
        w = fit(tfit, idx, a.l2)
        te_loss = sum(per_row_loss(w, test, idx)) / len(test)
        lik_loss = sum(lik_only_loss(test)) / len(test)
        sys.stderr.write("%s: held-out log-loss %.4f (QLIK alone %.4f); kept %s\n" %
                         (cls, te_loss, lik_loss, ",".join(FEATURES[i] for i in idx)))
        for j, fi in enumerate(idx):
            out.write("%s %s %s %.6g\n" % (a.tech, cls, FEATURES[fi], w[j]))


def read_weights(fn, tech):
    """{class: weights by feature index} for one technology."""
    m = {}
    with open(fn) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            t, c, f, w = line.split()
            if t == tech:
                m.setdefault(c, [0.0] * len(FEATURES))[FEATURES.index(f)] = float(w)
    return m


def qual_from_weights(r, m):
    if r["cls"] not in m:
        return min(999.0, round(r["qlik"]))
    z = sum(wi * xi for wi, xi in zip(m[r["cls"]], r["x"]))
    ln_p = z - math.log1p(math.exp(z)) if z < 0 else -math.log1p(math.exp(-z))
    return min(999.0, max(0.0, round(-10 * ln_p / math.log(10))))


# ----------------------------------------------------------------- report --

def auprc(scores_y):
    """Average precision, ranking by score descending (true = y == 0)."""
    s = sorted(scores_y, key=lambda t: -t[0])
    n_true = sum(1 for _, y in s if y == 0)
    if n_true == 0:
        return float("nan")
    tp = fp = 0
    ap = 0.0
    i = 0
    while i < len(s):
        j = i
        while j < len(s) and s[j][0] == s[i][0]:
            j += 1
        dtp = sum(1 for k in range(i, j) if s[k][1] == 0)
        tp += dtp
        fp += (j - i) - dtp
        ap += dtp / n_true * (tp / (tp + fp))
        i = j
    return ap


def best_cut(rows, key, n_truth):
    s = sorted(rows, key=lambda r: -r[key])
    best = (0.0, None, 0, 0)
    tp = fp = 0
    for i, r in enumerate(s):
        if r["y"] == 0:
            tp += 1
        else:
            fp += 1
        if i + 1 < len(s) and s[i + 1][key] == r[key]:
            continue
        rec = tp / n_truth if n_truth else 0
        prec = tp / (tp + fp)
        f1 = 2 * prec * rec / (prec + rec) if prec + rec else 0
        if f1 > best[0]:
            best = (f1, r[key], tp, fp)
    return best


def calibration(rows, key):
    print("\n%-9s  n      false  observed  predicted(mean)" % (key + " bin"))
    edges = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 1000]
    for lo, hi in zip(edges, edges[1:]):
        b = [r for r in rows if lo <= r[key] < hi]
        if not b:
            continue
        nf = sum(r["y"] for r in b)
        pred = sum(10 ** (-r[key] / 10) for r in b) / len(b)
        print("[%3d,%4d) %6d %6d   %.4f    %.4f" % (lo, hi, len(b), nf, nf / len(b), pred))


def cmd_report(a):
    rows = read_rows(a.rows, set(a.classes.split(",")) if a.classes else None)
    if a.pass_only:
        rows = [r for r in rows if r["filter"] == "PASS"]
    if a.half:
        rows = split(rows, "window", a.win)[1 if a.half == "even" else 0]
    if a.weights:
        m = read_weights(a.weights, a.tech)
        for r in rows:
            r["qual_fit"] = qual_from_weights(r, m)
    keys = ["qual"] + [k for k in (rows[0].keys() if rows else []) if k.startswith("qual_")]
    print("records %d  false %d" % (len(rows), sum(r["y"] for r in rows)))
    for k in keys:
        ap = auprc([(r[k], r["y"]) for r in rows])
        line = "%-10s AUPRC %.4f" % (k, ap)
        if a.n_truth:
            f1, cut, tp, fp = best_cut(rows, k, a.n_truth)
            line += "   best F1 %.4f at >= %s (TP %d FP %d)" % (f1, cut, tp, fp)
        print(line)
    for k in ["qual"] + (["qual_fit"] if a.weights else []):
        calibration(rows, k)
    if a.gq:
        # Genotype calibration, on records whose allele is real: GT wrong given GQ.
        t = [r for r in rows if r["y"] == 0]
        print("\nGQ bin     n      GT wrong  observed  predicted(mean)")
        for lo, hi in zip([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 99], [10, 20, 30, 40, 50, 60, 70, 80, 90, 99, 100]):
            b = [r for r in t if lo <= r["gq"] < hi]
            if not b:
                continue
            nw = sum(1 - r["gt_ok"] for r in b)
            pred = sum(10 ** (-r["gq"] / 10) for r in b) / len(b)
            print("[%2d,%3d) %7d %8d   %.4f    %.4f" % (lo, hi, len(b), nw, nw / len(b), pred))


def gq_table(t, scale):
    print("GQ bin     n      GT wrong  observed  predicted(mean)   [PL x %.2f]" % scale)
    for lo, hi in zip([0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 99], [10, 20, 30, 40, 50, 60, 70, 80, 90, 99, 100]):
        b = [r for r in t if lo <= min(99, round(scale * r["gq_raw"])) < hi]
        if not b:
            continue
        nw = sum(1 - r["gt_ok"] for r in b)
        pred = sum(10 ** (-min(99, round(scale * r["gq_raw"])) / 10) for r in b) / len(b)
        print("[%2d,%3d) %7d %8d   %.4f    %.4f" % (lo, hi, len(b), nw, nw / len(b), pred))


def cmd_gqfit(a):
    """Fit the PL tempering factor c (PL' = c PL) on records whose allele is real,
    so GQ' = c GQ predicts genotype error: P(wrong) = 1 / (1 + 10^(c GQ / 10))."""
    train = [r for r in read_rows(a.train) if r["y"] == 0]
    test = [r for r in read_rows(a.test) if r["y"] == 0] if a.test else train

    def nll(c, rows):
        s = 0.0
        for r in rows:
            z = c * r["gq_raw"] / 10 * math.log(10)
            lw = -math.log1p(math.exp(z)) if z < 50 else -z  # ln P(wrong)
            lr = z - math.log1p(math.exp(z)) if z < 50 else 0.0  # ln P(right)
            s -= lw if not r["gt_ok"] else lr
        return s / len(rows)

    best = min((nll(c / 100.0, train), c / 100.0) for c in range(5, 151))
    print("train: c = %.2f (log-loss %.5f; c = 1: %.5f)" % (best[1], best[0], nll(1.0, train)))
    print("test:  log-loss %.5f at c = %.2f, %.5f at c = 1\n" % (nll(best[1], test), best[1], nll(1.0, test)))
    gq_table(test, 1.0)
    print()
    gq_table(test, best[1])


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("label")
    p.add_argument("vcf"); p.add_argument("truth")
    p.add_argument("--bed"); p.add_argument("--region"); p.add_argument("--out")
    p.add_argument("--mode", choices=["somatic", "germline"], required=True)
    p.add_argument("--dataset", default="ds")
    p.add_argument("--old", help="the same run from another binary; its QUAL becomes column qual_old")
    p.add_argument("--old-formula", action="store_true",
                   help="qual_old recomputed from AD with the pre-model formula (an approximation in paired "
                        "mode, where the old QUAL used the pooled candidate depth)")
    p = sub.add_parser("fit")
    p.add_argument("--train", nargs="+", required=True)
    p.add_argument("--test", nargs="+", help="held-out rows; default: the even --win windows of --train")
    p.add_argument("--tech", choices=["hifi", "ont"], required=True)
    p.add_argument("--classes", default="som_snv,som_indel")
    p.add_argument("--win", type=int, default=25_000_000)
    p.add_argument("--l2", type=float, default=1.0)
    p.add_argument("--min-false", type=int, default=20)
    p.add_argument("--max-rows", type=int, default=20000)
    p.add_argument("--exclude", default="")
    p.add_argument("--out")
    p = sub.add_parser("report")
    p.add_argument("rows", nargs="+")
    p.add_argument("--classes")
    p.add_argument("--n-truth", type=int)
    p.add_argument("--pass-only", action="store_true")
    p.add_argument("--half", choices=["odd", "even"])
    p.add_argument("--win", type=int, default=25_000_000)
    p.add_argument("--gq", action="store_true")
    p.add_argument("--weights", help="recompute QUAL from these weights as column qual_fit")
    p.add_argument("--tech", choices=["hifi", "ont"], default="hifi")
    p = sub.add_parser("gqfit")
    p.add_argument("--train", nargs="+", required=True)
    p.add_argument("--test", nargs="+")
    a = ap.parse_args()
    {"label": cmd_label, "fit": cmd_fit, "report": cmd_report, "gqfit": cmd_gqfit}[a.cmd](a)


if __name__ == "__main__":
    main()
