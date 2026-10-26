#!/usr/bin/env python3
"""Label, fit and report the longcallD quality model. Standard library only.

  label   VCF (run with LONGCALLD_QUALCHK=1, normalised with bcftools norm)
          against a truth VCF and BED, one TSV row per scored record.
  fit     the per-technology artifact model from labelled rows; writes weight
          lines for src/var_qual_weights.inc or --qual-weights.
  report  ranking (AUPRC), the F1-optimal QUAL cut and calibration, optionally
          against a second QUAL column (e.g. the previous binary's).

Labels
  somatic   SOMATIC records only; true if the allele is in the somatic truth.
  germline  every record; true if the allele is in the truth (genotype ignored),
            and gt_ok says whether the genotype matches too.

The QF field is the quality class, P_art and the artifact features in the order
of src/var_qual.h (LONGCALLD_VQF_NAMES).
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
            "log_alt", "vaf"]
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
                         "dataset"] + FEATURES + (["qual_old"] if a.old else [])) + "\n")
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
            cls = CLASSES[int(qf[0])]
            if a.mode == "somatic" and "SOMATIC" not in info:
                continue
            key = (c, p, ref, alt)
            y_true = key in truth
            # Germline: a record is right if the allele is; gt_ok asks the genotype too.
            gt_ok = 0
            if a.mode == "germline" and y_true:
                fmt = f[8].split(":")
                smp = dict(zip(fmt, f[9].split(":")))
                gt_ok = int(dosage(smp.get("GT", "./.")) == truth[key])
            n["tp" if y_true else "fp"] += 1
            row = [c, str(p), ref, alt, f[6], cls, str(0 if y_true else 1), str(gt_ok), f[5],
                   info.get("QLIK", "0"), "%g" % qf[1], a.dataset] + ["%g" % x for x in qf[2:]]
            if a.old:
                row.append(old.get(key, "0"))
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
                     "qlik": float(f[ix["qlik"]]), "dataset": f[ix["dataset"]],
                     "filter": f[ix["filter"]],
                     "x": [float(f[ix[n]]) for n in FEATURES]}
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


def p_false(w, x, q):
    s = sigmoid(sum(wi * xi for wi, xi in zip(w, x)))
    return 1.0 - (1.0 - q) * (1.0 - s), s


def q_lik(r):
    return 10.0 ** (-r["qlik"] / 10.0)


def loglik(w, rows, idx):
    ll = 0.0
    for r in rows:
        x = [r["x"][i] for i in idx]
        p, _ = p_false(w, x, q_lik(r))
        p = min(max(p, 1e-12), 1 - 1e-12)
        ll += math.log(p) if r["y"] else math.log(1 - p)
    return ll


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


def fit(rows, idx, l2=1.0, iters=30):
    """Fisher scoring on the combined model P(false) = 1 - (1-q)(1-sigmoid(w.x)).
    The intercept starts at the log-odds of the false rate and is not penalised."""
    d = len(idx)
    rate = (sum(r["y"] for r in rows) + 0.5) / (len(rows) + 1.0)
    w = [0.0] * d
    if 0 in idx:
        w[idx.index(0)] = math.log(rate / (1 - rate))
    ll = loglik(w, rows, idx) - 0.5 * l2 * sum(wi * wi for j, wi in enumerate(w) if idx[j] != 0)
    for _ in range(iters):
        g = [0.0] * d
        H = [[0.0] * d for _ in range(d)]
        for r in rows:
            x = [r["x"][i] for i in idx]
            q = q_lik(r)
            p, s = p_false(w, x, q)
            p = min(max(p, 1e-12), 1 - 1e-12)
            dp = (1 - q) * s * (1 - s)
            gi = (r["y"] - p) / (p * (1 - p)) * dp
            hi = dp * dp / (p * (1 - p))
            for j in range(d):
                if x[j] == 0.0:
                    continue
                g[j] += gi * x[j]
                hx = hi * x[j]
                Hj = H[j]
                for k in range(j, d):
                    if x[k] != 0.0:
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
            w2 = [wi + t * si for wi, si in zip(w, step)]
            ll2 = loglik(w2, rows, idx) - 0.5 * l2 * sum(wi * wi for j, wi in enumerate(w2) if idx[j] != 0)
            if ll2 >= ll - 1e-9:
                break
            t /= 2
        if t <= 1e-4:
            break
        improved = ll2 - ll
        w, ll = w2, ll2
        if improved < 1e-6:
            break
    return w


def per_row_loss(w, rows, idx):
    out = []
    for r in rows:
        x = [r["x"][i] for i in idx]
        p, _ = p_false(w, x, q_lik(r))
        p = min(max(p, 1e-12), 1 - 1e-12)
        out.append(-(math.log(p) if r["y"] else math.log(1 - p)))
    return out


def split(rows, how, win):
    if how == "dataset":
        return rows, rows
    tr = [r for r in rows if (r["pos"] // win) % 2 == 1]
    te = [r for r in rows if (r["pos"] // win) % 2 == 0]
    return tr, te


def cmd_fit(a):
    out = open(a.out, "w") if a.out else sys.stdout
    for cls in a.classes.split(","):
        train = read_rows(a.train, {cls})
        test = read_rows(a.test, {cls}) if a.test else None
        if test is None:
            train, test = split(train, "window", a.win)
        n_false = sum(r["y"] for r in train)
        sys.stderr.write("%s: train %d (%d false), test %d (%d false)\n" %
                         (cls, len(train), n_false, len(test), sum(r["y"] for r in test)))
        if n_false < a.min_false or len(train) - n_false < a.min_false:
            sys.stderr.write("%s: too few of one label to fit; left unset\n" % cls)
            continue
        if a.max_rows and len(train) > a.max_rows:
            # Keep every false record; subsample the true ones and reweight
            # nothing -- the intercept is refit on the full set below.
            rnd = random.Random(1)
            fl = [r for r in train if r["y"]]
            tr = [r for r in train if not r["y"]]
            rnd.shuffle(tr)
            train_fit = fl + tr[:max(a.max_rows - len(fl), len(fl))]
        else:
            train_fit = train
        idx = [i for i, n in enumerate(FEATURES) if n not in set(a.exclude.split(","))]
        # Drop features that do not earn their place on held-out records.
        while True:
            w = fit(train_fit, idx, a.l2)
            base = per_row_loss(w, test, idx)
            worst, worst_gain = None, None
            for j, fi in enumerate(idx):
                if fi == 0:
                    continue
                idx2 = idx[:j] + idx[j+1:]
                w2 = fit(train_fit, idx2, a.l2)
                l2_ = per_row_loss(w2, test, idx2)
                diff = [b2 - b1 for b1, b2 in zip(base, l2_)]  # >0: removing hurts
                m = sum(diff) / len(diff)
                sd = math.sqrt(sum((d - m) ** 2 for d in diff) / max(1, len(diff) - 1))
                se = sd / math.sqrt(len(diff))
                gain = m / se if se > 0 else (1e9 if m > 0 else -1e9)
                if worst is None or gain < worst_gain:
                    worst, worst_gain = j, gain
            if worst is None or worst_gain > 1.0:
                break
            sys.stderr.write("%s: drop %s (held-out gain %.2f SE)\n" % (cls, FEATURES[idx[worst]], worst_gain))
            idx = idx[:worst] + idx[worst+1:]
        w = fit(train_fit, idx, a.l2)
        if train_fit is not train:
            # Restore the intercept to the full training set's false rate.
            for _ in range(50):
                ll_obs = sum(r["y"] for r in train)
                ll_exp = sum(p_false(w, [r["x"][i] for i in idx], q_lik(r))[0] for r in train)
                if abs(ll_obs - ll_exp) < 0.5:
                    break
                w[idx.index(0)] += math.log((ll_obs + 0.5) / (ll_exp + 0.5))
        te_loss = sum(per_row_loss(w, test, idx)) / len(test)
        lik_loss = sum(per_row_loss([-30.0] + [0.0] * (len(idx) - 1), test, idx)) / len(test)
        sys.stderr.write("%s: held-out log-loss %.4f (likelihood only %.4f)\n" % (cls, te_loss, lik_loss))
        for j, fi in enumerate(idx):
            out.write("%s %s %s %.6g\n" % (a.tech, cls, FEATURES[fi], w[j]))


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


def cmd_report(a):
    rows = read_rows(a.rows, set(a.classes.split(",")) if a.classes else None)
    if a.pass_only:
        rows = [r for r in rows if r["filter"] == "PASS"]
    if a.half:
        rows = split(rows, "window", a.win)[1 if a.half == "even" else 0]
    keys = ["qual"] + [k for k in (rows[0].keys() if rows else []) if k.startswith("qual_")]
    print("records %d  false %d" % (len(rows), sum(r["y"] for r in rows)))
    for k in keys:
        ap = auprc([(r[k], r["y"]) for r in rows])
        line = "%-10s AUPRC %.4f" % (k, ap)
        if a.n_truth:
            f1, cut, tp, fp = best_cut(rows, k, a.n_truth)
            line += "   best F1 %.4f at >= %s (TP %d FP %d)" % (f1, cut, tp, fp)
        print(line)
    # Calibration: observed false rate per QUAL bin against 10^(-Q/10).
    print("\nQUAL bin   n      false  observed  predicted(mean)")
    edges = [0, 5, 10, 15, 20, 25, 30, 40, 50, 60, 80, 100, 1000]
    for lo, hi in zip(edges, edges[1:]):
        b = [r for r in rows if lo <= r["qual"] < hi]
        if not b:
            continue
        nf = sum(r["y"] for r in b)
        pred = sum(10 ** (-r["qual"] / 10) for r in b) / len(b)
        print("[%3d,%4d) %6d %6d   %.4f    %.4f" % (lo, hi, len(b), nf, nf / len(b), pred))
    if a.gq:
        print("\nGT accuracy by GQ is in gq_report; QUAL rows here use the record label")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("label")
    p.add_argument("vcf"); p.add_argument("truth")
    p.add_argument("--bed"); p.add_argument("--region"); p.add_argument("--out")
    p.add_argument("--mode", choices=["somatic", "germline"], required=True)
    p.add_argument("--dataset", default="ds")
    p.add_argument("--old", help="the same run from another binary; its QUAL becomes column qual_old")
    p = sub.add_parser("fit")
    p.add_argument("--train", nargs="+", required=True)
    p.add_argument("--test", nargs="+")
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
    a = ap.parse_args()
    {"label": cmd_label, "fit": cmd_fit, "report": cmd_report}[a.cmd](a)


if __name__ == "__main__":
    main()
