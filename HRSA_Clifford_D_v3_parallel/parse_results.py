#!/usr/bin/env python3
"""
Parse output files of the form:
    out_a<angle>_e<epsilon>_f<something>

The angle may use either a dot or an underscore as the decimal separator:
    out_a1.2475730576_e1_f6        (dot form)
    out_a0_7107987091_e0.01_f12   (underscore form → 0.7107987091)

Extracts from each file:
  - angle        : rotation angle parsed from filename
  - target_eps   : epsilon parsed from filename
  - final_f      : last 'f = N' value enumerated in the file
  - eps_diff_val : the matrix Frobenius distance ||X_(0,1) H - R||_F
                   (falls back to the vector-level "Epsilon Diff. Val." if the
                   matrix line is absent, e.g. in truncated/in-progress files)
  - result       : "Success", "Failure", or None if still running

Usage:
    python parse_results.py [data_dir] [--output results.csv]

    data_dir defaults to ./data
"""

import os
import re
import csv
import sys
import argparse


# ── filename parsing ──────────────────────────────────────────────────────────

# The angle field may use '.' or '_' as its decimal separator.
# We match everything between 'out_a' and '_e' as the raw angle string,
# then normalise underscores to dots.
FNAME_RE = re.compile(
    r"out_a(?P<angle>-?[0-9]+(?:[._][0-9]+)?)_e(?P<target_eps>[0-9]+(?:\.[0-9]+)?)_f"
)

def parse_filename(name):
    """
    Return (angle_float, target_eps_str) from the base filename, or (None, None).
    Handles both 'out_a1.234_e0.1_f4' and 'out_a1_234_e0.1_f4' angle formats.
    """
    m = FNAME_RE.search(os.path.basename(name))
    if not m:
        return None, None
    # Normalise underscore decimal separator to dot
    angle_str = m.group("angle").replace("_", ".")
    return angle_str, m.group("target_eps")


# ── file-content parsing ──────────────────────────────────────────────────────

F_LINE_RE    = re.compile(r"f\s*=\s*(\d+)")
MATRIX_RE    = re.compile(r"Matrix Frobenius distance.*?=\s*([0-9eE+\-.]+)")
EPS_DIFF_RE  = re.compile(r"Epsilon Diff\. Val\.\s*:\s*([0-9eE+\-.]+)")
SUCCESS_RE   = re.compile(r"\bSuccess\b", re.IGNORECASE)
FAILURE_RE   = re.compile(r"\bFailure\b", re.IGNORECASE)

def parse_file(path):
    """
    Returns a dict with keys: final_f, eps_diff_val, result.

    eps_diff_val is taken from the matrix Frobenius distance line
    (||X_(0,1) H - R_(0,1)^Z(theta)||_F) in preference to the vector-level
    Epsilon Diff. Val., because the matrix norm is the true end-to-end check.
    The vector value is used as a fallback for truncated/in-progress files
    that don't yet have a matrix Frobenius line.

    Any field that cannot be found is set to None.
    """
    final_f      = None
    matrix_frob  = None
    vec_eps_diff = None
    result       = None

    with open(path, "r", errors="replace") as fh:
        for line in fh:
            # Track every 'f = N' line; last one wins
            for m in F_LINE_RE.finditer(line):
                final_f = int(m.group(1))

            # Matrix Frobenius distance (preferred)
            m = MATRIX_RE.search(line)
            if m:
                matrix_frob = float(m.group(1))

            # Vector-level epsilon diff (fallback)
            m = EPS_DIFF_RE.search(line)
            if m:
                vec_eps_diff = float(m.group(1))

            # Success / failure
            if SUCCESS_RE.search(line):
                result = "Success"
            elif FAILURE_RE.search(line):
                result = "Failure"

    # Use matrix norm if available, otherwise fall back to vector eps_diff
    eps_diff_val = matrix_frob if matrix_frob is not None else vec_eps_diff

    return {
        "final_f":      final_f,
        "eps_diff_val": eps_diff_val,
        "result":       result,
    }


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="Parse simulation output files.")
    parser.add_argument("data_dir", nargs="?", default="data",
                        help="Directory containing output files (default: ./data)")
    parser.add_argument("--output", default="results.csv",
                        help="CSV output path (default: results.csv)")
    args = parser.parse_args()

    data_dir = args.data_dir
    if not os.path.isdir(data_dir):
        print(f"Error: '{data_dir}' is not a directory.", file=sys.stderr)
        sys.exit(1)

    rows = []
    skipped = []

    for fname in sorted(os.listdir(data_dir)):
        fpath = os.path.join(data_dir, fname)
        if not os.path.isfile(fpath):
            continue

        angle, target_eps = parse_filename(fname)
        if angle is None:
            skipped.append(fname)
            continue

        info = parse_file(fpath)
        rows.append({
            "filename":     fname,
            "angle":        angle,
            "target_eps":   target_eps,
            "final_f":      info["final_f"],
            "eps_diff_val": info["eps_diff_val"],
            "result":       info["result"],
        })

    fieldnames = ["angle", "target_eps", "final_f", "eps_diff_val", "result"]
    with open(args.output, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter=" ",
                                extrasaction="ignore")
        csvfile.write("# ")
        writer.writeheader()
        for row in rows:
            missing = [f for f in fieldnames if row[f] is None]
            if missing:
                csvfile.write(
                    f"# error - missing: {', '.join(missing)} | "
                    f"angle={row['angle']} target_eps={row['target_eps']} "
                    f"final_f={row['final_f']} eps_diff_val={row['eps_diff_val']} "
                    f"result={row['result']}\n"
                )
            else:
                writer.writerow(row)

    print(f"Parsed {len(rows)} file(s) → {args.output}")
    if skipped:
        print(f"Skipped {len(skipped)} file(s) with unrecognised names: {skipped}")


if __name__ == "__main__":
    main()
