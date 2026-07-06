#!/usr/bin/python3
"""
Aggregate SIP forward-solver run logs into LaTeX tables.

Each file in LOGDIR holds one JSON run record per line; the last line is
taken as the record for that log. For every contrast, two tables are
written to OUTFILE: solver iterations and total timing, with mesh node
counts as rows and the anomaly `ratio` as columns. Progress and warnings
go to stderr.

In the timing tables a second number is appended after "/": the timing
relative to the smallest mesh (first row) in the same column. In the
iterations tables, if the directory LOGDIR0 exists, the matching value
from it is appended after "/" instead (e.g. "3/4").
"""
import os
import sys
import json
import numpy as np

LOGDIR = "logs/"
LOGDIR0 = "logs0/"
OUTFILE = "stats.tex"


def read_logs(logdir):
    """Return the run records (last JSON line of each file in `logdir`)."""
    records = []
    for name in sorted(os.listdir(logdir)):
        path = os.path.join(logdir, name)
        if not os.path.isfile(path):
            continue
        with open(path) as f:
            lines = [ln for ln in f.read().splitlines() if ln.strip()]
        if not lines:
            print(f"warning: {path} is empty, skipped.", file=sys.stderr)
            continue
        try:
            records.append(json.loads(lines[-1]))
        except json.JSONDecodeError as e:
            print(f"warning: cannot parse last line of {path}: {e}", file=sys.stderr)
    return records


def build_index(records):
    """Map (mesh, contrast, ratio) -> record; first record wins on duplicates."""
    index = {}
    for rec in records:
        index.setdefault((rec["mesh"], rec["contrast"], rec["ratio"]), rec)
    return index


def fit_power_law(xs, ys):
    """Least-squares fit y = a * x**b in log-log space; returns (a, b)."""
    b, loga = np.polyfit(np.log(xs), np.log(ys), 1)
    return float(np.exp(loga)), float(b)


def node_count(flyfile):
    """Number of 3D nodes declared in a finley .fly mesh file, or None."""
    with open(flyfile) as f:
        for row in f:
            if row.startswith("3D-Nodes"):
                return int(row[len("3D-Nodes"):])
    return None


def make_table(index, meshes, node_counts, contrast, ratios, field, fmt,
               index0=None, relative=False):
    """Build one LaTeX table for `field` at `contrast`.

    `index` maps (mesh, contrast, ratio) -> record. Missing cells are blank.
    A second number may be appended after "/":
      - if `relative`, the value divided by the value for the smallest mesh
        (`meshes[0]`) in the same column, i.e. the scaling relative to it;
      - else if `index0` is given, its value for the same key.
    """
    header = " # nodes & $N/N_0$" + "".join(f" & {r}" for r in ratios) + "\\\\"
    lines = ["\\toprule", header, "\\midrule"]
    base_nodes = node_counts[meshes[0]]
    for mesh in meshes:
        cells = []
        for r in ratios:
            key = (mesh, contrast, r)
            rec = index.get(key)
            cell = format(rec[field], fmt) if rec is not None else ""
            if relative:
                base = index.get((meshes[0], contrast, r))
                if rec is not None and base is not None and base[field]:
                    cell += "/" + format(rec[field] / base[field], fmt)
                else:
                    cell += "/"
            elif index0 is not None:
                rec0 = index0.get(key)
                cell += "/" + (format(rec0[field], fmt) if rec0 is not None else "")
            cells.append(cell)
        node_ratio = format(node_counts[mesh] / base_nodes, ".3g")
        lines.append(f"{node_counts[mesh]} & {node_ratio}" + "".join(f" & {c}" for c in cells) + " \\\\")
    if relative:
        # per column, fit absolute time  T = a * N**b  in log-log space
        coeffs, exps = [], []
        for r in ratios:
            xs, ys = [], []
            for mesh in meshes:
                rec = index.get((mesh, contrast, r))
                if rec is not None and rec[field] > 0:
                    xs.append(node_counts[mesh])
                    ys.append(rec[field])
            if len(xs) >= 2:
                a, b = fit_power_law(xs, ys)
                coeffs.append(format(a, ".3g"))
                exps.append(format(b, ".3g"))
            else:
                coeffs.append("")
                exps.append("")
        lines.append("\\midrule")
        lines.append("fit $a$ & & " + " & ".join(coeffs) + " \\\\")
        lines.append("fit $b$ & & " + " & ".join(exps) + " \\\\")
    lines.append("\\bottomrule")
    return "\n".join(lines)


def main():
    records = read_logs(LOGDIR)
    print(f"{len(records)} logs read from {LOGDIR}.", file=sys.stderr)
    if not records:
        return

    # optional second set of logs, appended per cell as "value/value0"
    index0 = None
    if os.path.isdir(LOGDIR0):
        records0 = read_logs(LOGDIR0)
        print(f"{len(records0)} logs read from {LOGDIR0}.", file=sys.stderr)
        index0 = build_index(records0)

    # node count per mesh (drop meshes we cannot size), rows ordered by size
    node_counts = {}
    for mesh in sorted({rec["mesh"] for rec in records}):
        nn = node_count(mesh)
        if nn is None:
            print(f"warning: no '3D-Nodes' line in {mesh}, skipped.", file=sys.stderr)
            continue
        node_counts[mesh] = nn
    meshes = sorted(node_counts, key=node_counts.get)

    contrasts = sorted({rec["contrast"] for rec in records})
    ratios = sorted({rec["ratio"] for rec in records})
    index = build_index(records)

    blocks = []
    for contrast in contrasts:
        blocks.append(f"% solver iterations, contrast = {contrast}")
        blocks.append(make_table(index, meshes, node_counts, contrast, ratios, "iterations", "", index0))
        blocks.append(f"% total timing [s] (rel. to smallest mesh; fit T = a N^b), contrast = {contrast}")
        blocks.append(make_table(index, meshes, node_counts, contrast, ratios, "timing_total", ".3g", relative=True))

    with open(OUTFILE, "w") as f:
        f.write("\n".join(blocks) + "\n")
    print(f"tables written to {OUTFILE}.", file=sys.stderr)


if __name__ == "__main__":
    main()
