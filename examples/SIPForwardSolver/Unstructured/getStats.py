#!/usr/bin/python3
"""
Aggregate SIP forward-solver run logs into LaTeX tables.

Each file in LOGDIR holds one JSON run record per line; the last line is
taken as the record for that log. For every contrast (excluding the first,
baseline value) two tables are written to OUTFILE: solver iterations and
total timing, with mesh node counts as rows and the anomaly `ratio` as
columns. Progress and warnings go to stderr.
"""
import os
import sys
import json

LOGDIR = "logs/"
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


def node_count(flyfile):
    """Number of 3D nodes declared in a finley .fly mesh file, or None."""
    with open(flyfile) as f:
        for row in f:
            if row.startswith("3D-Nodes"):
                return int(row[len("3D-Nodes"):])
    return None


def make_table(index, meshes, node_counts, contrast, ratios, field, fmt):
    """Build one LaTeX table for `field` at `contrast`.

    `index` maps (mesh, contrast, ratio) -> record. Missing cells are blank.
    """
    header = " # nodes" + "".join(f" & {r}" for r in ratios) + "\\\\"
    lines = ["\\toprule", header, "\\midrule"]
    for mesh in meshes:
        cells = []
        for r in ratios:
            rec = index.get((mesh, contrast, r))
            cells.append(format(rec[field], fmt) if rec is not None else "")
        lines.append(f"{node_counts[mesh]} " + "".join(f" & {c}" for c in cells) + " \\\\")
    lines.append("\\bottomrule")
    return "\n".join(lines)


def main():
    records = read_logs(LOGDIR)
    print(f"{len(records)} logs read.", file=sys.stderr)
    if not records:
        return

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

    # (mesh, contrast, ratio) -> record; first record wins on duplicates
    index = {}
    for rec in records:
        index.setdefault((rec["mesh"], rec["contrast"], rec["ratio"]), rec)

    blocks = []
    for contrast in contrasts[1:]:
        blocks.append(f"% solver iterations, contrast = {contrast}")
        blocks.append(make_table(index, meshes, node_counts, contrast, ratios, "iterations", ""))
        blocks.append(f"% total timing [s], contrast = {contrast}")
        blocks.append(make_table(index, meshes, node_counts, contrast, ratios, "timing_total", ".3g"))

    with open(OUTFILE, "w") as f:
        f.write("\n".join(blocks) + "\n")
    print(f"tables written to {OUTFILE}.", file=sys.stderr)


if __name__ == "__main__":
    main()
