#!/usr/bin/python3
"""
Plot number of cores vs. solver timing for the CoalMine SIP forward runs.

Each file in LOGDIR is named "log<i>" where i is the number of cores used.
The last line of each file is a JSON run record holding "timing_solver" [s].
"""
import os
import re
import sys
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

LOGDIR = "logs/"
OUTFILE = "scaling.png"


def read_last_record(path):
    """Return the last non-empty JSON line of `path` as a dict, or None."""
    with open(path) as f:
        lines = [ln for ln in f.read().splitlines() if ln.strip()]
    if not lines:
        return None
    return json.loads(lines[-1])


def main():
    cores, timings = [], []
    for name in os.listdir(LOGDIR):
        m = re.fullmatch(r"log(\d+)", name)
        if not m:
            continue
        rec = read_last_record(os.path.join(LOGDIR, name))
        if rec is None or "timing_solver" not in rec:
            print(f"warning: no timing_solver in {name}, skipped.", file=sys.stderr)
            continue
        cores.append(int(m.group(1)))
        timings.append(rec["timing_solver"])

    if not cores:
        print("no usable logs found.", file=sys.stderr)
        return

    # sort by number of cores
    cores, timings = zip(*sorted(zip(cores, timings)))

    # least-squares power-law fit  T = a * n**b  in log-log space
    b, loga = np.polyfit(np.log(cores), np.log(timings), 1)
    a = np.exp(loga)

    fig, ax = plt.subplots()
    ax.loglog(cores, timings, "o", label="measured")
    fit = a * np.power(cores, b)
    ax.loglog(cores, fit, "-", label=f"fit: $T = {a:.3g}\\,n^{{{b:.3g}}}$")
    ax.legend()
    ax.set_xlabel("number of cores")
    ax.set_ylabel("solver time [s]")
    ax.set_title("CoalMine SIP forward solver scaling")
    ax.set_xticks(cores)
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.grid(True, which="both", ls=":")
    fig.tight_layout()
    fig.savefig(OUTFILE, dpi=150)
    print(f"plot written to {OUTFILE}.", file=sys.stderr)


if __name__ == "__main__":
    main()
