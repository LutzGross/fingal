"""
Post-process a TestDamageModel results directory into plots.

The test drivers stream their plot data to CSV files as the run proceeds
(history.csv, and slice.csv for the multi-element cases), so a plot can be
produced from whatever was written even if the run was terminated early. This
script reads those CSVs and writes the PNGs; it needs only numpy + matplotlib
(no esys-escript), so post-processing is cheap and decoupled from the run.

Usage:  python3 plot.py <results_dir> [<results_dir> ...]
        (default results_dir: the current directory)
"""
import csv
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


def read_csv(path):
    """read a headed CSV into {column_name: np.array(float)}."""
    with open(path) as f:
        rows = list(csv.DictReader(f))
    if not rows:
        return {}
    keys = rows[0].keys()
    return {k: np.array([float(r[k]) for r in rows]) for k in keys}


def plot_response(d, H):
    """stress-strain (and, if present, damage-evolution) from history.csv."""
    if "eqstrain" in H and "damage" in H:                # single-element: 2 panels
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        ax[0].plot(H["strain"], H["stress"] / 1e6, "b.-")
        ax[0].set(xlabel="axial compressive strain",
                  ylabel="axial compressive stress [MPa]", title="stress-strain")
        ax[0].grid(True, alpha=0.3)
        ax[1].plot(H["eqstrain"], H["damage"], "r.-")
        ax[1].set(xlabel=r"equivalent strain  $\tilde{\varepsilon}$",
                  ylabel="damage  $D$", title="damage evolution")
        ax[1].grid(True, alpha=0.3)
    else:
        plt.figure(figsize=(6, 4))
        plt.plot(H["strain"] * 1e3, H["stress"] / 1e6, "b.-")
        plt.xlabel("axial compressive strain  [1e-3]")
        plt.ylabel("axial compressive stress [MPa]")
        plt.title("stress-strain")
        plt.grid(True, alpha=0.3)
    plt.tight_layout()
    out = os.path.join(d, "response.png")
    plt.savefig(out, dpi=120)
    plt.close()
    print("wrote", out)

    # damage evolution vs strain, if peak/mean damage were recorded
    if "Dmax" in H or "Dmean" in H:
        plt.figure(figsize=(6, 4))
        if "Dmax" in H:
            plt.plot(H["strain"] * 1e3, H["Dmax"], "r.-", label="peak $D$")
        if "Dmean" in H:
            plt.plot(H["strain"] * 1e3, H["Dmean"], "b.-", label="mean $D$")
        plt.xlabel("axial compressive strain  [1e-3]")
        plt.ylabel("damage  $D$")
        plt.title("damage evolution")
        plt.ylim(0., 1.)
        plt.grid(True, alpha=0.3)
        plt.legend()
        plt.tight_layout()
        out = os.path.join(d, "damage_evolution.png")
        plt.savefig(out, dpi=120)
        plt.close()
        print("wrote", out)


def plot_slice(d, S):
    """damage on the structured mid-plane slice from slice.csv (imshow grid)."""
    xs = np.unique(S["x_mm"])
    zs = np.unique(S["z_mm"])
    grid = np.full((zs.size, xs.size), np.nan)
    ix = {v: i for i, v in enumerate(xs)}
    iz = {v: i for i, v in enumerate(zs)}
    for xv, zv, Dv in zip(S["x_mm"], S["z_mm"], S["D"]):
        grid[iz[zv], ix[xv]] = Dv
    dx = (xs[-1] - xs[0]) / max(xs.size - 1, 1)
    dz = (zs[-1] - zs[0]) / max(zs.size - 1, 1)
    # autoscale the colour range so a mild (pre-peak) concentration shows up as
    # well as a fully-developed band.
    plt.figure(figsize=(4.2, 6))
    im = plt.imshow(grid, origin="lower", cmap="inferno",
                    vmin=np.nanmin(grid), vmax=np.nanmax(grid),
                    aspect="equal", extent=[xs[0] - dx / 2, xs[-1] + dx / 2,
                                            zs[0] - dz / 2, zs[-1] + dz / 2])
    plt.colorbar(im, label="damage  $D$")
    plt.xlabel("x [mm]")
    plt.ylabel("z [mm]")
    plt.title("damage slice")
    plt.tight_layout()
    out = os.path.join(d, "damage_slice.png")
    plt.savefig(out, dpi=120)
    plt.close()
    print("wrote", out)


def process(d):
    hist = os.path.join(d, "history.csv")
    if os.path.exists(hist):
        H = read_csv(hist)
        if H:
            plot_response(d, H)
    sl = os.path.join(d, "slice.csv")
    if os.path.exists(sl):
        S = read_csv(sl)
        if S:
            plot_slice(d, S)


if __name__ == "__main__":
    dirs = sys.argv[1:] or ["."]
    for d in dirs:
        process(d)
