"""Radial ionization and temperature structure of the Lexington HII40 and
HII20 benchmarks, in the panel layout of Wood, Mathis & Ercolano (2004,
MNRAS 348, 1337), their Figures 1 and 2.

Wood et al. plot, against radius, the neutral fractions H0/H and He0/He, the
ion fractions of O, C, N and Ne, and the electron temperature, one point for
each cell of their Cartesian grid, and overlay the Cloudy solution of the
same model as a line.  The panels here carry the same quantities on the same
axes for MoCHII (points) and Cloudy c25.00 (black lines), so that the two
can be laid side by side.  Nitrogen stops at N++ because the MoCHII registry
carries three nitrogen stages for these models; Wood et al. also show
N+3/N.  Following their convention, the temperature panel shows only cells
that are less than 95 per cent neutral, since the temperature of an almost
neutral cell is not a temperature of the H II region; the same cut is
applied to the Cloudy zones.

Data (do not change without re-running the benchmark):

  tests/g2_hii/hii40_bench_rates.h5   from tests/g2_hii/hii40_bench.in
  tests/g2_hii/hii20_bench_rates.h5   from tests/g2_hii/hii20_bench.in

Both are 4e7-packet, level-6 octree (262144 leaves) runs with the thermal
balance and the five metals, the explicit Monte Carlo diffuse field
including the He I channels (case A), and the
threshold-aligned 32-bin ionizing grid (par%ion_align_edges, the default),
converged on the volume-integrated criterion.  These are the SAME two runs
that tables/make_wood_tables.py reads for the MoCHII column of the two
benchmark tables, so the tables and these figures describe one calculation.

The point clouds are a random thinning of the full leaf list for legibility;
every number quoted in the text is computed from all leaves.

The Cloudy c25.00 lines come from the runs that also feed the C25 column of
the two benchmark tables (paper/tables/make_wood_tables.py),

  tests/cloudy_c25/hii40_c25.in    tests/cloudy_c25/hii20_c25.in

whose physics blocks are the unmodified c25.00 test-suite decks
tsuite/auto/hii_paris.in and tsuite/auto/hii_coolstar.in.  The radial
structure is read from the save files those inputs write,

  hii{40,20}_c25.ovr        depth [cm] and electron temperature [K]
  hii{40,20}_c25.ele_h      H0/H, H+/H
  hii{40,20}_c25.ele_he     He0/He, He+/He, He+2/He
  hii{40,20}_c25.ele_c      C0/C ... C+6/C
  hii{40,20}_c25.ele_n      N0/N ... N+7/N
  hii{40,20}_c25.ele_o      O0/O ... O+8/O
  hii{40,20}_c25.ele_ne     Ne0/Ne ... Ne+10/Ne

Cloudy tabulates the depth to the center of each zone, so the radius is the
inner radius of the model (10^18.477121 = 3.0e18 cm, the `radius' command of
both decks) plus that depth; every save file of one run carries the same
depth column, which is checked on read.  The Cloudy fractions are fractions
of the element, the same definition MoCHII writes, and the curves are cut at
the outer radius of the MoCHII shell so that both codes cover one volume.

Output: figures/wood_hii40.pdf, figures/wood_hii20.pdf

Run:  python3 make_wood_benchmark.py
"""

import os
import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "figures")
TESTS = os.path.join(HERE, "..", "..", "tests", "g2_hii")
CLOUDY = os.path.join(HERE, "..", "..", "tests", "cloudy_c25")

# shell of the benchmark grids (tests/g2_hii/make_hii_grid.py)
RIN = 0.97223                      # pc, = 3e18 cm
ROUT = {"hii40": 4.73154, "hii20": 3.07874}

PC_CM = 3.0856776e18               # pc2cm of src/define.f90
RIN_CM = 3.0e18                    # Cloudy `radius = 18.477121'

THIN = 6000                        # points drawn for each curve
RNG = np.random.default_rng(20260725)

plt.rcParams.update({"font.size": 10, "axes.linewidth": 0.9})
# Paper figures: axis, tick, and title fonts sized to ~85% of the body text.
_S = 1.31
_b = plt.rcParams["font.size"]
plt.rcParams.update({
    "axes.labelsize":   _b * _S,
    "axes.titlesize":   _b * 1.2 * _S,
    "xtick.labelsize":  _b * _S,
    "ytick.labelsize":  _b * _S,
    "figure.titlesize": _b * 1.2 * _S,
})

BLUE, RED, GREEN = "#1f77b4", "#d62728", "#2ca02c"
ORANGE = "#c1440e"
CLOUDY_KW = dict(color="black", lw=1.1, zorder=5)
CLOUDY_LABEL = r'\textsc{Cloudy}~c25.00'


def read_run(tag):
    with h5py.File(os.path.join(TESTS, tag + "_bench_rates.h5"), "r") as f:
        d = {k: f[k + "/data"][:] for k in
             ("LeafXYZ", "T_e", "n_e", "x_HI", "x_HeI", "x_HeII")}
        for el in ("c", "n", "o", "ne", "s"):
            a = f["x_%s_stages/data" % el][:]
            d["x_" + el] = a if a.shape[0] == d["T_e"].size else a.T
    r = np.sqrt((d["LeafXYZ"] ** 2).sum(axis=0))
    gas = (r >= RIN) & (r < ROUT[tag])
    d = {k: (v[:, gas] if k == "LeafXYZ" else v[gas]) for k, v in d.items()}
    d["r"] = r[gas]
    return d


def read_cloudy_save(path):
    """Columns of one Cloudy save file, keyed by its tab-separated header."""
    with open(path) as fh:
        head = fh.readline().lstrip("#").rstrip("\n").split("\t")
    a = np.loadtxt(path, skiprows=1)
    return {name.strip(): a[:, i] for i, name in enumerate(head)}


def stage_column(symbol, stage):
    """Cloudy column name of an ionization stage: H, H+, He+2, C+3, ..."""
    if stage == 0:
        return symbol
    return symbol + "+" + ("" if stage == 1 else str(stage))


def read_cloudy(tag):
    """Radial structure of the Cloudy c25.00 model of the same benchmark.

    The radius is the inner radius of the model plus the depth to the zone
    center that Cloudy tabulates; the same depth column appears in every
    save file of one run, which is asserted here.  Stage fractions are
    fractions of the element, as MoCHII writes them.  The profile is cut at
    the outer radius of the MoCHII shell.
    """
    ovr = read_cloudy_save(os.path.join(CLOUDY, tag + "_c25.ovr"))
    keep = (RIN_CM + ovr["depth"]) / PC_CM <= ROUT[tag]
    d = {"r": ((RIN_CM + ovr["depth"]) / PC_CM)[keep],
         "T_e": ovr["Te"][keep]}
    for el, symbol, nstage in (("h", "H", 2), ("he", "He", 3), ("c", "C", 4),
                               ("n", "N", 3), ("o", "O", 3), ("ne", "Ne", 3)):
        t = read_cloudy_save(os.path.join(CLOUDY,
                                          "%s_c25.ele_%s" % (tag, el)))
        assert np.array_equal(t["depth"], ovr["depth"])
        d["x_" + el] = np.column_stack(
            [t[stage_column(symbol, i)][keep] for i in range(nstage)])
    d["x_HI"] = d["x_h"][:, 0]
    d["x_HeI"] = d["x_he"][:, 0]
    return d


def cloud(ax, r, q, color, label, floor):
    """Thinned point cloud of a leaf quantity against radius."""
    m = np.isfinite(q) & (q > floor)
    idx = np.flatnonzero(m)
    if idx.size > THIN:
        idx = RNG.choice(idx, THIN, replace=False)
    ax.plot(r[idx], q[idx], ".", ms=1.6, color=color, alpha=0.35,
            rasterized=True)
    ax.plot([], [], "s", ms=5, color=color, label=label)   # legend proxy


def cloudy_curve(ax, r, q, floor):
    """Cloudy profile of the same quantity, as a line over the point cloud.

    Zones below the floor are blanked rather than dropped, so that a gap in
    the profile is not bridged by a straight segment.
    """
    ax.plot(r, np.where(np.isfinite(q) & (q > floor), q, np.nan),
            "-", **CLOUDY_KW)


def ion_panel(ax, d, c, series, ylim, ylabel, loc):
    floor = 0.2 * ylim[0]
    for (q, cq, color, label) in series:
        if not (np.any(q > ylim[0]) or np.any(cq > ylim[0])):
            continue                 # no legend entry for an absent stage
        cloud(ax, d["r"], q, color, label, floor)
        cloudy_curve(ax, c["r"], cq, floor)
    ax.plot([], [], "-", label=CLOUDY_LABEL, **CLOUDY_KW)   # legend proxy
    ax.set_yscale("log")
    ax.set_ylim(*ylim)
    ax.set_ylabel(ylabel)
    ax.grid(True, which="major", ls=":", lw=0.5, alpha=0.5)
    leg = ax.legend(loc=loc, fontsize=_b * _S * 0.85, handletextpad=0.3,
                    borderaxespad=0.4, labelspacing=0.25, framealpha=0.88,
                    edgecolor="none")
    leg.get_frame().set_facecolor("white")


def make_figure(tag, ylim_neutral, ylim_temp, loc_n, loc_i, loc_t):
    d = read_run(tag)
    c = read_cloudy(tag)
    ylim_ion = (1e-3, 2.0)
    fig, axes = plt.subplots(2, 3, figsize=(11.0, 6.6), sharex=True)
    for ax in axes[1]:
        ax.set_xlabel(r'$r$~[pc]')
    for ax in axes.ravel():
        ax.set_xlim(0.9, ROUT[tag] * 1.02)

    ion_panel(axes[0, 0], d, c,
              [(d["x_HI"], c["x_HI"], BLUE, r'H$^{0}$/H'),
               (d["x_HeI"], c["x_HeI"], RED, r'He$^{0}$/He')],
              ylim_neutral, r'neutral fraction', loc_n)
    ion_panel(axes[0, 1], d, c,
              [(d["x_o"][:, 1], c["x_o"][:, 1], BLUE, r'O$^{+}$/O'),
               (d["x_o"][:, 2], c["x_o"][:, 2], RED, r'O$^{+2}$/O')],
              ylim_ion, r'ion fraction', loc_i)
    ion_panel(axes[0, 2], d, c,
              [(d["x_c"][:, 1], c["x_c"][:, 1], BLUE, r'C$^{+}$/C'),
               (d["x_c"][:, 2], c["x_c"][:, 2], RED, r'C$^{+2}$/C'),
               (d["x_c"][:, 3], c["x_c"][:, 3], GREEN, r'C$^{+3}$/C')],
              ylim_ion, r'ion fraction', loc_i)
    ion_panel(axes[1, 0], d, c,
              [(d["x_n"][:, 1], c["x_n"][:, 1], BLUE, r'N$^{+}$/N'),
               (d["x_n"][:, 2], c["x_n"][:, 2], RED, r'N$^{+2}$/N')],
              ylim_ion, r'ion fraction', loc_i)
    ion_panel(axes[1, 1], d, c,
              [(d["x_ne"][:, 1], c["x_ne"][:, 1], BLUE, r'Ne$^{+}$/Ne'),
               (d["x_ne"][:, 2], c["x_ne"][:, 2], RED, r'Ne$^{+2}$/Ne')],
              ylim_ion, r'ion fraction', loc_i)

    ax = axes[1, 2]
    ionized = d["x_HI"] < 0.95
    ri, ti = d["r"][ionized], d["T_e"][ionized]
    idx = RNG.choice(ri.size, min(THIN, ri.size), replace=False)
    ax.plot(ri[idx], ti[idx], ".", ms=1.6, color=ORANGE, alpha=0.35,
            rasterized=True)
    ion = c["x_HI"] < 0.95           # the same convention on the Cloudy zones
    ax.plot(c["r"][ion], c["T_e"][ion], "-", **CLOUDY_KW)
    ax.plot([], [], "s", ms=5, color=ORANGE, label=r'$T_{\rm e}$')
    ax.plot([], [], "-", label=CLOUDY_LABEL, **CLOUDY_KW)
    ax.set_ylim(*ylim_temp)
    ax.set_ylabel(r'$T_{\rm e}$~[K]')
    ax.grid(True, which="major", ls=":", lw=0.5, alpha=0.5)
    leg = ax.legend(loc=loc_t, fontsize=_b * _S * 0.85, handletextpad=0.3,
                    borderaxespad=0.4, labelspacing=0.25, framealpha=0.88,
                    edgecolor="none")
    leg.get_frame().set_facecolor("white")

    fig.tight_layout()
    out = os.path.join(FIGDIR, "wood_%s.pdf" % tag)
    fig.savefig(out, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print("wrote %s   nleaf(gas) = %d,  T_e(x_HI<0.95) = %.0f-%.0f K"
          "   Cloudy: %d zones, r = %.3f-%.3f pc, T_e(x_HI<0.95) = %.0f-%.0f K"
          % (out, d["r"].size, ti.min(), ti.max(), c["r"].size,
             c["r"].min(), c["r"].max(),
             c["T_e"][ion].min(), c["T_e"][ion].max()))


# The figure title lives in the paper caption, so none is drawn here.
# The temperature limits bracket both plotted populations of each run, so
# nothing that passes the x_HI < 0.95 cut is clipped out of the panel: neither
# a MoCHII cell nor a Cloudy zone.  Both floors are set by the front.  At HII40
# the Cloudy profile reaches 7331 K in the middle of the nebula, well below the
# coolest MoCHII cell.  At HII20 the MoCHII cells just outside the front reach
# the te_min floor of the run (3000 K, 0.6% of the cells that pass the cut) and
# the Cloudy profile plunges to 2919 K over the last 4e-4 pc before its own
# front, where its zoning is far finer than a MoCHII cell.
make_figure("hii40", (1e-5, 2.0), (7200, 10800), loc_n="lower right",
            loc_i="lower left", loc_t="upper left")
make_figure("hii20", (1e-5, 2.0), (2800, 9200), loc_n="center left",
            loc_i="center", loc_t="lower left")
