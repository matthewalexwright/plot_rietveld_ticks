r"""
Rietveld refinement plot (TOPAS output) -- command-line friendly.

The script can live anywhere (e.g. a PyScripts folder). It looks in the
target folder for TOPAS output:

    <stem>.txt            x, Yobs, [SigmaYobs,] Ycalc, Diff  (comma-separated)
    *_2Th_Ip.txt          hkl tick positions - one file per phase

Usage examples (Windows shown; same on Mac/Linux with forward slashes):

    cd "C:\XRD data\NaRhO2 run 3"
    python C:\PyScripts\plot_rietveld_ticks.py            -> finds data here, asks 2Th or Q
    python C:\PyScripts\plot_rietveld_ticks.py Q          -> no prompt, Q axis

    python C:\PyScripts\plot_rietveld_ticks.py "C:\XRD data\NaRhO2 run 3"      (folder as argument)
    python C:\PyScripts\plot_rietveld_ticks.py "C:\XRD data\NaRhO2 run 3" Q
    python C:\PyScripts\plot_rietveld_ticks.py Stoich Q     (part of a file name picks the dataset)

If a folder holds several refinements you get a numbered menu; typing part
of a file name on the command line (as above) skips the menu. If several
*_2Th_Ip.txt files belong to the dataset they are treated as separate
phases: the script lists them and asks for the phase order, then for a
legend label for each. Every phase gets its own colour and its own row of
tick marks.

Options:  --xlim MIN MAX   --ylim MIN MAX   --wavelength 0.7093
          --label "NaRhO2" ["Rh2O3" ...]    (one label per phase, in order)
          --size 3 3       (axes size in inches; one number = square)
          --data FILE --ticks FILE [FILE ...]      (pick files explicitly)

Prompts (Enter accepts the bracketed default): the axis, the wavelength when
Q is chosen, the phase order and labels, and the axes size in inches. Typing
a phase name at a label prompt sets digits as subscripts automatically,
including refined compositions with uncertainties (NaRhO2, Fe2O3,
Na0.96(1)RhO2).

Figures are saved next to the data as <data stem>_2Th.png/.pdf or _Q.png/.pdf.
Also still runs in Jupyter: prompts appear as input boxes, or preset the
DEFAULT_* values below.
"""
import argparse
import glob
import os
import re
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator, MultipleLocator
from mpl_toolkits.axes_grid1 import Divider, Size

# ---------- defaults (handy to edit when running in Jupyter) ----------
DEFAULT_AXIS       = None    # "2theta" or "Q"; None = ask
DEFAULT_XLIM       = None    # e.g. (30, 65); units follow the axis choice
DEFAULT_YLIM       = None    # e.g. (-15, 42), in 10^3 counts
DEFAULT_WAVELENGTH = None    # angstrom; None = ask when plotting in Q
DEFAULT_LABELS     = None    # e.g. ["NaRhO2", "Rh2O3"]; None = ask
DEFAULT_SIZE       = None    # axes size in inches, e.g. (3, 3) or 3; None = ask
AXES_SIZE_DEFAULT  = (3.0, 3.0)   # offered at the prompt
CU_KA1 = 1.5405980           # angstrom; Cu Ka1 (Holzer et al. 1997) - offered as the default
PLOT_ONLY_FIT_RANGE = True   # Ycalc = 0 outside the refined window; mask those points
SCALE        = 1e-3          # intensities handled in 10^3 counts
TICK_SUFFIX  = "_2Th_Ip.txt"

# vertical layout, in typographic points (independent of the intensity scale,
# so the spacing looks the same whatever the counts happen to be)
MAJOR_TICK_PT = 10  # x-axis major ticks, drawn inward from the top and bottom spines
PAD_TOP_PT = MAJOR_TICK_PT + 4   # above the tallest peak: clears the top-spine ticks
PAD_BOT_PT = MAJOR_TICK_PT + 8   # below the difference curve: clears the bottom ticks
GAP1_PT    = 5      # pattern baseline -> top of the first row of hkl ticks
ROW_PT_GAP = 2      # between one row of ticks and the next
GAP2_PT    = 4      # bottom of the last tick row -> top of the difference curve
TICK_MS    = 7      # hkl tick mark height, points
# PAD_TOP/PAD_BOT are held fixed so the pattern and the difference curve always
# clear the inward-pointing axis ticks; only the gaps between them are squeezed
# if the axes are too short to hold everything at full size.
# ----------------------------------------------------------------------


def _norm(stem):
    """Filename stem normalized for matching: casefolded, repeated _ collapsed."""
    return re.sub(r"_+", "_", stem).strip("_").lower()


def tick_files(folder):
    return sorted(glob.glob(os.path.join(folder, "*" + TICK_SUFFIX)))


def data_files(folder):
    return sorted(f for f in glob.glob(os.path.join(folder, "*.txt"))
                  if not f.endswith(TICK_SUFFIX))


def ticks_for(data_file, ticks):
    """Tick files belonging to *data_file*, matched on how many leading
    underscore-separated name parts they share with it. So EXAMPLE_lab_data
    claims EXAMPLE_lab_PhaseA/B/C (two parts shared) in preference to
    EXAMPLE_synchrotron_PhaseA (one), and <stem>_2Th_Ip.txt always wins."""
    dt = _norm(os.path.splitext(os.path.basename(data_file))[0]).split("_")
    best, hits = 0, []
    for t in ticks:
        tt = _norm(os.path.basename(t)[: -len(TICK_SUFFIX)]).split("_")
        n = 0
        while n < min(len(dt), len(tt)) and dt[n] == tt[n]:
            n += 1
        if n > best:
            best, hits = n, [t]
        elif n == best and n > 0:
            hits.append(t)
    return hits if best else []


def pretty_label(text):
    """Set digits as chemical-formula subscripts: NaRhO2 -> NaRhO$_2$,
    Na0.75CoO2 -> Na$_{0.75}$CoO$_2$. A parenthesized uncertainty rides
    along: Na0.96(1)RhO2 -> Na$_{0.96(1)}$RhO$_2$. Text already
    containing mathtext ($...$) is left untouched."""
    if "$" in text:
        return text
    return re.sub(r"(?<=[A-Za-z)])(\d+(?:\.\d+)?(?:\(\d+\))?)", r"$_{\1}$", text)


def ask_int_list(prompt, n_max, default):
    """Read a list of 1-based indices, e.g. '2 1 3'. Enter accepts *default*."""
    while True:
        s = input(prompt).strip()
        if s == "":
            return default
        try:
            picked = [int(v) for v in s.replace(",", " ").split()]
        except ValueError:
            print("  enter numbers from the list, e.g. 2 1 3")
            continue
        if any(not 1 <= v <= n_max for v in picked) or len(set(picked)) != len(picked):
            print(f"  use each number once, from 1 to {n_max}")
            continue
        return picked


# ---------------- command line (skipped inside Jupyter) ----------------
VERSION = "2026-08-22"
IN_JUPYTER = "ipykernel" in sys.modules
parser = argparse.ArgumentParser(
    description="Plot a TOPAS Rietveld refinement (see the file header for examples).")
parser.add_argument("where", nargs="*",
                    help="optionally: a data folder, an axis keyword (2Th or Q), and/or "
                         "part of a file name to pick a dataset - any order")
parser.add_argument("--data", help="main data file (x, Yobs, [Sigma,] Ycalc, Diff)")
parser.add_argument("--ticks", nargs="+", help="hkl tick file(s), in phase order")
parser.add_argument("--wavelength", type=float, default=None,
                    help="wavelength in angstrom for the Q axis (skips the prompt; "
                         "prompt default is Cu Ka1)")
parser.add_argument("--size", nargs="+", type=float, metavar=("W", "H"), default=None,
                    help="axes size in inches: one number for square, or width height "
                         "(skips the prompt; default 3 3)")
parser.add_argument("--label", nargs="+", default=None,
                    help='legend label per phase, in order, e.g. "NaRhO2" "Rh2O3"')
parser.add_argument("--xlim", nargs=2, type=float, metavar=("MIN", "MAX"), default=DEFAULT_XLIM)
parser.add_argument("--ylim", nargs=2, type=float, metavar=("MIN", "MAX"), default=DEFAULT_YLIM)
args = parser.parse_args([] if IN_JUPYTER else None)

X_AXIS, folder, name_filters = DEFAULT_AXIS, None, []
for tok in args.where:
    key = tok.lower()
    if key in ("2th", "2theta", "tth"):
        X_AXIS = "2theta"
    elif key == "q":
        X_AXIS = "Q"
    elif os.path.isdir(tok):
        folder = tok
    else:
        name_filters.append(tok)        # part of a file name, used to pick a dataset
folder = folder or os.getcwd()

# ---------------- pick the data file ----------------
if args.data:
    DATA_FILE = args.data
else:
    cands = data_files(folder)
    if not cands:
        sys.exit(f"No data .txt files found in {folder}\n"
                 "(files can also be named explicitly with --data and --ticks)")
    paired = [d for d in cands if ticks_for(d, tick_files(folder))]
    cands = paired or cands
    if name_filters:
        matches = [d for d in cands
                   if all(s.lower() in os.path.basename(d).lower() for s in name_filters)]
        if not matches:
            sys.exit("No dataset matches " + " ".join(repr(s) for s in name_filters)
                     + "\nAvailable:\n  " + "\n  ".join(os.path.basename(d) for d in cands))
        cands = matches
    if len(cands) == 1:
        DATA_FILE = cands[0]
    else:
        print("Found these datasets:")
        for i, d in enumerate(cands, 1):
            print(f"  {i}) {os.path.basename(d)}")
        while True:
            s = input(f"Which one? [1-{len(cands)}]: ").strip()
            if s.isdigit() and 1 <= int(s) <= len(cands):
                break
            print("  please enter a number from the list")
        DATA_FILE = cands[int(s) - 1]
print("plot_rietveld_ticks", VERSION)
print("Data :", os.path.basename(DATA_FILE))

# ---------------- pick the phases (hkl tick files) ----------------
if args.ticks:
    TICK_LIST = list(args.ticks)
else:
    found = ticks_for(DATA_FILE, tick_files(folder)) or tick_files(folder)
    if not found:
        sys.exit(f"No '*{TICK_SUFFIX}' file found in {folder}")
    if len(found) == 1:
        TICK_LIST = found
    else:
        print(f"Found {len(found)} hkl tick files:")
        for i, t in enumerate(found, 1):
            print(f"  {i}) {os.path.basename(t)}")
        order = ask_int_list(
            "Phase order - main phase first, e.g. 2 1 3 "
            f"[{' '.join(str(i) for i in range(1, len(found) + 1))}]: ",
            len(found), list(range(1, len(found) + 1)))
        TICK_LIST = [found[i - 1] for i in order]
for i, t in enumerate(TICK_LIST, 1):
    print(f"Phase {i}: {os.path.basename(t)}")
N_PHASES = len(TICK_LIST)

# ---------------- axis choice ----------------
while X_AXIS is None:
    s = input("x-axis - 2Th or Q? [2Th]: ").strip().lower()
    if s in ("", "2th", "2theta", "tth"):
        X_AXIS = "2theta"
    elif s == "q":
        X_AXIS = "Q"
    else:
        print("  please answer 2Th or Q")

# ---------------- wavelength: confirm when plotting in Q ----------------
WAVELENGTH = args.wavelength if args.wavelength is not None else DEFAULT_WAVELENGTH
if X_AXIS == "Q":
    while WAVELENGTH is None or WAVELENGTH <= 0:
        s = input(f"wavelength in angstrom [{CU_KA1} = Cu Ka1]: ").strip()
        if s == "":
            WAVELENGTH = CU_KA1
        else:
            try:
                WAVELENGTH = float(s)
            except ValueError:
                print("  please enter a number, e.g. 0.7093")
        if WAVELENGTH is not None and WAVELENGTH <= 0:
            print("  wavelength must be positive")
            WAVELENGTH = None

# ---------------- legend label for each phase ----------------
supplied = args.label if args.label is not None else DEFAULT_LABELS
LABELS = []
for i in range(N_PHASES):
    fallback = r"$hkl$" if N_PHASES == 1 else f"phase {i + 1}"
    if supplied is not None and i < len(supplied):
        LABELS.append(pretty_label(supplied[i]))
        continue
    shown = "hkl" if N_PHASES == 1 else f"phase {i + 1}"
    tag = "" if N_PHASES == 1 else f" ({os.path.basename(TICK_LIST[i])})"
    s = input(f"legend label for phase {i + 1}{tag} [{shown}]: ").strip()
    LABELS.append(fallback if s == "" else pretty_label(s))

# ---------------- axes size in inches ----------------
def _parse_size(vals):
    """One number -> square; two -> width, height."""
    if len(vals) == 1:
        return (float(vals[0]), float(vals[0]))
    if len(vals) == 2:
        return (float(vals[0]), float(vals[1]))
    raise ValueError

AXES_SIZE = None
if args.size is not None:
    AXES_SIZE = _parse_size(args.size)
elif DEFAULT_SIZE is not None:
    AXES_SIZE = _parse_size(np.atleast_1d(DEFAULT_SIZE))
while AXES_SIZE is None or min(AXES_SIZE) <= 0:
    w0, h0 = AXES_SIZE_DEFAULT
    s = input(f"axes size in inches - one number, or width height [{w0:g} {h0:g}]: ").strip()
    if s == "":
        AXES_SIZE = AXES_SIZE_DEFAULT
    else:
        try:
            AXES_SIZE = _parse_size(s.replace(",", " ").split())
        except ValueError:
            print("  enter e.g. 3   or   6 3")
            continue
    if min(AXES_SIZE) <= 0:
        print("  sizes must be positive")
        AXES_SIZE = None

X_LIMITS = tuple(args.xlim) if args.xlim else None
Y_LIMITS = tuple(args.ylim) if args.ylim else None

# figure files take their name from the data file
STEM = os.path.splitext(os.path.basename(DATA_FILE))[0]

# ---------------- global style ----------------
plt.rcParams.update({
    "font.size": 15,
    "axes.linewidth": 1,
    "xtick.direction": "in",
    "xtick.top": True,
    "xtick.major.size": MAJOR_TICK_PT, "xtick.minor.size": MAJOR_TICK_PT / 2,
    "xtick.major.width": 1, "xtick.minor.width": 1,
})

# ---------------- data ----------------
# TOPAS writes a comment header naming the columns, e.g.
#   'x,<datafile>.xy,Ycalc,Diff                  (4 columns)
#   'x,<datafile>.xye,SigmaYobs,Ycalc,Diff       (5 columns, with uncertainties)
# so the layout is read from that line rather than assumed.
with open(DATA_FILE) as fh:
    header = [ln for ln in (fh.readline(), fh.readline()) if ln.startswith("'")]
cols = header[-1].lstrip("'").strip().split(",") if header else []

def _col(*names):
    for i, c in enumerate(cols):
        if c.strip().lower() in names:
            return i
    return None

i_calc, i_diff = _col("ycalc"), _col("diff")
i_sig = _col("sigmayobs", "sigma", "esd", "error")
if i_calc is None or i_diff is None:            # unrecognised header: fall back to position
    ncol = np.loadtxt(DATA_FILE, delimiter=",", comments="'", max_rows=1).size
    i_calc, i_diff = (3, 4) if ncol >= 5 else (2, 3)
    i_sig = 2 if ncol >= 5 else None
# Yobs is the remaining column (TOPAS names it after the source data file)
i_obs = next(i for i in range(1, max(i_calc, i_diff) + 1)
             if i not in (i_sig, i_calc, i_diff))

raw = np.loadtxt(DATA_FILE, delimiter=",", comments="'")
x, yobs, ycalc, ydiff = raw[:, 0], raw[:, i_obs], raw[:, i_calc], raw[:, i_diff]
print(f"Columns: {raw.shape[1]} (x=1, Yobs={i_obs+1}, Ycalc={i_calc+1}, Diff={i_diff+1})")

if PLOT_ONLY_FIT_RANGE:
    m = ycalc != 0
    x, yobs, ycalc, ydiff = x[m], yobs[m], ycalc[m], ydiff[m]
yobs, ycalc, ydiff = yobs * SCALE, ycalc * SCALE, ydiff * SCALE

hkl_list = [np.atleast_2d(np.loadtxt(t))[:, 0] for t in TICK_LIST]

if X_AXIS == "Q":                          # Q = 4 pi sin(theta) / lambda
    tth_to_q = lambda t: 4.0 * np.pi * np.sin(np.radians(t / 2.0)) / WAVELENGTH
    x = tth_to_q(x)
    hkl_list = [tth_to_q(h) for h in hkl_list]
    xlabel = r"$Q$ ($\mathrm{\AA}^{-1}$)"
else:
    xlabel = r"$2\theta\ (^{\circ})$"

# ---------------- vertical layout ----------------
# Spacing is set in points and converted to data units, so the gap between
# the pattern, the phase rows and the difference curve is always the same
# physical size no matter how large the intensities or the misfit are.
if X_LIMITS is not None:
    w = (x >= min(X_LIMITS)) & (x <= max(X_LIMITS))
    if not w.any():
        sys.exit("X_LIMITS lie outside the data range")
else:
    w = np.ones(x.size, dtype=bool)

ax_w, ax_h = AXES_SIZE
y_hi, y_lo = yobs[w].max(), yobs[w].min()
A = max(y_hi - y_lo, 1e-12)                        # pattern height, data units
d_up, d_dn = max(ydiff[w].max(), 0.0), min(ydiff[w].min(), 0.0)
D = d_up - d_dn                                    # difference curve height

axes_pt = 72.0 * ax_h                              # axes height in points
fixed_pt = PAD_TOP_PT + PAD_BOT_PT                 # must never be squeezed
flex_pt  = (GAP1_PT + GAP2_PT
            + N_PHASES * TICK_MS + (N_PHASES - 1) * ROW_PT_GAP)
f = 1.0                                            # squeeze factor for the flexible gaps
if fixed_pt + flex_pt > 0.85 * axes_pt:            # very short axes / many phases
    f = max((0.85 * axes_pt - fixed_pt) / flex_pt, 0.15)
g1, g2, tms = f * GAP1_PT, f * GAP2_PT, f * TICK_MS
rowsp = tms + f * ROW_PT_GAP        # centre-to-centre spacing of the phase rows

gap_pt = fixed_pt + f * flex_pt
T = (A + D) / (1.0 - gap_pt / axes_pt)             # total data range of the y axis
p = T / axes_pt                                    # data units per point

rows  = [y_lo - (g1 + 0.5 * tms + i * rowsp) * p for i in range(N_PHASES)]
d_off = rows[-1] - (0.5 * tms + g2) * p - d_up
y_top = y_hi + PAD_TOP_PT * p
y_bot = d_off + d_dn - PAD_BOT_PT * p

# ---------------- figure: axes exactly AXES_SIZE inches ----------------
pad_l, pad_r, pad_b, pad_t = 0.2, 0.3, 0.8, 0.2      # margins in inches
fig = plt.figure(figsize=(pad_l + ax_w + pad_r, pad_b + ax_h + pad_t))
h = [Size.Fixed(pad_l), Size.Fixed(ax_w), Size.Fixed(pad_r)]
v = [Size.Fixed(pad_b), Size.Fixed(ax_h), Size.Fixed(pad_t)]
div = Divider(fig, (0, 0, 1, 1), h, v, aspect=False)
ax = fig.add_axes(div.get_position(), axes_locator=div.new_locator(nx=1, ny=1))

# ---------------- colours: viridis, with a gray difference curve ----------------
cmap = plt.get_cmap("viridis")
c_obs, c_calc = cmap(0.0), cmap(0.55)
c_diff = "0.45"                                    # medium gray
if N_PHASES == 1:
    phase_colors = [cmap(0.8)]
else:
    # stops chosen to stay clear of the Ycalc teal (~0.55) and of each other
    stops = [0.90, 0.10, 0.72, 0.30, 0.98, 0.20, 0.80, 0.38]
    phase_colors = [cmap(stops[i % len(stops)]) for i in range(N_PHASES)]

ax.plot(x, yobs, "o", ms=2.2, mfc="none", mec=c_obs, mew=0.7,
        label=r"$Y_\mathrm{obs}$", zorder=3)
ax.plot(x, ycalc, "-", lw=1.2, color=c_calc, label=r"$Y_\mathrm{calc}$", zorder=4)
ax.plot(x, ydiff + d_off, "-", lw=0.9, color=c_diff,
        label=r"$Y_\mathrm{obs}-Y_\mathrm{calc}$", zorder=2)
for hkl, y_row, col, lab in zip(hkl_list, rows, phase_colors, LABELS):
    ax.plot(hkl, np.full_like(hkl, y_row), "|", ms=tms, mew=1.1, color=col,
            ls="none", label=lab, zorder=2)

# ---------------- axis limits: automatic unless given ----------------
if X_AXIS == "Q":
    # ticks on clean fractions of Q; in auto mode the outermost ticks sit
    # at equal distances from the two spines
    lo, hi = (min(X_LIMITS), max(X_LIMITS)) if X_LIMITS else (x.min(), x.max())
    max_ticks = max(3, round(2 * ax_w))          # tick density follows the axes width
    for q_step in (0.05, 0.1, 0.2, 0.25, 0.5, 1.0, 2.0, 2.5, 5.0):
        if (hi - lo) / q_step <= max_ticks:
            break
    if X_LIMITS is None:
        first = np.ceil(lo / q_step) * q_step
        last  = np.floor(hi / q_step) * q_step
        margin = max(first - lo, hi - last) + 0.02 * (hi - lo)
        ax.set_xlim(first - margin, last + margin)
    else:
        ax.set_xlim(*X_LIMITS)
    ax.xaxis.set_major_locator(MultipleLocator(q_step))
elif X_LIMITS is None:
    pad = 0.01 * (x.max() - x.min())
    ax.set_xlim(x.min() - pad, x.max() + pad)
else:
    ax.set_xlim(*X_LIMITS)

ax.set_ylim(*(Y_LIMITS if Y_LIMITS else (y_bot, y_top)))

# ---------------- axes cosmetics ----------------
ax.xaxis.set_minor_locator(AutoMinorLocator(4 if X_AXIS == "Q" else None))
ax.tick_params(axis="y", which="both", left=False, right=False,
               labelleft=False, labelright=False)      # bare y axis: spine only
ax.set_xlabel(xlabel)

leg = ax.legend(loc="upper right", fontsize=15, frameon=True,
                facecolor="white", edgecolor="black", framealpha=1,
                handlelength=1.1, handletextpad=0.4,
                borderpad=0.3, labelspacing=0.25, borderaxespad=0.25)
leg.get_frame().set_linewidth(1.0)
leg.set_zorder(10)

out = os.path.join(os.path.dirname(os.path.abspath(DATA_FILE)),
                   f"{STEM}_{'2Th' if X_AXIS == '2theta' else 'Q'}")
fig.savefig(out + ".png", dpi=600, bbox_inches="tight", pad_inches=0.05)
fig.savefig(out + ".pdf", bbox_inches="tight", pad_inches=0.05)
print("Saved:", out + ".png / .pdf")
