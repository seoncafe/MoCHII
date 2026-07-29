#!/usr/bin/env python3
"""Report the rendered text size of every paper figure as a fraction of the
body text.

A figure's labels are typeset at the size set in its generating script, but
\\includegraphics then scales the whole figure down, so what the reader sees is
(script font size) x (include width / natural figure width).  This script
measures the natural width of each figure PDF, applies the include width taken
from mochii.tex, and prints the resulting size next to the 10 pt body text.
"""
import os
import re
import subprocess

HERE = os.path.dirname(os.path.abspath(__file__))
PAPER = os.path.join(HERE, "..")
TEX = os.path.join(PAPER, "mochii.tex")
FIGDIR = os.path.join(PAPER, "figures")

BODY_PT = 10.0            # \f@size in the two-column body
COLW_TEX = 242.26653      # \columnwidth  [TeX pt]
TEXTW_TEX = 513.11743     # \textwidth    [TeX pt]
TEXPT_PER_IN = 72.27
BP_PER_IN = 72.0

# Font size used by each generating script: (base font.size, scale _S).
SCRIPT_FONT = {
    "g0_gamma_check.pdf":     (10, 1.78),
    "g1_stromgren_check.pdf": (10, 1.58),
    "g2a_thermal_check.pdf":  (10, 1.58),
    "g4_refine_check.pdf":    (10, 1.78),
    "g4_tng_maps.pdf":        (10, 2.02),
    "dust_stromgren.pdf":     (11, 1.63),
    "pahlive_maps.pdf":       (10, 1.35),
    "dustemis_maps.pdf":      (10, 1.35),
    "pdr_profiles.pdf":       (11, 1.58),
    "wood_hii40.pdf":         (10, 1.31),
    "wood_hii20.pdf":         (10, 1.31),
}


def include_widths():
    """Include width [TeX pt] of every figure, read from the paper source."""
    out = {}
    pat = re.compile(r"\\includegraphics\[width=([^\]]+)\]\{figures/([^}]+)\}")
    for m in pat.finditer(open(TEX).read()):
        spec, fig = m.group(1), m.group(2)
        if "columnwidth" in spec:
            frac = spec.replace("\\columnwidth", "").strip()
            out[fig] = (float(frac) if frac else 1.0) * COLW_TEX
        elif "textwidth" in spec:
            frac = spec.replace("\\textwidth", "").strip()
            out[fig] = (float(frac) if frac else 1.0) * TEXTW_TEX
    return out


def natural_width_in(fig):
    out = subprocess.run(["pdfinfo", os.path.join(FIGDIR, fig)],
                         capture_output=True, text=True).stdout
    m = re.search(r"Page size:\s+([\d.]+) x ([\d.]+)", out)
    return float(m.group(1)) / BP_PER_IN


def main():
    inc = include_widths()
    print(f"body text = {BODY_PT:.0f} pt; target for figure text = "
          f"{0.8*BODY_PT:.1f}-{0.9*BODY_PT:.1f} pt\n")
    print(f"{'figure':24s} {'scale':>6s} {'label':>7s} {'title':>7s} "
          f"{'% of body':>10s}  {'':4s}")
    worst = []
    for fig, w_tex in sorted(inc.items()):
        if fig not in SCRIPT_FONT:
            continue
        base, S = SCRIPT_FONT[fig]
        scale = (w_tex / TEXPT_PER_IN) / natural_width_in(fig)
        lab = base * S * scale
        tit = base * 1.2 * S * scale
        frac = 100.0 * lab / BODY_PT
        flag = "ok" if 80.0 <= frac <= 90.0 else "OUT"
        worst.append((abs(frac - 85.0), fig, frac))
        print(f"{fig:24s} {scale:6.3f} {lab:7.2f} {tit:7.2f} {frac:9.1f}%  {flag}")
    worst.sort(reverse=True)
    print(f"\nlargest deviation from 85%: {worst[0][1]} at {worst[0][2]:.1f}%")


if __name__ == "__main__":
    main()
