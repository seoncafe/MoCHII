"""The MoCHII icon: a soft mochi that happens to be an H II region.

Regenerates docs/mochii_icon.png.  Run from the repository root:
    python3 tools/python/make_mochii_icon.py
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams['text.usetex'] = False
from matplotlib.patches import Ellipse, Circle, FancyArrowPatch


def mochi(fname, face=True, with_text=True):
    fig, ax = plt.subplots(figsize=(4.2, 4.2), dpi=220)
    bg = "white"                         # outside the front
    fig.patch.set_facecolor(bg)
    ax.set_facecolor(bg)
    ax.set_xlim(-1.3, 1.3)
    ax.set_ylim(-1.35, 1.25)
    ax.set_aspect("equal")
    ax.axis("off")

    # --- inside the ionization front: warm ionized zone (cream)
    ax.add_patch(Circle((0, -0.05), 1.08, color="#fdf3e3", zorder=0))

    # --- soft shadow under the mochi
    ax.add_patch(Ellipse((0, -0.72), 1.5, 0.22, color="#e8d8c6", alpha=0.9))

    # --- mochi body: squishy flattened bun, pale pink
    body = Ellipse((0, -0.10), 1.72, 1.32, color="#f9d9de", zorder=2)
    ax.add_patch(body)
    # bottom blush shading
    ax.add_patch(Ellipse((0, -0.34), 1.50, 0.78, color="#f3c2cb",
                         alpha=0.8, zorder=3))
    # inner warm glow: the ionizing star inside the mochi
    for i in range(30):
        r = 0.52 * (1.0 - i / 30)**0.8
        ax.add_patch(Ellipse((0, -0.02), 1.7 * r, 1.25 * r,
                             color="#ffb36b", alpha=0.045, zorder=4))
    for i in range(18):
        r = 0.22 * (1.0 - i / 18)**0.8
        ax.add_patch(Ellipse((0, -0.02), 1.7 * r, 1.25 * r,
                             color="#fff1b8", alpha=0.10, zorder=5))
    # specular highlight (soft, glossy)
    ax.add_patch(Ellipse((-0.42, 0.28), 0.42, 0.20, angle=20,
                         color="white", alpha=0.85, zorder=6))

    # --- the star: peeking sparkle at the glow center
    x0, y0 = 0.0, -0.02
    for ang, L in [(0, 0.13), (90, 0.13), (180, 0.13), (270, 0.13),
                   (45, 0.075), (135, 0.075), (225, 0.075), (315, 0.075)]:
        t = np.deg2rad(ang)
        ax.plot([x0, x0 + L * np.cos(t)], [y0, y0 + L * np.sin(t)],
                color="white", lw=2.0, alpha=0.95, zorder=7,
                solid_capstyle="round")
    ax.plot(x0, y0, "o", ms=5.5, color="white", zorder=8)

    # --- Monte Carlo photons: dotted paths escaping the mochi, one kink
    paths = [
        [(0.10, 0.06), (0.55, 0.42), (0.95, 0.60)],
        [(-0.12, 0.02), (-0.75, 0.45), (-1.02, 0.78)],
        [(0.05, -0.12), (0.75, -0.55)],
    ]
    for p in paths:
        p = np.asarray(p)
        ax.plot(p[:, 0], p[:, 1], ls=(0, (1, 2.4)), lw=2.0,
                color="#e8843c", alpha=0.9, zorder=9,
                solid_capstyle="round")
        d = p[-1] - p[-2]
        d = d / np.hypot(*d)
        ax.annotate("", xy=p[-1] + d * 0.05, xytext=p[-1] - d * 0.04,
                    arrowprops=dict(arrowstyle="-|>", color="#e8843c",
                                    lw=1.8, alpha=0.95), zorder=9)
    # scattering grains
    ax.plot(0.55, 0.42, "o", ms=4.5, color="#a9805e", zorder=10)
    ax.plot(-0.75, 0.45, "o", ms=4.5, color="#a9805e", zorder=10)

    # --- cute face
    if face:
        ax.plot([-0.30], [-0.28], "o", ms=6, color="#5a4a42", zorder=11)
        ax.plot([0.30], [-0.28], "o", ms=6, color="#5a4a42", zorder=11)
        th = np.linspace(np.deg2rad(200), np.deg2rad(340), 40)
        ax.plot(0.0 + 0.09 * np.cos(th), -0.38 + 0.06 * np.sin(th),
                color="#5a4a42", lw=2.0, zorder=11,
                solid_capstyle="round")
        # blush
        ax.add_patch(Ellipse((-0.52, -0.36), 0.16, 0.09,
                             color="#f2a0ae", alpha=0.9, zorder=11))
        ax.add_patch(Ellipse((0.52, -0.36), 0.16, 0.09,
                             color="#f2a0ae", alpha=0.9, zorder=11))

    if with_text:
        ax.text(0.0, -1.14, "MoCHII", ha="center", va="center",
                fontsize=24, fontweight="bold", color="#c46a4a",
                family="DejaVu Sans")
    fig.savefig(fname, bbox_inches="tight", pad_inches=0.05,
                facecolor=fig.get_facecolor())
    plt.close(fig)
    print("written:", fname)


mochi("docs/mochii_icon.png", face=True, with_text=False)
