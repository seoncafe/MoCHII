"""Shared axis styling for the paper's image maps.

The image panels are square in physical units, so the two axes should carry
the same tick spacing and the same number of decimals.  Matplotlib chooses
the locators independently for x and y, which gives mismatched ticks (for
example 0.1 steps on one axis and 0.05 on the other, or 0.1 against 0.10).
`square_ticks` picks one step for both axes and formats both with the same
number of decimals.
"""

import numpy as np
from matplotlib.ticker import FormatStrFormatter, MultipleLocator

# candidate tick steps, in units of the axis span
_STEPS = np.array([0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.25,
                   0.5, 1.0, 2.0, 2.5, 5.0, 10.0, 20.0, 25.0, 50.0, 100.0,
                   200.0, 500.0, 1000.0])


def square_ticks(ax, target=5):
    """Give the x and y axes one common tick step and decimal count.

    The step is chosen from a round-number ladder so that the wider of the
    two axes carries about `target` intervals; both axes then use it, and
    both are formatted with the decimals that step requires.
    """
    x0, x1 = ax.get_xlim()
    y0, y1 = ax.get_ylim()
    span = max(abs(x1 - x0), abs(y1 - y0))
    if not np.isfinite(span) or span <= 0.0:
        return
    step = _STEPS[np.argmin(np.abs(_STEPS - span / float(target)))]
    # decimals needed to write the step exactly (capped at 3)
    dec = 0
    while dec < 3 and abs(step * 10 ** dec - round(step * 10 ** dec)) > 1e-9:
        dec += 1
    fmt = FormatStrFormatter("%%.%df" % dec)
    for axis in (ax.xaxis, ax.yaxis):
        axis.set_major_locator(MultipleLocator(step))
        axis.set_major_formatter(fmt)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
