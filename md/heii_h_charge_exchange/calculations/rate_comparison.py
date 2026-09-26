#!/usr/bin/env python3
"""Evaluate the He+ + H formulas inspected on 2026-09-25.

This is a formula comparison, not a call into MoCHII or a validation of
the KF96 spin normalization. The quarter-scaled result is hypothetical.
"""

import csv
import math
import sys


def radiative_HeII_H(temperature):
    """Glover & Jappsen (2007), R26, in cm^3 s^-1."""
    return 1.25e-15 * (temperature / 300.0) ** 0.25


def nonradiative_HeII_H(temperature):
    """KF96 term with the upper temperature clamp used by MoCHII."""
    t4 = min(temperature, 1.0e5) / 1.0e4
    return 7.47e-15 * t4**2.06 * (1.0 + 9.93 * math.exp(-3.89 * t4))


def main():
    writer = csv.writer(sys.stdout, lineterminator="\n")
    writer.writerow(
        [
            "T_K",
            "kf96_usage",
            "radiative_cm3_s",
            "nonradiative_as_implemented_cm3_s",
            "total_as_implemented_cm3_s",
            "total_with_quarter_nonradiative_cm3_s",
            "quarter_total_ratio",
        ]
    )
    for temperature in (100, 1000, 3000, 5000, 6000, 8000, 10000,
                        20000, 30000, 100000, 200000):
        if temperature < 6000:
            usage = "extrapolated_below_fit_range"
        elif temperature > 100000:
            usage = "clamped_to_upper_fit_limit"
        else:
            usage = "within_fit_range"
        radiative = radiative_HeII_H(temperature)
        nonradiative = nonradiative_HeII_H(temperature)
        total = radiative + nonradiative
        quarter_total = radiative + 0.25 * nonradiative
        writer.writerow(
            [temperature, usage]
            + [f"{value:.10e}" for value in
               (radiative, nonradiative, total, quarter_total)]
            + [f"{quarter_total / total:.10f}"]
        )


if __name__ == "__main__":
    main()
