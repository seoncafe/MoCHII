#!/usr/bin/env python3
"""Operation-count benchmark for threshold-sorted metal cross sections."""

from pathlib import Path

import numpy as np


REPO = Path(__file__).resolve().parents[2]
ATOMIC = REPO / "data" / "atomic"
KB_EV = 8.617333262145e-5
EMIN = 13.598
EMAX = 100.0
HII_ELEMENTS = ("c", "n", "o", "ne", "s")
ALL_ELEMENTS = ("c", "n", "o", "ne", "s", "ar", "mg", "fe", "si", "cl", "ca")


def thresholds(element: str) -> list[float]:
    values = []
    path = ATOMIC / f"element_{element}.txt"
    for line in path.read_text(encoding="ascii").splitlines():
        fields = line.split()
        if fields and fields[0] == "ETH":
            values.append(float(fields[1]))
    return values


def expected_calls(elements: tuple[str, ...], temperature: float) -> tuple[int, float]:
    edges = sorted(value for element in elements for value in thresholds(element))
    energy = np.linspace(EMIN, EMAX, 400_001)
    density = energy**3 / np.expm1(energy / (KB_EV * temperature))
    total = np.trapz(density, energy)
    calls = 0.0
    for threshold in edges:
        if threshold <= EMIN:
            calls += 1.0
        elif threshold < EMAX:
            mask = energy >= threshold
            calls += np.trapz(density[mask], energy[mask]) / total
    return len(edges), calls


def main() -> int:
    print("Metal packet-cache cross-section operation benchmark")
    print("Spectrum: band-limited Planck energy PDF, 13.598--100 eV")
    for temperature in (20_000.0, 40_000.0):
        for label, elements in (
            ("HII C/N/O/Ne/S", HII_ELEMENTS),
            ("all registry elements", ALL_ELEMENTS),
        ):
            full, sorted_calls = expected_calls(elements, temperature)
            reduction = 100.0 * (1.0 - sorted_calls / full)
            print(
                f"T={temperature:7.0f} K  {label:22s}: "
                f"all={full:2d}, sorted={sorted_calls:6.3f}, "
                f"cross-section calls reduced by {reduction:5.1f}%"
            )
    print(
        "Decision: use the exact threshold-sorted compact loop; "
        "do not introduce stochastic transition subsampling."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
