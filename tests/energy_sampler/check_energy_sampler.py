#!/usr/bin/env python3
"""Statistical and analytic gate for the isolated continuous-energy samplers."""

import math
import pathlib
import subprocess
import sys
import tempfile

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import zeta

ROOT = pathlib.Path(__file__).resolve().parents[2]
EXE = ROOT / "tests/energy_sampler/energy_sampler_dump.x"
MPI_EXE = ROOT / "tests/energy_sampler/energy_sampler_mpi.x"
KB_EV = 8.617333262145e-5
EMIN, EMAX = 13.598, 100.0


def planck_shape(e, temperature):
    x = e / (KB_EV * temperature)
    return e**3 / np.expm1(x)


def planck_stats(temperature):
    norm = quad(lambda e: planck_shape(e, temperature), EMIN, EMAX,
                epsabs=0.0, epsrel=2e-12, points=[24.383, 24.587, 54.418])[0]
    mean = quad(lambda e: e * planck_shape(e, temperature), EMIN, EMAX,
                epsabs=0.0, epsrel=2e-12)[0] / norm
    second = quad(lambda e: e * e * planck_shape(e, temperature), EMIN, EMAX,
                  epsabs=0.0, epsrel=2e-12)[0] / norm

    def cdf(e):
        return quad(lambda x: planck_shape(x, temperature), EMIN, e,
                    epsabs=0.0, epsrel=2e-12)[0] / norm

    quantiles = np.array([brentq(lambda e: cdf(e)-q, EMIN, EMAX)
                          for q in (0.25, 0.5, 0.75)])
    window = cdf(24.587) - cdf(24.383)
    return mean, second, quantiles, window, cdf


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def check_distribution(name, values, expected, qmc, conditional_qmc=False):
    mean, second, quantiles, window, cdf = expected
    got_mean = values.mean()
    got_q = np.quantile(values, [0.25, 0.5, 0.75])
    got_window = np.mean((values >= 24.383) & (values < 24.587))
    variance = second - mean * mean
    if qmc:
        # A component-selected subset is a two-dimensional Sobol projection
        # rather than a complete one-dimensional net; retain a tight but
        # realistic deterministic tolerance for both cases.
        mean_tol = max(3e-5 * mean, 3e-4)
        q_tol = 2.5e-3
        window_tol = (12.0 if conditional_qmc else 1.1) / len(values)
    else:
        mean_tol = 5.0 * math.sqrt(variance / len(values))
        density_scale = max(window, 1.0 / len(values))
        window_tol = 5.0 * math.sqrt(density_scale / len(values))
        q_tol = 0.12
    require(abs(got_mean-mean) < mean_tol,
            f"{name}: mean {got_mean} vs {mean} (tol {mean_tol})")
    require(np.max(np.abs(got_q-quantiles)) < q_tol,
            f"{name}: quantiles {got_q} vs {quantiles}")
    require(abs(got_window-window) < window_tol,
            f"{name}: CII/HeI window {got_window} vs {window}")
    for threshold in (24.383, 24.587, 54.418):
        expected_above = 1.0 - cdf(threshold)
        got_above = np.mean(values >= threshold)
        if qmc:
            threshold_tol = (12.0 if conditional_qmc else 1.1) / len(values)
        else:
            threshold_tol = 5.0 * math.sqrt(
                max(expected_above*(1.0-expected_above), 1.0/len(values))
                / len(values))
        require(abs(got_above-expected_above) < threshold_tol,
                f"{name}: fraction above {threshold} eV failed")
    print(f"PASS {name}: mean={got_mean:.9g}, window={got_window:.9g}")


def main():
    if not EXE.exists():
        raise SystemExit(f"missing {EXE}; run tests/energy_sampler/run_test.sh")
    expected20 = planck_stats(2.0e4)
    expected40 = planck_stats(4.0e4)
    with tempfile.TemporaryDirectory(prefix="mochii-energy-") as tmp:
        subprocess.run([str(EXE), tmp], cwd=ROOT, check=True)
        rnd = np.loadtxt(pathlib.Path(tmp) / "random.dat")
        full_planck = np.loadtxt(pathlib.Path(tmp) / "full_planck.dat")
        qmc = np.loadtxt(pathlib.Path(tmp) / "qmc.dat")
        diffuse = np.loadtxt(pathlib.Path(tmp) / "diffuse.dat")
        tab = np.loadtxt(pathlib.Path(tmp) / "tabulated.dat")
        mpi_arrays = []
        for ranks in (1, 2, 3, 5, 8):
            outfile = pathlib.Path(tmp) / f"mpi_{ranks}.bin"
            subprocess.run(["mpirun", "-np", str(ranks), str(MPI_EXE), str(outfile)],
                           cwd=ROOT, check=True)
            mpi_arrays.append(np.fromfile(outfile, dtype=np.float64))

    check_distribution("random Planck 20 kK", rnd[:, 0], expected20, False)
    check_distribution("random Planck 40 kK", rnd[:, 1], expected40, False)
    check_distribution("Sobol Planck 20 kK", qmc[:, 1], expected20, True)
    check_distribution("Sobol Planck 40 kK", qmc[:, 2], expected40, True)

    # Full Planck Barnett--Canfield: dimensionless analytic moments and
    # independent band-acceptance probabilities.
    x = full_planck / (KB_EV * 2.0e4)
    mean_x = 4.0 * zeta(5) / zeta(4)
    second_x = 20.0 * zeta(6) / zeta(4)
    require(abs(x.mean()-mean_x) < 5.0*math.sqrt((second_x-mean_x**2)/len(x)),
            "Barnett--Canfield mean failed")
    require(abs(np.mean(x*x)-second_x) < 0.16, "Barnett--Canfield second moment failed")
    full_norm = math.pi**4 / 15.0

    def full_cdf(value):
        return quad(lambda y: y**3/np.expm1(y), 0.0, value,
                    epsabs=1e-13, epsrel=2e-12)[0] / full_norm

    full_quantiles = np.array([
        brentq(lambda value: full_cdf(value)-q, 1e-12, 30.0)
        for q in (0.25, 0.5, 0.75)
    ])
    require(np.max(np.abs(np.quantile(x, [0.25, 0.5, 0.75])-full_quantiles)) < 0.025,
            "Barnett--Canfield quantiles failed")
    for temperature in (2.0e4, 4.0e4):
        lo, hi = EMIN/(KB_EV*temperature), EMAX/(KB_EV*temperature)
        expected_accept = quad(lambda y: y**3/np.expm1(y), lo, hi,
                               epsabs=0.0, epsrel=2e-12)[0] / (math.pi**4/15.0)
        got = np.mean((x >= lo) & (x <= hi))
        tol = 5.0*math.sqrt(expected_accept*(1.0-expected_accept)/len(x))
        require(abs(got-expected_accept) < tol,
                f"Barnett--Canfield {temperature:g} K band acceptance failed")
        accepted_energy = x[(x >= lo) & (x <= hi)] * KB_EV * temperature
        expected = expected20 if temperature == 2.0e4 else expected40
        check_distribution(f"rejected full-Planck {temperature/1000:g} kK",
                           accepted_energy, expected, False)
    print("PASS full-Planck Barnett--Canfield moments and band acceptance")

    # The inverse is monotone and is numerically inverse to its stored CDF.
    order = np.argsort(qmc[:, 0])
    require(np.all(np.diff(qmc[order, 1]) > 0.0), "Planck inverse is not monotone")
    require(np.max(np.abs(qmc[:, 0]-qmc[:, 7])) < 2e-14,
            "20 kK CDF round trip failed")
    require(np.max(np.abs(qmc[:, 0]-qmc[:, 8])) < 2e-14,
            "40 kK CDF round trip failed")
    probe = np.linspace(0, len(qmc)-1, 65, dtype=int)
    true20 = np.array([expected20[4](value) for value in qmc[probe, 1]])
    true40 = np.array([expected40[4](value) for value in qmc[probe, 2]])
    require(np.max(np.abs(true20-qmc[probe, 0])) < 5e-8,
            "20 kK adaptive Planck CDF accuracy failed")
    require(np.max(np.abs(true40-qmc[probe, 0])) < 5e-8,
            "40 kK adaptive Planck CDF accuracy failed")

    # Constant, globally linear, and sharply broken/threshold-adjacent tables.
    require(np.max(np.abs(tab[:, 0]-tab[:, 4])) < 2e-14,
            "constant table round trip failed")
    require(np.max(np.abs(tab[:, 0]-tab[:, 5])) < 2e-14,
            "linear table round trip failed")
    require(np.max(np.abs(tab[:, 0]-tab[:, 6])) < 2e-14,
            "broken table round trip failed")
    require(np.max(np.abs(tab[:, 1]-(EMIN+(EMAX-EMIN)*tab[:, 0]))) < 2e-13,
            "constant-spectrum analytic inverse failed")
    require(np.all(np.diff(tab[:, 2:4], axis=0) > 0.0),
            "tabulated inverse is not strictly monotone")
    print("PASS tabulated constant/linear/broken analytic inversion")

    # Existing coordinate contract: u1 component, u2 conditional energy;
    # external energy also u2.  Conditional samples must match their own CDFs.
    selected1 = qmc[:, 6] == 1.0
    require(abs(selected1.mean()-0.3) <= 1.0/len(qmc),
            "multi-source component selection is not Sobol-balanced")
    check_distribution("multi-source component 1", qmc[selected1, 5], expected20,
                       True, conditional_qmc=True)
    check_distribution("multi-source component 2", qmc[~selected1, 5], expected40,
                       True, conditional_qmc=True)
    check_distribution("external u2 Planck 40 kK", qmc[:, 4], expected40, True)
    check_distribution("diffuse-stream u8 Planck 40 kK", diffuse[:, 1], expected40, True)
    require(abs(np.corrcoef(qmc[:, 0], diffuse[:, 0])[0, 1]) < 0.03,
            "stellar and diffuse scrambled energy streams are correlated")
    print("PASS QMC component/external/diffuse coordinate mappings")
    require(all(np.array_equal(mpi_arrays[0], values) for values in mpi_arrays[1:]),
            "QMC energies depend on MPI rank count")
    print("PASS QMC energy bitwise identity at 1/2/3/5/8 MPI ranks")

    print("ENERGY SAMPLER GATE: PASS")


if __name__ == "__main__":
    try:
        main()
    except (AssertionError, subprocess.CalledProcessError) as exc:
        print(f"ENERGY SAMPLER GATE: FAIL: {exc}", file=sys.stderr)
        raise SystemExit(1)
