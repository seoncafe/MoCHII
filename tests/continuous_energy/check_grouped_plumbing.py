#!/usr/bin/env python3
"""Fast Phase-1 gate for grouped metadata and continuous-mode fail-fast."""

import os
import re
import subprocess
import sys
import tempfile
from pathlib import Path

import h5py
import numpy as np


HERE = Path(__file__).resolve().parent
REPO = HERE.parent.parent


def decoded_attr(group, name):
    value = np.asarray(group.attrs[name]).reshape(-1)[0]
    if isinstance(value, bytes):
        return value.decode().strip()
    return str(value).strip()


def run_case(executable, input_text, input_path, work, run_env, ranks=None):
    input_path.write_text(input_text)
    command = [str(executable), str(input_path)]
    if ranks is not None:
        command = ["mpirun", "-np", str(ranks)] + command
    return subprocess.run(
        command,
        cwd=work,
        text=True,
        errors="replace",
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        timeout=60,
        check=False,
        env=run_env,
    )


def check_shadow(label, result, require_scattered, require_metals=True):
    if result.returncode != 0:
        print(result.stdout)
        print(f"{label}: FAIL (run failed)")
        return False
    marker = "ION: grouped H/He direct-rate shadow: PASS"
    metal_marker = "ION: grouped metal direct-rate shadow: PASS"
    opacity_marker = "ION: grouped dynamic-opacity shadow: PASS"
    match = re.search(
        r"scored path segments direct/scattered =\s*(\d+)\s+(\d+)",
        result.stdout,
    )
    direct = int(match.group(1)) if match else -1
    scattered = int(match.group(2)) if match else -1
    passed = (
        marker in result.stdout
        and opacity_marker in result.stdout
        and direct > 0
    )
    if require_metals:
        passed &= metal_marker in result.stdout
    if require_scattered:
        passed &= scattered > 0
    else:
        passed &= scattered == 0
    errors = next(
        (
            line.split(":", 1)[1].strip()
            for line in result.stdout.splitlines()
            if "H/He shadow relative errors" in line
        ),
        "missing",
    )
    metal_errors = next(
        (
            line.split(":", 1)[1].strip()
            for line in result.stdout.splitlines()
            if "metal shadow relative errors" in line
        ),
        "missing",
    )
    opacity_result = next(
        (
            line.split("ION:", 1)[1].strip()
            for line in result.stdout.splitlines()
            if "dynamic-opacity shadow checks" in line
        ),
        "missing",
    )
    print(
        f"{label}: {'PASS' if passed else 'FAIL'} "
        f"(segments={direct}/{scattered}, H/He errors={errors}, "
        f"metal errors={metal_errors}, opacity={opacity_result})"
    )
    return passed


def main():
    executable = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else REPO / "MoCHII.x"
    if not executable.exists():
        print(f"missing executable: {executable}")
        return 2

    source = (HERE / "grouped_plumbing.in").read_text()
    run_env = os.environ.copy()
    # A one-rank smoke needs no network fabric.  Pinning Intel MPI to shared
    # memory also makes the gate usable in restricted CI/sandbox environments.
    run_env.setdefault("I_MPI_FABRICS", "shm")
    with tempfile.TemporaryDirectory(prefix="mochii-energy-plumbing-") as tmp:
        work = Path(tmp)
        grouped = run_case(
            executable, source, work / "grouped.in", work, run_env
        )
        ok = check_shadow("AMR direct shadow", grouped, require_scattered=False)

        with h5py.File(work / "grouped_plumbing_rates.h5", "r") as handle:
            header = handle["Gamma_HI"]
            checks = {
                "IONEMODE": "grouped",
                "NUBINUSE": "solver_and_diagnostic",
                "ENRGSAMP": "random",
            }
            for key, expected in checks.items():
                actual = decoded_attr(header, key)
                passed = actual == expected
                print(f"{key:8s} = {actual!r}: {'PASS' if passed else 'FAIL'}")
                ok &= passed
            cdf_rtol = float(np.asarray(header.attrs["CDFRTOL"]).reshape(-1)[0])
            passed = cdf_rtol == 1.0e-8
            print(f"CDFRTOL  = {cdf_rtol:.9e}: {'PASS' if passed else 'FAIL'}")
            ok &= passed
            shadow_attr = bool(np.asarray(header.attrs["SHDWRATE"]).reshape(-1)[0])
            print(f"SHDWRATE = {shadow_attr}: {'PASS' if shadow_attr else 'FAIL'}")
            ok &= shadow_attr

        no_shadow_source = source.replace(
            "par%ion_shadow_rates = .true.",
            "par%ion_shadow_rates = .false.",
        ).replace("grouped_plumbing.h5", "grouped_no_shadow.h5")
        no_shadow = run_case(
            executable,
            no_shadow_source,
            work / "no_shadow.in",
            work,
            run_env,
        )
        no_shadow_ok = (
            no_shadow.returncode == 0
            and "direct-rate shadow: PASS" not in no_shadow.stdout
        )
        with h5py.File(work / "grouped_plumbing_rates.h5", "r") as enabled, h5py.File(
            work / "grouped_no_shadow_rates.h5", "r"
        ) as disabled:
            for name in enabled:
                no_shadow_ok &= np.array_equal(
                    enabled[f"{name}/data"][:],
                    disabled[f"{name}/data"][:],
                    equal_nan=True,
                )
            disabled_attr = bool(
                np.asarray(disabled["Gamma_HI"].attrs["SHDWRATE"]).reshape(-1)[0]
            )
            no_shadow_ok &= not disabled_attr
        print(
            "shadow off/on solver arrays bitwise identical:"
            f" {'PASS' if no_shadow_ok else 'FAIL'}"
        )
        ok &= no_shadow_ok

        hhe_source = source.replace(
            "par%use_metals       = .true.",
            "par%use_metals       = .false.",
        ).replace("grouped_plumbing.h5", "grouped_hhe.h5")
        hhe = run_case(executable, hhe_source, work / "hhe.in", work, run_env)
        hhe_ok = check_shadow(
            "H/He-only opacity",
            hhe,
            require_scattered=False,
            require_metals=False,
        )
        hhe_ok &= "grouped metal direct-rate shadow: PASS" not in hhe.stdout
        ok &= hhe_ok

        diffuse_source = (
            source.replace("par%xHI_init         = 1.0", "par%xHI_init         = 0.01")
            .replace("par%xHeI_init        = 1.0", "par%xHeI_init        = 0.10")
            .replace("par%xHeII_init       = 0.0", "par%xHeII_init       = 0.90")
            .replace(
                "par%source_geometry  = 'point'",
                "par%case_ab          = 'A'\n"
                " par%diffuse_field    = .true.\n"
                " par%source_geometry  = 'point'",
            )
            .replace("grouped_plumbing.h5", "grouped_diffuse.h5")
        )
        diffuse = run_case(
            executable, diffuse_source, work / "diffuse.in", work, run_env
        )
        diffuse_ok = check_shadow(
            "grouped diffuse energy", diffuse, require_scattered=False
        )
        match = re.search(
            r"diffuse field: L =\s*[0-9.Ee+-]+\s*erg/s,\s*(\d+)\s*packets",
            diffuse.stdout,
        )
        diffuse_ok &= match is not None and int(match.group(1)) > 0
        print(
            "grouped diffuse packets launched:"
            f" {'PASS' if match is not None and int(match.group(1)) > 0 else 'FAIL'}"
        )
        ok &= diffuse_ok

        dda_source = source.replace(
            "par%car_walk         = 'amr'",
            "par%car_walk         = 'dda'",
        ).replace("grouped_plumbing.h5", "grouped_dda.h5")
        dda = run_case(executable, dda_source, work / "dda.in", work, run_env)
        ok &= check_shadow("DDA direct shadow", dda, require_scattered=False)

        dust_table = REPO / "data/kext_albedo_WD_MW_3.1_60_D03.all_2003"
        dust_source = source.replace(
            "par%dust_model       = 'none'",
            "par%dust_model       = 'global_dgr'",
        ).replace(
            "par%ion_add_dust     = .false.",
            "par%ion_add_dust     = .true.\n"
            " par%ion_dust_scatter = .true.\n"
            f" par%ion_dust_kext    = '{dust_table}'\n"
            " par%cext_dust        = 4.868e-22\n"
            " par%DGR              = 1.0",
        ).replace("grouped_plumbing.h5", "grouped_scattered.h5")
        scattered = run_case(
            executable, dust_source, work / "scattered.in", work, run_env
        )
        ok &= check_shadow(
            "AMR scattered shadow", scattered, require_scattered=True
        )

        refine_source = (REPO / "tests/g4_refine/g4_refine.in").read_text()
        refine_source = (
            refine_source.replace(
                "par%no_photons    = 4.0e7", "par%no_photons    = 2000"
            )
            .replace("par%gas_niter     = 60", "par%gas_niter     = 2")
            .replace("par%refine_iter   = 15", "par%refine_iter   = 1")
            .replace("par%refine_lmax   = 7", "par%refine_lmax   = 6")
            .replace(
                "par%amr_file      = 'amr_L5.fits'",
                f"par%amr_file      = '{REPO / 'tests/g4_refine/amr_L5.fits'}'",
            )
            .replace("par%out_file        = 'g4_refine.h5'", "par%out_file = 'refine.h5'")
            .replace(
                "par%source_geometry = 'point'",
                "par%ion_shadow_rates = .true.\n"
                " par%use_metals = .true.\n"
                " par%source_geometry = 'point'",
            )
        )
        refined = run_case(
            executable, refine_source, work / "refine.in", work, run_env
        )
        refine_ok = (
            "AMR: re-refined on the I-front" in refined.stdout
            and check_shadow("re-refined AMR cache", refined, require_scattered=False)
        )
        ok &= refine_ok

        continuous_source = (
            source.replace(
                "par%ion_energy_mode  = 'grouped'",
                "par%ion_energy_mode  = 'continuous'",
            )
            .replace(
                "par%ion_shadow_rates = .true.",
                "par%ion_shadow_rates = .false.",
            )
            .replace(
                "par%use_metals       = .true.",
                "par%use_metals       = .false.",
            )
            .replace(
                "par%no_photons       = 1000",
                "par%no_photons       = 8192",
            )
            .replace(
                "par%no_print         = 1000",
                "par%no_print         = 8192",
            )
            .replace(
                "par%source_geometry  = 'point'",
                "par%launch_sequence  = 'sobol'\n"
                " par%qmc_seed         = 20260728\n"
                " par%source_geometry  = 'point'",
            )
            .replace("grouped_plumbing.h5", "continuous_8.h5")
        )
        continuous_8 = run_case(
            executable, continuous_source, work / "continuous_8.in", work, run_env
        )
        continuous_32_source = continuous_source.replace(
            "par%nnu_ion          = 8", "par%nnu_ion          = 32"
        ).replace("continuous_8.h5", "continuous_32.h5")
        continuous_32 = run_case(
            executable,
            continuous_32_source,
            work / "continuous_32.in",
            work,
            run_env,
        )
        continuous_ok = (
            continuous_8.returncode == 0
            and continuous_32.returncode == 0
            and "ION: continuous H/He direct rates: ACTIVE" in continuous_8.stdout
            and "ION: continuous H/He direct rates: ACTIVE" in continuous_32.stdout
        )
        solver_fields = (
            "Gamma_HI", "Gamma_HeI", "Gamma_HeII",
            "Heat_HI", "Heat_HeI", "Heat_HeII",
        )
        if continuous_ok:
            with h5py.File(work / "continuous_8_rates.h5", "r") as coarse, h5py.File(
                work / "continuous_32_rates.h5", "r"
            ) as fine:
                for name in solver_fields:
                    continuous_ok &= np.array_equal(
                        coarse[f"{name}/data"][:],
                        fine[f"{name}/data"][:],
                        equal_nan=True,
                    )
                    continuous_ok &= np.all(np.isfinite(coarse[f"{name}/data"][:]))
                header = coarse["Gamma_HI"]
                continuous_ok &= decoded_attr(header, "IONEMODE") == "continuous"
                continuous_ok &= decoded_attr(header, "NUBINUSE") == "diagnostic_only"
                continuous_ok &= decoded_attr(header, "ENRGSAMP") == "sobol"
                emitted = float(np.asarray(header.attrs["L_EMIT"]).reshape(-1)[0])
                absorbed = float(np.asarray(header.attrs["L_ABS"]).reshape(-1)[0])
                escaped = float(np.asarray(header.attrs["L_ESC"]).reshape(-1)[0])
                continuous_ok &= emitted > 0.0
                continuous_ok &= abs(absorbed + escaped - emitted) <= 5.0e-13*emitted
        print(
            "continuous H/He 8/32-bin solver arrays bitwise identical:"
            f" {'PASS' if continuous_ok else 'FAIL'}"
        )
        if not continuous_ok:
            print(continuous_8.stdout)
            print(continuous_32.stdout)
        ok &= continuous_ok

        continuous_metal_source = (
            continuous_source.replace(
                "par%use_metals       = .false.",
                "par%use_metals       = .true.\n"
                " par%ion_metal_abs    = .false.",
            )
            .replace("continuous_8.h5", "continuous_metal_8.h5")
        )
        continuous_metal_8 = run_case(
            executable,
            continuous_metal_source,
            work / "continuous_metal_8.in",
            work,
            run_env,
        )
        continuous_metal_32_source = continuous_metal_source.replace(
            "par%nnu_ion          = 8", "par%nnu_ion          = 32"
        ).replace("continuous_metal_8.h5", "continuous_metal_32.h5")
        continuous_metal_32 = run_case(
            executable,
            continuous_metal_32_source,
            work / "continuous_metal_32.in",
            work,
            run_env,
        )
        metal_ok = (
            continuous_metal_8.returncode == 0
            and continuous_metal_32.returncode == 0
            and "ION: continuous metal direct rates: ACTIVE"
            in continuous_metal_8.stdout
            and "ION: continuous metal direct rates: ACTIVE"
            in continuous_metal_32.stdout
        )
        if metal_ok:
            with h5py.File(
                work / "continuous_metal_8_rates.h5", "r"
            ) as coarse, h5py.File(
                work / "continuous_metal_32_rates.h5", "r"
            ) as fine:
                metal_fields = sorted(
                    name for name in coarse.keys()
                    if name.startswith("x_") and name.endswith("_stages")
                )
                metal_ok &= bool(metal_fields)
                for name in metal_fields:
                    coarse_frac = coarse[f"{name}/data"][:]
                    fine_frac = fine[f"{name}/data"][:]
                    metal_ok &= np.array_equal(
                        coarse_frac, fine_frac, equal_nan=True
                    )
                    metal_ok &= np.all(np.isfinite(coarse_frac))
                    metal_ok &= np.allclose(
                        np.sum(coarse_frac, axis=1),
                        1.0,
                        rtol=5.0e-14,
                        atol=5.0e-14,
                    )
                for prefix in ("Gamma_", "Heat_"):
                    rate_fields = sorted(
                        name for name in coarse.keys()
                        if name.startswith(prefix) and name.endswith("_stages")
                    )
                    metal_ok &= len(rate_fields) == len(metal_fields)
                    for name in rate_fields:
                        coarse_rate = coarse[f"{name}/data"][:]
                        fine_rate = fine[f"{name}/data"][:]
                        metal_ok &= np.array_equal(
                            coarse_rate, fine_rate, equal_nan=True
                        )
                        metal_ok &= np.all(np.isfinite(coarse_rate))
                        metal_ok &= np.all(coarse_rate >= 0.0)
                        metal_ok &= np.any(coarse_rate > 0.0)
        print(
            "continuous metal 8/32-bin rates/heating/stages bitwise identical:"
            f" {'PASS' if metal_ok else 'FAIL'}"
        )
        if not metal_ok:
            print(continuous_metal_8.stdout)
            print(continuous_metal_32.stdout)
        ok &= metal_ok

        continuous_metal_opacity_source = (
            continuous_metal_source.replace(
                "par%ion_metal_abs    = .false.",
                "par%ion_metal_abs    = .true.",
            )
            .replace(
                "par%gas_niter        = 0",
                "par%gas_niter        = 3",
            )
            .replace("continuous_metal_8.h5", "continuous_metal_opacity_8.h5")
        )
        continuous_metal_opacity_8 = run_case(
            executable,
            continuous_metal_opacity_source,
            work / "continuous_metal_opacity_8.in",
            work,
            run_env,
        )
        continuous_metal_opacity_32_source = (
            continuous_metal_opacity_source.replace(
                "par%nnu_ion          = 8", "par%nnu_ion          = 32"
            ).replace(
                "continuous_metal_opacity_8.h5",
                "continuous_metal_opacity_32.h5",
            )
        )
        continuous_metal_opacity_32 = run_case(
            executable,
            continuous_metal_opacity_32_source,
            work / "continuous_metal_opacity_32.in",
            work,
            run_env,
        )
        metal_opacity_ok = (
            continuous_metal_opacity_8.returncode == 0
            and continuous_metal_opacity_32.returncode == 0
        )
        opacity_fields = ()
        if metal_opacity_ok:
            with h5py.File(
                work / "continuous_metal_opacity_8_rates.h5", "r"
            ) as coarse, h5py.File(
                work / "continuous_metal_opacity_32_rates.h5", "r"
            ) as fine:
                opacity_fields = solver_fields + (
                    "x_HI", "x_HeI", "x_HeII", "n_e", "T_e"
                ) + tuple(
                    sorted(
                        name for name in coarse.keys()
                        if (
                            name.endswith("_stages")
                            and name.startswith(("x_", "Gamma_", "Heat_"))
                        )
                    )
                )
                for name in opacity_fields:
                    metal_opacity_ok &= np.array_equal(
                        coarse[f"{name}/data"][:],
                        fine[f"{name}/data"][:],
                        equal_nan=True,
                    )
        print(
            "continuous exact metal opacity 8/32-bin state/rate identity:"
            f" {'PASS' if metal_opacity_ok else 'FAIL'}"
        )
        if not metal_opacity_ok:
            print(continuous_metal_opacity_8.stdout)
            print(continuous_metal_opacity_32.stdout)
        ok &= metal_opacity_ok

        continuous_metal_opacity_np3_source = (
            continuous_metal_opacity_source.replace(
                "continuous_metal_opacity_8.h5",
                "continuous_metal_opacity_np3.h5",
            )
        )
        continuous_metal_opacity_np3 = run_case(
            executable,
            continuous_metal_opacity_np3_source,
            work / "continuous_metal_opacity_np3.in",
            work,
            run_env,
            ranks=3,
        )
        metal_opacity_mpi_ok = (
            metal_opacity_ok and continuous_metal_opacity_np3.returncode == 0
        )
        if metal_opacity_mpi_ok:
            with h5py.File(
                work / "continuous_metal_opacity_8_rates.h5", "r"
            ) as one, h5py.File(
                work / "continuous_metal_opacity_np3_rates.h5", "r"
            ) as three:
                for name in opacity_fields:
                    metal_opacity_mpi_ok &= np.allclose(
                        one[f"{name}/data"][:],
                        three[f"{name}/data"][:],
                        rtol=5.0e-14,
                        atol=0.0,
                        equal_nan=True,
                    )
        print(
            "continuous exact metal opacity 1/3-rank agreement:"
            f" {'PASS' if metal_opacity_mpi_ok else 'FAIL'}"
        )
        if not metal_opacity_mpi_ok:
            print(continuous_metal_opacity_np3.stdout)
        ok &= metal_opacity_mpi_ok

        continuous_iter_source = (
            continuous_source.replace(
                "par%gas_niter        = 0", "par%gas_niter        = 3"
            )
            .replace("continuous_8.h5", "continuous_iter_8.h5")
        )
        continuous_iter_8 = run_case(
            executable,
            continuous_iter_source,
            work / "continuous_iter_8.in",
            work,
            run_env,
        )
        continuous_iter_32_source = continuous_iter_source.replace(
            "par%nnu_ion          = 8", "par%nnu_ion          = 32"
        ).replace("continuous_iter_8.h5", "continuous_iter_32.h5")
        continuous_iter_32 = run_case(
            executable,
            continuous_iter_32_source,
            work / "continuous_iter_32.in",
            work,
            run_env,
        )
        iter_ok = (
            continuous_iter_8.returncode == 0
            and continuous_iter_32.returncode == 0
        )
        if iter_ok:
            with h5py.File(
                work / "continuous_iter_8_rates.h5", "r"
            ) as coarse, h5py.File(
                work / "continuous_iter_32_rates.h5", "r"
            ) as fine:
                for name in solver_fields + (
                    "x_HI", "x_HeI", "x_HeII", "n_e", "T_e"
                ):
                    iter_ok &= np.array_equal(
                        coarse[f"{name}/data"][:],
                        fine[f"{name}/data"][:],
                        equal_nan=True,
                    )
        print(
            "continuous H/He iteration 8/32-bin state identity:"
            f" {'PASS' if iter_ok else 'FAIL'}"
        )
        if not iter_ok:
            print(continuous_iter_8.stdout)
            print(continuous_iter_32.stdout)
        ok &= iter_ok

        table_path = REPO / "tests/continuous_energy/threshold_spectrum.txt"
        continuous_table_source = (
            continuous_source.replace(
                "par%tstar            = 4.0e4",
                f"par%ion_spectrum     = '{table_path}'",
            )
            .replace("continuous_8.h5", "continuous_table_8.h5")
        )
        continuous_table_8 = run_case(
            executable,
            continuous_table_source,
            work / "continuous_table_8.in",
            work,
            run_env,
        )
        continuous_table_32_source = continuous_table_source.replace(
            "par%nnu_ion          = 8", "par%nnu_ion          = 32"
        ).replace("continuous_table_8.h5", "continuous_table_32.h5")
        continuous_table_32 = run_case(
            executable,
            continuous_table_32_source,
            work / "continuous_table_32.in",
            work,
            run_env,
        )
        table_ok = (
            continuous_table_8.returncode == 0
            and continuous_table_32.returncode == 0
        )
        if table_ok:
            with h5py.File(
                work / "continuous_table_8_rates.h5", "r"
            ) as coarse, h5py.File(
                work / "continuous_table_32_rates.h5", "r"
            ) as fine:
                for name in solver_fields:
                    table_ok &= np.array_equal(
                        coarse[f"{name}/data"][:],
                        fine[f"{name}/data"][:],
                        equal_nan=True,
                    )
        print(
            "continuous tabulated spectrum 8/32-bin solver identity:"
            f" {'PASS' if table_ok else 'FAIL'}"
        )
        if not table_ok:
            print(continuous_table_8.stdout)
            print(continuous_table_32.stdout)
        ok &= table_ok

        continuous_np3_source = continuous_source.replace(
            "continuous_8.h5", "continuous_np3.h5"
        )
        continuous_np3 = run_case(
            executable,
            continuous_np3_source,
            work / "continuous_np3.in",
            work,
            run_env,
            ranks=3,
        )
        mpi_ok = continuous_np3.returncode == 0
        if mpi_ok:
            with h5py.File(work / "continuous_8_rates.h5", "r") as one, h5py.File(
                work / "continuous_np3_rates.h5", "r"
            ) as three:
                for name in solver_fields:
                    mpi_ok &= np.allclose(
                        one[f"{name}/data"][:],
                        three[f"{name}/data"][:],
                        rtol=5.0e-14,
                        atol=0.0,
                        equal_nan=True,
                    )
        print(
            "continuous H/He 1/3-rank solver agreement:"
            f" {'PASS' if mpi_ok else 'FAIL'}"
        )
        if not mpi_ok:
            print(continuous_np3.stdout)
        ok &= mpi_ok

        continuous_random_source = (
            continuous_source.replace(
                "par%launch_sequence  = 'sobol'",
                "par%launch_sequence  = 'random'",
            )
            .replace("continuous_8.h5", "continuous_random.h5")
        )
        continuous_random = run_case(
            executable,
            continuous_random_source,
            work / "continuous_random.in",
            work,
            run_env,
        )
        random_ok = (
            continuous_random.returncode == 0
            and "ION: continuous H/He direct rates: ACTIVE"
            in continuous_random.stdout
        )
        if random_ok:
            with h5py.File(work / "continuous_random_rates.h5", "r") as handle:
                for name in solver_fields:
                    random_ok &= np.all(np.isfinite(handle[f"{name}/data"][:]))
                random_ok &= decoded_attr(handle["Gamma_HI"], "ENRGSAMP") == "random"
        print(
            "continuous H/He pseudorandom production path:"
            f" {'PASS' if random_ok else 'FAIL'}"
        )
        if not random_ok:
            print(continuous_random.stdout)
        ok &= random_ok

        continuous_diffuse_source = (
            continuous_source.replace(
                "par%xHI_init         = 1.0", "par%xHI_init         = 0.01"
            )
            .replace(
                "par%xHeI_init        = 1.0", "par%xHeI_init        = 0.10"
            )
            .replace(
                "par%xHeII_init       = 0.0", "par%xHeII_init       = 0.90"
            )
            .replace(
                "par%source_geometry  = 'point'",
                "par%case_ab          = 'A'\n"
                " par%diffuse_field    = .true.\n"
                " par%source_geometry  = 'point'",
            )
            .replace("continuous_8.h5", "continuous_diffuse_8.h5")
        )
        continuous_diffuse_8 = run_case(
            executable,
            continuous_diffuse_source,
            work / "continuous_diffuse_8.in",
            work,
            run_env,
        )
        continuous_diffuse_32_source = continuous_diffuse_source.replace(
            "par%nnu_ion          = 8", "par%nnu_ion          = 32"
        ).replace("continuous_diffuse_8.h5", "continuous_diffuse_32.h5")
        continuous_diffuse_32 = run_case(
            executable,
            continuous_diffuse_32_source,
            work / "continuous_diffuse_32.in",
            work,
            run_env,
        )
        diffuse_cont_ok = (
            continuous_diffuse_8.returncode == 0
            and continuous_diffuse_32.returncode == 0
        )
        diffuse_match = re.search(
            r"diffuse field: L =\s*[0-9.Ee+-]+\s*erg/s,\s*(\d+)\s*packets",
            continuous_diffuse_8.stdout,
        )
        diffuse_cont_ok &= diffuse_match is not None and int(diffuse_match.group(1)) > 0
        if diffuse_cont_ok:
            with h5py.File(
                work / "continuous_diffuse_8_rates.h5", "r"
            ) as coarse, h5py.File(
                work / "continuous_diffuse_32_rates.h5", "r"
            ) as fine:
                for name in solver_fields:
                    diffuse_cont_ok &= np.array_equal(
                        coarse[f"{name}/data"][:],
                        fine[f"{name}/data"][:],
                        equal_nan=True,
                    )
        print(
            "continuous diffuse 8/32-bin solver arrays bitwise identical:"
            f" {'PASS' if diffuse_cont_ok else 'FAIL'}"
        )
        if not diffuse_cont_ok:
            print(continuous_diffuse_8.stdout)
            print(continuous_diffuse_32.stdout)
        ok &= diffuse_cont_ok

        continuous = run_case(
            executable,
            source.replace(
                "par%ion_energy_mode  = 'grouped'",
                "par%ion_energy_mode  = 'continuous'",
            ).replace(
                "par%use_metals       = .true.",
                "par%use_metals       = .true.\n"
                " par%metal_ne         = .true.",
            ),
            work / "continuous.in",
            work,
            run_env,
        )
        guard_text = "continuous H/He vertical slice requires metal_ne=.false."
        passed = continuous.returncode != 0 and guard_text in continuous.stdout
        print(
            "continuous fail-fast:"
            f" {'PASS' if passed else 'FAIL'} (exit={continuous.returncode})"
        )
        ok &= passed

    print("GROUPED ENERGY PLUMBING:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
