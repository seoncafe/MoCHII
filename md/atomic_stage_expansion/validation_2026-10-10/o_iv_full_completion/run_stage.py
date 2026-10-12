"""Run one gated private O IV full-bound source or cooling product."""

import argparse
import hashlib
import json
import math
from pathlib import Path
import subprocess
import sys

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
BASE = HERE / "candidate"
FITTING = ROOT / "tools/fitting"
sys.path.insert(0, str(FITTING))
import chianti_cooling as C  # noqa: E402
import fit_cooling as F  # noqa: E402
from validate_stage_cooling import coefficients  # noqa: E402
from write_excitation_atom import _source_files, _selection_code_files  # noqa: E402
from validate_source_excitation_atom import source_subset  # noqa: E402


OXYGEN_BUILD = HERE.parent / "oxygen_valence_convergence"
SOURCE_DRIVER = OXYGEN_BUILD / "optimized_build/source_check.x"
REGISTRY_DRIVER = OXYGEN_BUILD / "optimized_build/registered_check.x"
COOLING_BUILD = HERE.parent / "bt_floor_cooling_products"
COOLING_DRIVER = COOLING_BUILD / "build_mpi/check_cooling_fit.x"
COUNTS = (147, 2319, 10731, 160965)
DENSITY = [1e-6, 100.0, 1e6, 1e10]
FIT_TARGET = 0.003


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save_result(path, data):
    if path.exists():
        raise FileExistsError(f"Preserve existing O IV product: {path}")
    path.write_text(json.dumps(data, indent=2) + "\n")


def finite_nonnegative(value, bound):
    return isinstance(value, (int, float)) and math.isfinite(value) and 0 <= value < bound


def checked_source():
    levels, radiation, collisions = C.bound_excitation_atom("o", 4)
    counts = (len(levels), len(radiation), len(collisions),
              sum(len(row["xs"]) for row in collisions))
    if counts != COUNTS or list(levels) != list(range(1, 148)):
        raise ValueError(f"O IV full-bound source changed: {counts}")
    if (len({(row[0], row[1]) for row in radiation}) != COUNTS[1]
        or any(rate <= 0 or levels[upper][1] <= levels[lower][1]
               for lower, upper, _, _, rate in radiation)):
        raise ValueError("O IV source radiation has an invalid physical direction")
    if len({tuple(sorted((row["ll"], row["ul"]))) for row in collisions}) != COUNTS[2]:
        raise ValueError("O IV source collision pair identity changed")
    if any(row.get("bt_tmin_K") != 3200.0 for row in collisions):
        raise ValueError("O IV published BT lower bound is absent from a collision row")
    return levels, radiation, collisions


def exact_collision_nodes(atom_path):
    """Prove every serialized BT abscissa and ordinate round-trips exactly."""
    lines = [line.strip() for line in atom_path.read_text().splitlines()
             if line.strip() and not line.lstrip().startswith("#")]
    if lines[:2] != ["NLEVEL_VERSION 2", "NLEV 147"]:
        raise ValueError("Unexpected O IV native atom header")
    cursor = 2 + COUNTS[0]
    if lines[cursor] != "NRAD 2319":
        raise ValueError("O IV native radiation count changed")
    cursor += 1 + COUNTS[1]
    if lines[cursor] != "NODES 160965" or lines[cursor + 1] != "NUPS 10731":
        raise ValueError("O IV total BT nodes or collision rows changed")
    cursor += 2
    levels, _, collisions = source_subset("o", 4, list(range(1, 148)))
    checked = 0
    for row in collisions:
        header = lines[cursor].split()
        cursor += 1
        if len(header) != 8 or header[2] != "BT_SOURCE":
            raise ValueError("Invalid O IV serialized collision header")
        lower, upper = C.excitation_levels(levels, row)
        count = len(row["xs"])
        if ((int(header[0]), int(header[1]), int(header[3]), int(header[7]))
            != (lower, upper, row["ttype"], count)
            or float(header[4]) != row["cups"]
            or float(header[5]) != row["de"]):
            raise ValueError("O IV BT pair or scaling parameters changed")
        for abscissa, ordinate in zip(row["xs"], row["ys"]):
            fields = lines[cursor].split()
            cursor += 1
            if len(fields) != 2 or float(fields[0]) != abscissa or float(fields[1]) != ordinate:
                raise ValueError("O IV original BT node value changed")
        checked += count
    if checked != COUNTS[3]:
        raise ValueError("O IV original BT node count changed")
    expected_tail = ["DENSITY_RATIO 1"]
    source_floor = collisions[0].get("bt_tmin_K")
    if source_floor is not None:
        expected_tail.insert(0, f"BT_TMIN_K {source_floor:.17g}")
    if lines[cursor:] != expected_tail:
        raise ValueError("O IV atom has incomplete or extra numerical records")
    return checked


def population_profiles(archive_path, output_path):
    """Plot all 147 levels and native-source absolute population differences."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    with np.load(archive_path) as archive:
        temperature = archive["temperature_K"]
        if (temperature.shape != (401,) or not np.all(np.isfinite(temperature))
            or not np.all(temperature > 0)):
            raise ValueError("Incomplete O IV temperature profile")
        maxima = []
        with PdfPages(output_path) as pages:
            for index, density in enumerate(DENSITY):
                native = archive[f"populations_native_{index}"]
                source = archive[f"populations_source_{index}"]
                if (native.shape != (401, 147) or source.shape != (401, 147)
                    or not np.all(np.isfinite(native)) or not np.all(np.isfinite(source))
                    or np.min(native) < -1e-13 or np.min(source) < -1e-13):
                    raise ValueError("Incomplete O IV 147-level population profile")
                absolute = abs(native - source)
                maxima.append(float(np.max(absolute)))
                figure, axes = plt.subplots(2, 1, figsize=(10, 8),
                                             constrained_layout=True)
                extent = (math.log10(temperature[0]), math.log10(temperature[-1]),
                          147.5, 0.5)
                fields = ((np.log10(np.maximum(native.T, 1e-300)),
                           r"$\log_{10}$ native level fraction"),
                          (np.log10(np.maximum(absolute.T, 1e-300)),
                           r"$\log_{10}$ absolute native-source difference"))
                for axis, (field, label) in zip(axes, fields):
                    plotted = axis.imshow(field, aspect="auto", interpolation="nearest",
                                          extent=extent, origin="upper")
                    axis.set_ylabel("Source level ID")
                    axis.set_xlabel(r"$\log_{10}(T/{\rm K})$")
                    figure.colorbar(plotted, ax=axis, label=label)
                figure.suptitle(f"O IV, electron density {density:g} cm^-3")
                pages.savefig(figure)
                plt.close(figure)
    if any(not finite_nonnegative(value, 1e-10) for value in maxima):
        raise ValueError("O IV population profiles exceed the native/source tolerance")
    return maxima


def preflight():
    levels, radiation, collisions = checked_source()
    sources = [*_source_files("o", 4), *_selection_code_files("o", 4)]
    raw_hashes = {str(path.resolve()): digest(path) for path in sources}
    save_result(BASE / "preflight_result.json", dict(
        passed=True, stage="preflight", ion="OIV", counts=COUNTS,
        selected_ids=list(levels), source_sha256=raw_hashes,
        lowest_positive_gap_cm_inverse=min(
            levels[upper][1] - levels[lower][1]
            for lower, upper, *_ in radiation),
        adopted_BT_cold_floor_K=3200.0,
        published_BT_calculation_interval_K=[3200.0, 32000000.0],
        source_accuracy_verified=False, physical_model_accepted=False,
        numerical_interval_K=[1000.0, 100000.0]))


def native_profiles():
    source = json.loads((BASE / "preflight_result.json").read_text())
    if source.get("passed") is not True or tuple(source["counts"]) != COUNTS:
        raise ValueError("O IV source preflight is incomplete")
    checked_source()
    output = BASE / "native"
    if output.exists():
        raise FileExistsError(output)
    command = [sys.executable, str(FITTING / "validate_source_excitation_atom.py"),
               "--element", "o", "--stage", "4", "--levels", "147",
               "--points", "401", "--temperature-min", "1000",
               "--temperature-max", "100000", "--driver", str(SOURCE_DRIVER),
               "--work-dir", str(output)]
    with (BASE / "native.log").open("w") as stream:
        subprocess.run(command, cwd=ROOT, stdout=stream,
                       stderr=subprocess.STDOUT, check=True, timeout=1800)
    result = json.loads((output / "result.json").read_text())
    expected = (result.get("nlev"), result.get("nrad"), result.get("nups"), result.get("nodes"))
    if (result.get("passed") is not True or expected != COUNTS
        or result.get("temperature_points") != 401
        or result.get("densities_cm3") != DENSITY
        or result.get("selected_ids") != list(range(1, 148))
        or result.get("full_source_cutoff_errors") != []
        or not finite_nonnegative(result.get("maximum_relative_strength_error"), 1e-10)
        or any(not finite_nonnegative(value, 1e-10)
               for key in ("maximum_relative_cooling_errors",
                           "maximum_absolute_population_errors")
               for value in result.get(key, []))
        or any(len(result.get(key, [])) != 4 for key in (
            "maximum_relative_cooling_errors", "maximum_absolute_population_errors"))):
        raise ValueError("Full O IV native/source numerical comparison failed")
    products = [output / name for name in (
        "nlevel_o_4.txt", "native_strengths.txt", "native_profiles.txt",
        "source_native_profiles.npz", "cooling_profiles.pdf",
        "native_differences.pdf", "result.json")]
    if any(not path.is_file() or path.stat().st_size == 0 for path in products):
        raise ValueError("O IV source profile products are incomplete")
    if result["product_sha256"] != digest(output / "nlevel_o_4.txt"):
        raise ValueError("O IV generated source atom changed")
    nodes = exact_collision_nodes(output / "nlevel_o_4.txt")
    population_plot = output / "all_147_level_populations.pdf"
    population_errors = population_profiles(
        output / "source_native_profiles.npz", population_plot)
    products.append(population_plot)
    save_result(BASE / "native_result.json", dict(
        passed=True, stage="native", ion="OIV", source_result=result,
        original_BT_nodes_verified=nodes,
        maximum_absolute_population_errors_from_saved_profiles=population_errors,
        products_sha256={str(path): digest(path) for path in products},
        source_accuracy_verified=False, physical_model_accepted=False))


def fit_and_reader():
    native = json.loads((BASE / "native_result.json").read_text())
    if native.get("passed") is not True:
        raise ValueError("O IV full-bound source profile is required")
    atom = BASE / "native/nlevel_o_4.txt"
    if native["source_result"]["product_sha256"] != digest(atom):
        raise ValueError("O IV source atom changed before fitting")
    checked_source()
    output = BASE / "fit"
    if output.exists():
        raise FileExistsError(output)
    output.mkdir()
    original_target = F.ERR_TARGET
    try:
        F.ERR_TARGET = FIT_TARGET
        F.generate({"OIV": ("o", 4)}, str(output), plot_path=None)
    finally:
        F.ERR_TARGET = original_target
    fit_path = output / "cooling_o_4.txt"
    amplitude, excitation = coefficients(fit_path)
    temperature = np.geomspace(1000.0, 100000.0, 801)
    source = C.cooling_low_density_limit("o", 4, temperature)
    fitted = F.evaluate(amplitude, excitation, temperature)
    if (not np.all(np.isfinite(source)) or not np.all(np.isfinite(fitted))
        or np.any(source <= 0) or np.any(fitted <= 0)):
        raise ValueError("O IV source or fit has nonpositive cooling")
    difference = np.expm1(np.log(fitted) - np.log(source))
    error = float(np.max(np.abs(difference)))
    if not np.all(np.isfinite(difference)) or error >= FIT_TARGET:
        raise ValueError(f"O IV serialized cooling fit exceeds 0.3%: {error}")
    profile_path = output / "fit_profiles.npz"
    np.savez_compressed(profile_path, temperature_K=temperature,
                        source_cooling=source, fitted_cooling=fitted,
                        relative_difference=difference)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    figure, axes = plt.subplots(2, 1, sharex=True, figsize=(6.5, 6.5),
                                gridspec_kw={"height_ratios": [3, 1]})
    axes[0].loglog(temperature, source, label="Full 147-state source")
    axes[0].loglog(temperature, fitted, "--", label="Serialized fit")
    axes[0].set_ylabel(r"$\Lambda$ [erg cm$^3$ s$^{-1}$]")
    axes[0].legend()
    axes[1].semilogx(temperature, 100 * difference)
    axes[1].set_xlabel(r"$T$ [K]")
    axes[1].set_ylabel(r"Fit difference [\%]")
    figure.tight_layout()
    plot_path = output / "OIV_fit.pdf"
    figure.savefig(plot_path)
    plt.close(figure)
    build = json.loads((COOLING_BUILD / "build_result.json").read_text())
    if (build.get("passed") is not True
        or build["executables"].get("check_cooling_fit.x") != digest(COOLING_DRIVER)
        or any(build["source_sha256"].get(str((ROOT / "src" / name).resolve()))
               != digest(ROOT / "src" / name)
               for name in ("define.f90", "cooling_mod.f90"))):
        raise ValueError("Actual GNU coefficient-reader build changed")
    native_dir = BASE / "native_reader"
    command = [sys.executable, str(FITTING / "validate_native_cooling_coefficients.py"),
               "--coefficients", str(fit_path), "--driver", str(COOLING_DRIVER),
               "--work-dir", str(native_dir)]
    with (output / "native_reader.log").open("w") as stream:
        subprocess.run(command, cwd=ROOT, stdout=stream,
                       stderr=subprocess.STDOUT, check=True, timeout=180)
    reader = json.loads((native_dir / "result.json").read_text())
    if (reader.get("passed") is not True or reader.get("points") != 801
        or not finite_nonnegative(reader.get("maximum_resolved_relative_error"), 1e-9)):
        raise ValueError("Actual GNU coefficient reader failed")
    save_result(BASE / "fit_result.json", dict(
        passed=True, stage="fit", ion="OIV", source_atom_sha256=digest(atom),
        temperature_points=801, temperature_K=[1000.0, 100000.0],
        maximum_relative_fit_error=error, fit_terms=len(amplitude),
        coefficient_sha256=digest(fit_path), source_fit_profiles_sha256=digest(profile_path),
        fit_plot_sha256=digest(plot_path), native_reader=reader,
        cooling_reader_executable_sha256=digest(COOLING_DRIVER),
        cooling_reader_build_sha256=digest(COOLING_BUILD / "build_result.json"),
        source_accuracy_verified=False, physical_model_accepted=False))


def register():
    native = json.loads((BASE / "native_result.json").read_text())
    ledger = json.loads((BASE / "ledger_result.json").read_text())
    fit = json.loads((BASE / "fit_result.json").read_text())
    if any(item.get("passed") is not True for item in (native, ledger, fit)):
        raise ValueError("O IV source, energy, and fit checks must pass before registration")
    atom = BASE / "native/nlevel_o_4.txt"
    fit_path = BASE / "fit/cooling_o_4.txt"
    if (native["source_result"]["product_sha256"] != digest(atom)
        or fit["coefficient_sha256"] != digest(fit_path)):
        raise ValueError("O IV atom or coefficient file changed")
    manifest_path = OXYGEN_BUILD / "optimized_build_manifest.json"
    build = json.loads(manifest_path.read_text())
    linked = [REGISTRY_DRIVER] + [REGISTRY_DRIVER.parent / name for name in (
        "define.o", "nlevel_mod.o", "nlevel_cooling_mod.o")]
    for path in linked:
        relative = str(path.relative_to(ROOT))
        if build["sha256"].get(relative) != digest(path):
            raise ValueError(f"Optimized Intel registry source differs: {path}")
    output = BASE / "registered"
    command = [sys.executable, str(FITTING / "validate_registered_excitation_atom.py"),
               "--source-dir", str(BASE / "native"), "--element", "o", "--stage", "4",
               "--driver", str(REGISTRY_DRIVER), "--work-dir", str(output),
               "--timeout", "1500"]
    with (BASE / "registered.log").open("w") as stream:
        subprocess.run(command, cwd=ROOT, stdout=stream,
                       stderr=subprocess.STDOUT, check=True, timeout=1560)
    data = json.loads((output / "result.json").read_text())
    errors = [data.get(name) for name in (
        "worst_grid_relative_error", "worst_middle_relative_error",
        "worst_apply_relative_error")]
    if (data.get("passed") is not True or data.get("source_sha256") != digest(atom)
        or data.get("samples") != 40 or data.get("table_cells") != 64561
        or tuple(data.get(key) for key in ("nlev", "nrad", "nups")) != COUNTS[:3]
        or any(not finite_nonnegative(value, 1) for value in errors)
        or errors[0] >= 1e-10 or max(errors[1:]) >= FIT_TARGET):
        raise ValueError("Actual O IV registered density table exceeds 0.3%")
    save_result(BASE / "result.json", dict(
        passed=True, ion="OIV", full_source_counts=COUNTS,
        native_result_sha256=digest(BASE / "native_result.json"),
        ledger_result_sha256=digest(BASE / "ledger_result.json"),
        fit_result_sha256=digest(BASE / "fit_result.json"),
        registered_result=data, registered_result_sha256=digest(output / "result.json"),
        atom_sha256=digest(atom), coefficient_sha256=digest(fit_path),
        registry_build_sha256=digest(manifest_path),
        linked_driver_and_object_sha256={str(path): digest(path) for path in linked},
        source_accuracy_verified=False, physical_model_accepted=False,
        production_installed=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("preflight", "native", "fit", "registry"))
    args = parser.parse_args()
    BASE.mkdir(parents=True, exist_ok=True)
    {"preflight": preflight, "native": native_profiles,
     "fit": fit_and_reader, "registry": register}[args.stage]()


if __name__ == "__main__":
    main()
