# Plan: Milne-Consistent Diffuse LyC Sampling

**Status (2026-07-29): implemented for the three ground continua and gated.**
`src/milne_recomb_spectrum_mod.f90` + `tests/milne`; selected by
`par%diffuse_energy_model = 'milne'`.  The default remains `exponential`
pending the HII20 half of the acceptance suite.  The converged HII40 result
and what it settles are in
[HII40_CLOUDY_RESIDUAL_INVESTIGATION.md](HII40_CLOUDY_RESIDUAL_INVESTIGATION.md):
the exact spectrum moves the He I front by +0.00063 pc, 0.35% of the Cloudy
residual, so it is a correctness fix and not the explanation of that residual.

**One deviation from the physical definition below.**  Section "Physical
definition" specifies sampling from the photon-number PDF `j_E/E`.  MoCHII
emits equal-energy packets whose count comes from the channel's ENERGY
luminosity, and with that pairing the packets must follow `E n(E)` instead:
the represented photon rate is then `R <E>_n (1/<E>_n) = R` exactly, whereas
drawing from `n(E)` inflates it by `<E>_n <1/E>_n > 1`.  That inflation was
measured at 1.002 (8000 K) and 1.010 (2x10^4 K) by reinstating the
number-PDF choice, which `tests/milne` rejects.  The mean energy used for the
luminosity is still the number-weighted `<E>_n` given below, so that part of
the specification is unchanged.

## Objective

Replace MoCHII's current ground-recombination diffuse-LyC approximation,

\[
E = E_{\rm th} + kT_e[-\ln(u)],
\]

with continuous sampling from a temperature-dependent, photon-number Milne
PDF for H II → H I(1s), He II → He I(1¹S), and He III → He II(1s). The
sampled energy must remain the energy used in all packet opacities, heating,
and photoionization scores.

The `threshold` and `exponential` models are retained as diagnostic modes;
the new model becomes the production default only after validation.

## Physical definition

For each recombining ion `i`, tabulate the photon-number PDF above its edge:

\[
p_i(E\mid T_e) = \frac{j^{\rm fb}_{E,i}(T_e)/E}
 {\int_{E_{\rm th,i}}^{E_{\max}} j^{\rm fb}_{E',i}(T_e)/E'\,dE'}.
\]

Use the Milne free-bound emissivity with the same ground-state
photoionization cross section used by transport. Include the free-bound Gaunt
factor where available. The normalization determines the shape only;
emission rates remain controlled by the adopted ground-level recombination
coefficient `alpha_1 = alpha_A - alpha_B` unless rate data are deliberately
upgraded in a separate change.

The packet-energy luminosity normalization must use the PDF-consistent mean:

\[
\langle E\rangle_i(T_e) = \int E\,p_i(E\mid T_e)\,dE.
\]

This keeps the number of emitted ground-recombination photons and their
energy luminosity mutually consistent.

## Implementation sequence

1. Add a small `milne_recomb_spectrum_mod` module.
   - Define a common energy coordinate, preferably dimensionless
     `x = (E - E_th)/(kT)` with a safe high-energy tail cutoff.
   - Build separate H I, He I, and He II ground-continuum tables.
   - Tabulate the photon PDF, CDF, mean energy, and normalization on a
     temperature grid covering the code's permitted HII-region range.
   - Make all arrays immutable after initialization, so MPI packet transport
     remains thread-safe and deterministic.

2. Implement the Milne kernel and numerical checks before connecting it to
   packet launch.
   - Evaluate cross sections through the existing H I, He I, and He II
     cross-section routines; do not duplicate threshold constants or fits.
   - Add the free-bound Gaunt factor consistently, or explicitly document a
     first validated `g_fb = 1` implementation if no compatible routine is
     available.
   - Integrate with stable quadrature and require normalized CDF endpoints,
     monotonicity, positive PDF values, and finite mean energy.
   - Compare tabulated `alpha_1` implied by the spectral integral against the
     coefficient used for packet-number emission. Do not silently rescale one
     to hide a disagreement.

3. Add a `par%diffuse_energy_model = 'milne'` option.
   - Preserve `exponential` as the current baseline and `threshold` as the
     sensitivity limit.
   - Validate the input in `setup.f90`; record the selected mode in the
     start-up summary and HDF5 metadata.
   - Initialise the immutable Milne tables once after atomic-data setup, not
     during each iteration.

4. Connect the sampler to `diffuse_mod`.
   - In `diffuse_build`, use the table-interpolated mean energy for each
     ground channel's luminosity.
   - In both `gen_diffuse_photon` and `gen_diffuse_photon_qmc`, invert the
     CDF at the local `gas_Te(il)` to obtain a continuous energy.
   - Interpolate the PDF/CDF in temperature first, then invert once. Do not
     interpolate two independently sampled frequencies, and do not choose a
     nearest temperature row.
   - Keep channel 4 (He I excited cascade) separate in this change; it needs
     level-resolved cascade data and a non-flat two-photon spectrum in a
     follow-up task.

5. Preserve QMC semantics.
   - Continue using diffuse Sobol dimension 8 as the ground-continuum energy
     variate. The dimensionality and stream assignment must not change.
   - Use dimension 9 only for the existing He I two-photon branch.
   - Add deterministic tests that pass fixed uniforms through both random and
     QMC sampler entry points and require identical energy results.

6. Add unit and numerical tests.
   - Threshold behavior: CDF begins at each physical edge and never emits
     below it.
   - Moments: numerical samples reproduce table normalization and mean energy
     for several temperatures.
   - Shape: H I and He I PDFs match independent Milne numerical integrals;
     compare their quantiles, not only their mean energies.
   - Continuity: sampled quantiles and mean energy vary smoothly across every
     temperature-table boundary.
   - Regression: `exponential` and `threshold` modes reproduce their current
     fixed-seed behavior bitwise where possible.

7. Run staged integration tests.
   - Start with a one-cell/controlled H II and He II recombination emitter.
   - Compare packet-energy histograms, mean energies, and emitted channel
     luminosities against the independent integrals.
   - Run HII20 and HII40 with `exponential`, `threshold`, and `milne` using
     the same RQMC seed ensemble and 68 ranks.
   - Apply the existing diagnostic-bin-independence test to the Milne mode.
   - Check MPI decomposition invariance and energy closure.

8. Reassess the Cloudy residual only after convergence.
   - Compare `r(He I=0.5)`, `r(C2+=0.5)`, `r(H I=0.1)`, radial C2+ tails,
     temperature, and line diagnostics with the existing HII40 table.
   - Report the Milne change separately from spatial-resolution, recombination
     rate, and He I excited-cascade effects.
   - Update Figure 14/15 and paper tables only if the production configuration
     is replaced and the full benchmark suite passes.

## Acceptance criteria

- All three ground-continuum CDFs are normalized, monotonic, and use the
  exact transport cross section at the sampled energy.
- Mean energies used in `diffuse_build` agree with the sampler's PDF moments.
- Pseudorandom and QMC paths implement the same distribution.
- `exponential` and `threshold` retain their documented behavior.
- HII20/HII40 outputs converge under the existing physical tolerances, pass
  energy closure and MPI tests, and have recorded reproducible seeds.
- The comparison report distinguishes a physically improved sampler from a
  demonstrated improvement in agreement with Cloudy.

## Deliberately deferred work

- Replacing recombination-rate coefficients or changing `alpha_1` data.
- A full level-resolved He I cascade calculation and a physical He I
  two-photon spectral CDF.
- Non-Maxwellian electron distributions, density-dependent continuum lowering,
  and line-transfer effects.

Those are separate physical changes and should not be mixed with validation
of the ground-continuum Milne sampler.
