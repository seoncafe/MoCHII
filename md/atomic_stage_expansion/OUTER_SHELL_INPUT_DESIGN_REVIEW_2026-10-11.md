> Current status: approved on October 11, 2026 and implemented with scoped
> validation. The proposal below is historical. See the
> [current implementation](validation_2026-10-11/approved_excitation_models/README.md)
> and [comparison](MOCASSIN_V2_MOCHII_BIDIRECTIONAL_COMPARISON.md).

# Consistent shell selection: proposed public input

## Concrete proposal

Add `par%metal_inner_shell_absorption`, default `.false.`. The default selects
the validated outer-shell component for the packet opacity, path ionization
estimator, path heat estimator, and direct field integration. It must not
leave additional-shell absorption in transport while discarding that
absorption's ionization and energy consequences.

Explicit `.true.` selects additional shells only within physically accepted
vacancy, charge-product, and energy domains. Unsupported Auger regimes must
raise a diagnostic instead of silently incrementing the charge by one.
This proposal is not a claim that those missing domains are implemented.
The parameter must not expose the known incorrect adjacent-product expression
as an unrestricted physical option.

No existing public argument list or primary atomic record layout changes are
proposed in this step. The current files identify a unique PH1 outer row by
the explicitly assigned physical ionization threshold. A separate full
vacancy implementation would still require an accepted physical subshell map,
charge probabilities, fluorescence and residual excitation, and electron
energy accounting.

## Ready code and evidence

`species_mod.outer_shell_photoionization` is private and included in the
rebuilt production module. Its below-PHINNER calculation is already invoked
by the existing cross-section routine. Its above-PHINNER result remains a
validated component rather than a selected production default.
The [actual linked test](validation_2026-10-11/outer_shell_photoionization/README.md)
passed 1,764 coefficient samples across 32 existing and 66 staged transitions.
The Python expression comparison has maximum relative error 4.44089e-16.
All sampled below-PHINNER values retain the former arithmetic.

## Scope consequences and remaining acceptance tests

This changes an existing default's absorption scope, as requested for
additional physics. It preserves valid outer absorption tails and corrects
the reference's mismatch between multi-shell opacity and outer-window cell
rates. It is therefore not literal numerical reproduction of every v2
stored frequency window. Current production contains only 43 metal stages;
this choice does not activate the expanded files.

Before accepting the selector, test actual omitted/false/true namelist choices,
packet opacity and scored ionization/heat using the same cross section,
threshold and transported-band limits, direct versus cached consumers,
and explicit rejection of each unsupported vacancy regime. Inspect the
actual code that scores and reduces those arrays, including continuum and
dust-emission packets. No matched spatial calculation is implied by the
present coefficient check.

Changing the public parameter type awaits user confirmation under the global
AGENTS.md instruction to ask before changing public inputs. A design question
has been sent. The parameter is not present in production at this checkpoint.
