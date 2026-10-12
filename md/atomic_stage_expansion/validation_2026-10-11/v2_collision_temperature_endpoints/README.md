# MOCASSIN v2 collision-strength temperature endpoints

## Verified source behavior

The actual `emission_mod::equilibrium` selects three cases. One temperature
uses a constant effective collision strength. Two temperatures use a linear
expression in log10(T) without an explicit endpoint clamp in that branch.
Three or more temperatures call the natural cubic spline and `splint`.
The latter explicitly returns the first strength below the first temperature
and the final strength above the last temperature.

This behavior was checked in both the invoked local v2 tree and the original
`mocassin-mocassin.2.02.73.2.tar.gz`. Exact original source bytes are preserved
here, with hashes in [the source review](source_review.json). These `.f90`
files are external source retained for direct comparison and were not edited.
No original-source endpoint probe or complete v2 calculation was run in this
source review. The archive confirms that the endpoint rule is not introduced
by local interpolation optimizations.

## Consequence for the reference model

The measured supplied-table maxima below 100000 K are C I and S III at
20000 K, Ar VI at 40000 K, and Mg I at 10000 K. The current coherent MoCHII
candidates reject temperatures above their supplied final nodes, except for
unaccepted separately prepared constant-strength extensions. An identical
level count therefore does not establish identical usable temperature scope.
The user accepted cold constant-strength continuation; hot continuation has
not been inferred from that acceptance. A v2-style upper endpoint extension
would need explicit approximation bounds and numerical checks and would not
establish high-temperature atomic accuracy.

Mg I's 1000-100000 K constant-strength candidate is already numerically
checked, but its hot physical choice remains pending. No production source
bounds or default input changed in this review. See
[the direct inventory](../shared_excitation_model_inventory/inventory.json).
