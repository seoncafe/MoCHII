# Email draft to B. T. Draine — He I 1.0833 um optical-depth factor of 3

Date: 2026-07-24.  Companion memo (to attach):
`docs/Draine2026_HeI10830_tau.pdf`.  Background: the paper is
Draine (2026), ApJ 999, 3 (arXiv:2601.11710); the manuscript had been
shared with the author of MoCHII before publication.  Full analysis in the
memo; the MoCHII/LaRT gate records are in `docs/HEI_10830_PLAN.md` and
`examples/hei_metastable/lart/`.

---

Subject: A possible factor-of-3 in the He I 1.0833 um optical depths
(Table 1 statistical weights)

Dear Bruce,

Thank you again for sharing your manuscript on resonant scattering of the
He I 1.0833 um triplet — I enjoyed it, and I am sorry I did not catch the
following point when I first read it.

While implementing a metastable 2^3S diagnostic in my Monte Carlo
photoionization code (MoCHII) and cross-checking against the published
version (ApJ 999, 3), I found that the population physics agrees very well
with your Eqs. (1)-(4) (alpha_trip, k_d, and the resulting critical
density and column densities all match my implementation to a few
percent), but the optical-depth normalization appears to be a factor of 3
too large.  Table 1 lists g_u = 3, 9, 15 for the 2^3P^o_J levels, whereas
the fine-structure degeneracies are 2J+1 = 1, 3, 5; correspondingly
Eq. (8) carries (2J+1) where the standard line-center cross section
sigma_0 = (g_u/g_l) A lambda^3 / (8 pi^{3/2} b) gives (2J+1)/3 with
g_l(2^3S_1) = 3.  This propagates to tau_tot = 262 in Eq. (11), where the
standard value is 87.5 per 10^14 cm^-2 at b = 10 km/s, and to 448 -> 149
in Eq. (12).

The check that convinced me is internal to Table 1 itself: the tabulated
A = 1.0216e7 s^-1 and the NIST oscillator strengths
f = 0.0599/0.1797/0.2996 (ratio 1:3:5, sum 0.539) are mutually consistent
only with g_u = 2J+1; with g_u = 3(2J+1) the same A-values would imply
oscillator strengths three times larger.  I also confirmed numerically
with two independent codes (my photoionization code and my resonance-line
transfer code LaRT reproduce each other's 10830 line-center tau on the
same model nebula to 0.5%), and the resulting sigma_0 agrees with the
opacity used in the exoplanet metastable-helium community (Oklopcic &
Hirata 2018; p-winds), which has been validated against transit
observations.

Since the H II regions of interest are optically thick under either
normalization, I believe all of your qualitative conclusions are
unaffected — the tau_tot values quoted for M51 NE-Strip, M17-B, and
NGC 3603 would come down by a factor of 3, as would the tau contours in
the (n_3, Q_48) plane.  Of course, please let me know if I am missing a
convention that reconciles the two.  I have written the cross-check up as
a short memo (attached) and would be happy to send more details.

With best regards,
Kwang-Il
