# Mg III pilot path correction

The initial controller stopped before numerical work because the private script pointed to `mg/mg_3.elvlc` instead of the existing `mg/mg_3/mg_3.elvlc`. Its failed state, log, and plan remain in the parent directory. The corrected raw path was read and yielded the expected 167 levels, 3,364 radiative rows, and 13,861 collision rows. This new plan freezes the corrected script and all current inputs; no previous step or output is silently reused.
