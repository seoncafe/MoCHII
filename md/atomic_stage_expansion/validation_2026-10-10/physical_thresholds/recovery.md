# Threshold Validation Recovery

The first controller stopped before promotion because the diagnostic format omitted the exponent marker for fully stripped carbon (`huge(real64)`). All 32 absorption checks passed before the Python diagnostic reader failed. The format now explicitly reserves three exponent digits. Existing staged products and the completed debug build are reused; the failed controller and its logs are preserved. The corrected controller validates the compiled entries and rejection, then promotes the 11 current-stage files and checks their consumers.
