# Original upper-endpoint diagnostic sequence

This controller is terminal failed. Its direct three-model profiles passed,
but cooling registration failed when a strict numerical T guard rejected
exp(log(100000)) rounding. Sources and products are preserved. The original
analysis tool is saved as `profiles/frozen_analysis_tool.py`; its hash matches
the startup source record. The [corrected continuation](../upper_temperature_models_checked/README.md)
verifies complete profiles and reruns only the unfinished registration checks.
Do not relabel this failed controller as completed or mistake candidate
verification for physical acceptance or production implementation.
