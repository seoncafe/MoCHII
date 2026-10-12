# Charge-transfer Validation Recovery

The first controller passed 25 Python tests, then stopped before staging or production changes because the new validation driver imported the published-rate reader from the grid module rather than charge_exchange_expansion. The import is corrected. A new controller preserves the failed logs and continues with preparation, build, compiled rates, promotion, and consumer checks.
