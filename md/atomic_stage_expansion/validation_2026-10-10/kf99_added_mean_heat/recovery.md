# Added Mean-heat Validation Recovery

The first controller passed 28 Python tests and completed all staged element files. Its Fortran test used the reactant stage index for an output array that the actual API defines by transition index; the first charge-conservation assertion stopped. The test now reads cx_out_stage(it,ie), matching the inspected implementation. Production code was correct. Failed outputs are preserved. A new controller reuses complete staged files through a link, writes fresh compiled outputs, checks all thermal and rate samples, and produces the comparison PDF.
