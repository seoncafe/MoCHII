# Native Build Recovery

The first build failed because unchanged IPO objects (first reported: random_mt.o) remained in a native debug link. All original logs are preserved. The Makefile now records compiler settings as an object prerequisite, so changing build mode recompiles every object. This is an inspected build failure, not an uncertain interrupted calculation. The successor workflow rebuilds and then validates the resulting executable and focused checks.
