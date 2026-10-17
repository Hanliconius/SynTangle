# GENESPACE interface

The supported R package is in [r/SynTangleR](../../r/SynTangleR). Install it with `R CMD INSTALL r/SynTangleR`, then `library(SynTangleR)`.

The convenience loader `flip_assist.R` can also be sourced directly from a complete checkout. Both routes expose `syntangle_flip_assist()`, `syntangle_run()` and `syntangle_export_genespace()`.

See the [user guide](../../docs/user-guide/genespace.md) for installation and examples. `test_adapter.R` tests schema validation and, when `SYNTANGLE_TEST_EXECUTABLE` is set, both installed Python workflows from R.
