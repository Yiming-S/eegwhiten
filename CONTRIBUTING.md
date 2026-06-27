# Contributing to eegwhiten

Thanks for your interest in improving **eegwhiten**.

## Reporting issues

Please file bugs and feature requests at
<https://github.com/Yiming-S/eegwhiten/issues>. A minimal reproducible example
(ideally with `reprex::reprex()`) makes problems much faster to fix.

## Pull requests

1. Fork the repository and create a feature branch.
2. Keep code base-R only (the package imports only `stats`, `graphics`, and
   `grDevices`); optional functionality goes behind `requireNamespace()` guards.
3. Add or update `testthat` tests under `tests/testthat/`. New numerical code
   should pin values or invariants, not just "runs without error".
4. Document exported functions with roxygen2 and run
   `roxygen2::roxygenise()` so `man/` and `NAMESPACE` stay in sync.
5. Run `devtools::check()` (or `R CMD check --as-cran`) and make sure it is
   clean (0 errors / 0 warnings) before opening the PR.
6. Update `NEWS.md` with a user-facing summary of the change.

## Code style

Follow the existing style: 2-space indentation, `snake_case` for objects,
explicit argument names, and informative error messages. English in all code
and comments.

By contributing you agree that your contributions are licensed under the
package's MIT license.
