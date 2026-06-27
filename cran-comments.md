## R CMD check results

0 errors | 0 warnings | 2 notes

* This is a new submission.
* The package URL <https://yiming-s.github.io/eegwhiten/> currently returns a
  404; the pkgdown documentation site is deployed from CI on the first push and
  the URL resolves once GitHub Pages is enabled for the repository.

## Test environments

* local macOS, R 4.4.2
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macOS-latest (R-release), windows-latest (R-release)

## Notes

* `robustbase` is a Suggests dependency used only by the optional
  `cov_estimator = "mcd"` code path; it is guarded with `requireNamespace()`.
