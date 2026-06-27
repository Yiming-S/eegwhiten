## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* local macOS, R 4.4.2
* GitHub Actions: ubuntu-latest (R-devel, R-release, R-oldrel-1),
  macOS-latest (R-release), windows-latest (R-release)

## Notes

* `robustbase` is a Suggests dependency used only by the optional
  `cov_estimator = "mcd"` code path; it is guarded with `requireNamespace()`.
