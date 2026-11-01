# Isolated legacy multiplication characterization

This bounded slice records multiplication behavior for Neuro master
`87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c` and Bitmap development
`4df8a9cecdce77c9e0eb98b0b0897a0c473c0710`. The downloaded legacy files
are exact pinned inputs, verified against `sources.json` before compilation.
Never modify them or compile/link both definitions into one program.

Run from the canonical Matrix checkout:

```bash
python3 characterization/multiplication-v1/run.py --output work/multiplication-validation
```

The runner builds separate executables with GCC, without Google Test or edits
to either consumer checkout. Each execution contains 16 named test groups:
rectangular products, square products, asymmetric column vectors and incompatible
nonempty dimensions for float/double/int (12); cancellation-sensitive products
for float/double, both column and general kernels (4). Exact values and output
shapes are checked. The rectangular case also checks input preservation. The
cancellation control verifies that reassociation changes the answer, so it cannot
silently become a vacuous arithmetic-order test. Inputs avoid integer overflow.

The fixtures reuse the pinned Neuro rectangular/square/cancellation references
and Bitmap square/int reference. The asymmetric column case extends those
references with values representable in all three types. The imported upstream
Google Test files and Bitmap fuzz harness are read-only reference inputs; this
slice does not claim to run their complete suites or fuzz target.

## October 9 results

Source Matrix revision: `1500576c382036085bcd5b267c815e9b9117d1ba`.
Role: `shared-matrix-daily-implementation`.
Run: `matrix-20261009T221711Z`; UTC and America/Halifax work date: 2026-10-09.
Evidence: `work/matrix-20261009T221711Z/validation/results.json` and its command
logs. Exact downloaded source hashes live in `sources.json`; local compiler
library hashes are recorded in the evidence.

GCC 15.2.1 (SUSE build 20251006), Linux x86_64, libstdc++:

| Implementation | Standard | Configuration | Test groups |
| --- | --- | --- | --- |
| Neuro | C++17 | Debug, Release, ASan+UBSan, OpenMP (2 threads) | 16/16 each |
| Bitmap | C++20 | Debug, Release, ASan+UBSan | 16/16 each |

Total: 112 passed groups across seven executions, zero runtime test failures,
and no reported ASan/UBSan diagnostics. These are 16 cases repeated under
configurations, not 112 unique cases. Builds disable floating contraction with
`-ffp-contract=off`; no fast-math/reassociation flags are used. Benchmarking is
not enabled. No performance result is claimed.

A strict `-std=c++20 -pedantic-errors` native compile fails for Neuro at
`matrix.h:422` and `matrix.h:440`: template-id constructor declarations
`Matrix<T>()` and `Matrix<T>(...)` are invalid in C++20. GCC's default C++20
mode accepts these as extensions with warnings. Bitmap's strict native C++20
compile passes. Neuro remains an unmodified C++17 baseline in runtime tests;
this is not proof of Neuro C++20 qualification. No compatibility rewrite was
performed. The runner expects this precise strict diagnostic and fails if its
compatibility outcome changes.

The initial runner exited 1 because it incorrectly expected GCC's non-strict
C++20 mode to reject those declarations. All seven runtime executions passed.
The runner now requests strict mode; the two native probes were rerun and the
unchanged successful runtime/build results reused. `initial-results.json`
preserves that initial failure, and `results.json` explains the reconciliation.

## Contract decisions and remaining gates

Both pinned implementations agree on the tested product shapes, exact values,
left-to-right cancellation references and invalid_argument for incompatible
nonempty dimensions. Preserve these observed arithmetic contracts during future
extraction. This record supplies characterization evidence for independent
review; it does not authorize a new shared contract or consumer migration.

Still open: Pixel/generic lifetimes, copy/move/self-assignment, resize, empty and
bounds behavior, capacity/API differences, Sigmoid semantics, determinants and
integer-division limits, saved-model round trips, full upstream suites/fuzzing,
Windows/MSVC and supported CI qualification. Empty multiplication semantics are
not covered here. No combined legacy header/ODR test was attempted. Extraction,
installed-package tests and independent exact-version acceptance remain later
gates with current owners. No package version was produced or accepted.

All pre-existing tracked files are preserved byte-for-byte, including CRLF-only
working-tree differences. This run adds only this characterization directory and
local evidence. Scoped branch publication and a draft PR were authorized on November 1, 2026;
this creates no release or downstream dependency change.
