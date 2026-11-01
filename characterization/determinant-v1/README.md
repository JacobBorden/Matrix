# Isolated legacy determinant characterization

Task `matrix-characterize-determinant-v1`, implementation run
`matrix-20261013-daily-01`, October 13, 2026 UTC / America/Halifax.
Matrix source HEAD `1500576c382036085bcd5b267c815e9b9117d1ba`.

Run from the canonical Matrix checkout:

```bash
python3 characterization/determinant-v1/run.py --output work/determinant-validation
```

The runner verifies every immutable source hash from the multiplication-v1
manifest: Neuro master `87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c`
and Bitmap development `4df8a9cecdce77c9e0eb98b0b0897a0c473c0710`.
Each executable contains only one legacy Matrix definition. No consumer checkout,
legacy input or existing characterization slice is edited.

Fourteen checks for each of float/double/int and two extra floating checks give
44 checks per execution. Independent analytic constants cover 0x0 determinant 1,
negative scalar, 2x2 determinant 10, the existing Bitmap 3x3 fixture determinant 1,
4x4 upper-triangular determinant -120 with unchanged input, permutation matrices
with one/two row swaps (including an initially zero pivot), duplicate-row and
all-zero singular matrices, 6x6 diagonal product 720, 7x7 diagonal product 5040
or the documented integer rejection, and wide/tall invalid_argument rejection.
Floating tiny-pivot block fixtures have analytic determinant -6; they exercise
elimination stability without borrowing the inverse's singularity threshold.
Fixture integer arithmetic stays inside signed-int range. Checks remain active
in Release/NDEBUG. Floating comparison rejects nonfinite values; singular cases
require exact zero rather than a tolerance that could accept a tiny nonzero value.

## Observed source contracts and compatibility decisions for review

| Behavior | Neuro pin | Bitmap pin | Qualification requirement |
| --- | --- | --- | --- |
| 0x0 | 1 | 1 | Preserve this common convention |
| 1x1 / 2x2 | scalar / direct formula | same | Preserve typed arithmetic |
| 3x3 | recursive minors | direct cofactor formula | Numeric order remains implementation-specific |
| Floating dimension >=4 | recursive cofactor expansion | elimination with partial pivoting | Preserve Bitmap pivoted behavior; do not silently replace either numeric contract |
| Integer dimension 6 | recursive, supported | recursive, supported | Keep bounded arithmetic fixtures |
| Integer dimension 7 | recursive, supported | invalid_argument | Preserve/review adapter-specific limit |
| Nonsquare | invalid_argument | invalid_argument | Preserve rejection type |
| Input | unchanged | elimination uses a copy | Preserve value independence |

These are bounded snapshot observations and proposed adapter requirements,
awaiting independent fixture and contract review. Mathematically equal analytic
fixtures do not establish bitwise equivalence of the two floating algorithms.
No randomized ill-conditioning, NaN-input, overflow or large-dimension performance
contract is claimed. Inverse/cofactor size caps and suspected integer-zero division
are a separate open slice; this task never executes undefined integer division.

## Validation and acceptance limits

Evidence is in `work/matrix-20261013-daily-01/validation/results.json` with exact
commands, return codes, per-command logs and counts. Six executions are expected:
GCC Linux/libstdc++, Neuro C++17 and Bitmap C++20, each Debug, Release(NDEBUG),
and ASan+UBSan. No strict Neuro C++20 repair, OpenMP qualification, Windows/MSVC,
full upstream/fuzz suite or exact-candidate remote CI is included.

Actual Pixel, failure paths, other API/empty cases, inverse/division and saved-model
serialization remain open. Extraction, installed/versioned package qualification
and consumer migration require independent acceptance through their current owners.
Neuro drafts and Bitmap/Neuro ownership remain unchanged. Scoped branch publication and a draft PR were authorized on November 1, 2026;
package and downstream acceptance remain open.
