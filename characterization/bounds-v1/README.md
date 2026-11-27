# Isolated legacy checked-access characterization

Task `matrix-characterize-bounds-v1`, role `shared-matrix-daily-implementation`,
run `matrix-20261016-daily-01`, October 16, 2026 UTC and America/Halifax.
Matrix source HEAD: `1500576c382036085bcd5b267c815e9b9117d1ba`.

Run from the canonical Matrix checkout:

```bash
python3 characterization/bounds-v1/run.py --output work/bounds-validation
```

The runner verifies the existing multiplication-v1 source manifest and builds one
legacy header per executable. Pins remain Neuro master
`87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c` and Bitmap development
`4df8a9cecdce77c9e0eb98b0b0897a0c473c0710`. Neither consumer is modified.

Each execution has 24 checks for each of float, double and int (72 total).
Six isolated executables cover Neuro C++17 and Bitmap C++20 in Debug,
Release with NDEBUG, and ASan+UBSan. Checks stay active in Release.
Valid asymmetric values and return-type checks distinguish references from copies.
Endpoint, converted-negative, SIZE_MAX and empty-row cases verify row.at errors.
Neuro additionally covers mutable/const checked matrix and row subscripts,
including empty matrices; Bitmap covers mutable/const matrix.at row/column
endpoints, converted-negative rows and empty matrices. Wrong exception types
and absent exceptions fail the checks. Invalid unchecked Bitmap subscripts are
never evaluated. Converted -1 and SIZE_MAX intentionally represent the same
size_t value, documenting both the caller conversion and maximal-index case.

| Operation | Neuro snapshot | Bitmap snapshot | Proposed adapter requirement |
| --- | --- | --- | --- |
| Valid mutable / const subscript | T& / const T& | T& / const T& | Preserve reference and constness contracts |
| Mutable row.at | T copy | T& alias | Preserve mutation behavior explicitly |
| Const row.at | T copy | T copy | Preserve copy independence |
| Matrix.at | absent | checked mutable/const overloads | Qualify Bitmap adapter API |
| Invalid subscript | throws out_of_range even with NDEBUG | unchecked source access | Never infer portable behavior from invalid unchecked access |
| Invalid row.at | out_of_range | out_of_range | Preserve checked errors |
| Negative index conversion | size_t maximum, checked access rejects | checked at rejects | Do not confuse constructor clamping with indexing |

The plan's signed-negative clamping applies to Neuro's int dimension
constructor. Subscript and at indices use size_t and do not clamp to zero.
This slice does not exercise negative dimensions or attempt huge allocations.
The proposed adapter requirements need independent review; they do not establish
a unified package contract. Exception message strings are not qualified here.

Evidence lives in `work/matrix-20261016-daily-01/validation/results.json`, with
exact commands, counts, exits and logs. Baseline and candidate manifests record
byte preservation. GCC 15.2.1 Linux/libstdc++ is the local qualification only.
Strict Neuro C++20 syntax, Windows/MSVC, remote CI, Pixel and failure paths,
serialization and remaining numeric/API qualification stay open. No extraction,
package version acceptance or downstream migration is established. Existing
consumer owners and Neuro drafts are preserved. Scoped branch publication and a draft PR were authorized on November 1, 2026;
package and consumer acceptance remain open.
