# Isolated legacy storage characterization

Task `matrix-characterize-storage-v1`, role `shared-matrix-daily-implementation`,
run `matrix-20261010-daily-01`, October 10, 2026 UTC and America/Halifax.
Matrix source HEAD: `1500576c382036085bcd5b267c815e9b9117d1ba`.

Run from the canonical Matrix checkout:

```bash
python3 characterization/storage-v1/run.py --output work/storage-validation
```

The runner reuses immutable, hash-verified inputs and manifest from
`../multiplication-v1`: Neuro master `87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c`
and Bitmap development `4df8a9cecdce77c9e0eb98b0b0897a0c473c0710`.
It compiles each header in a separate executable. Consumer checkouts and the
existing candidate stay untouched; neither legacy definition is linked to the other.

Each execution checks 13 groups for each of float, double, int and a nonnumeric,
copyable instrumented Value, plus one destructor-balance group (53 total): default
empty shape, shape/values, capacity units, independent deep copy construction and
assignment, self-copy, move construction and source reuse, move assignment,
self-move, rectangular growth/default values, shrink, zero-column resize and
throwing row.at bounds. Move tests check transfer of the original element address;
self-assignment checks its preservation. Lifetime accounting covers all normal
scope exits and replaced destinations. Arithmetic APIs are never instantiated for
Value. Value is a generic lifetime probe, not Bitmap's actual Pixel class.

## Observed contracts and proposed compatibility decisions

| Behavior | Neuro baseline | Bitmap baseline | Extraction requirement |
| --- | --- | --- | --- |
| Matrix capacity | bytes = rows × cols × sizeof(T) | row count | Decide explicit adapter units; do not unify silently |
| Row capacity | bytes = elements × sizeof(T) | element count | Preserve each public adapter's units |
| resize(2,0) | clears to 0×0; capacity 0 | retains 2×0; capacity 2 | Review empty-shape adapter behavior |
| Copy | independent deep values | independent deep values | Preserve ordinary copy independence |
| Move | transfers element address, clears source | same | Preserve usable/destructible moved-from state |
| Self-copy/self-move | values and address retained | same | Preserve observed self-assignment behavior |
| Grow/shrink | common rectangular region retained; new values default | same | Preserve value behavior, no pointer stability promised |
| row.at out of bounds | out_of_range | out_of_range | Preserve checked access |

These decisions are proposals for independent review, not an accepted shared
package contract. Bitmap operator[] is unchecked at both matrix and row levels;
Neuro source contains explicit throws. Invalid Bitmap subscripts are not executed,
since undefined behavior cannot establish a portable compatibility requirement.
No contiguous whole-matrix layout is assumed or changed.

## Validation and limits

Evidence: `work/matrix-20261010-daily-01/validation/results.json`, with exact
commands, exit codes, nonzero group counts and per-command logs. GCC 15.2.1,
Linux x86_64/libstdc++, Neuro C++17 and Bitmap C++20: separate Debug, Release
(NDEBUG), and ASan+UBSan executables. Intentional self-move produces compiler
warnings. Neuro's already characterized strict C++20 constructor-syntax failure
is unchanged; this storage slice does not repair or requalify it.

Actual Pixel, throwing-copy/allocation failure and exception guarantees, negative
constructor/index conversion, other empty shapes/splitting, Sigmoid, determinant,
integer division and saved-model serialization remain open. No full upstream test
suite, fuzz execution, Windows/MSVC, OpenMP storage or remote CI acceptance is
claimed. Extraction, package install/export/version qualification and consumer
migration remain gated by independent exact-version acceptance and current owners.
Scoped branch publication and a draft PR were authorized on November 1, 2026;
this does not establish package or consumer acceptance.
