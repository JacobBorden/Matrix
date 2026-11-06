# Legacy matrix-wide resize exceptions

Separate executables characterize the unchanged, hash-verified Neuro master
`87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c` and Bitmap development
`4df8a9cecdce77c9e0eb98b0b0897a0c473c0710` inputs recorded in
`../multiplication-v1/sources.json`. Each binary includes one legacy definition.

Run from the Matrix checkout:

```sh
python3 characterization/matrix-resize-failure-v1/run.py --output /absolute/evidence/directory
```

The runner verifies all pinned input hashes, then builds separate Neuro C++17
and Bitmap C++20 executables in GCC Debug, Release with NDEBUG, and ASan+UBSan.
Checks remain active in Release. Each Neuro execution requires 316 checks;
each Bitmap execution requires 556: 2616 checks across six executions.
Exact commands, compiler/platform, exit codes and complete logs are saved.

A nonnumeric throwing Value starts in an asymmetric 2x3 matrix and attempts
resize(3,5). A literal independent operation schedule identifies every selected
default-constructor and copy-assignment site. Each case throws before mutation
at one site and requires the exact marker and both operation counts.

* Neuro processes rows in order, moving each original row before resizing it.
  Failures in row 0 destroy that row but preserve row 1, including its address
  and all three sentinel values. Later failures destroy both original rows.
  The schedule covers 15 constructor and six assignment failures (21 cases).
* Bitmap moves both original rows before constructing the new third row or
  resizing any row. Every selected failure leaves both original rows moved out.
  Its new-row initialization constructs five elements plus a temporary Value,
  then assigns that default value to five elements. It subsequently reallocates
  all three rows, including the already matching third row. The schedule covers
  21 constructor and 16 assignment failures (37 cases).

Both retain original matrix rows/cols/size/capacity metadata after these failures,
but individual rows can be empty. The harness reads row metadata and only accesses
storage belonging to a surviving nonempty row; it never indexes a moved-out row
using the stale matrix column count. Live counts require all temporary and moved
storage to be destroyed. With faults disabled, resize succeeds again, restores a
consistent 3x5 matrix, and preserves only surviving original values. Destroying
all matrices returns the live count to zero. This diagnoses data loss and stale
shape relationships; neither is a required shared-package compatibility behavior.
A safe resize/error contract or adapter must be reviewed before extraction.

Scope is this exact growth shape and throw-before-mutation Value, with no OpenMP
execution. It does not qualify actual Pixel, bad_alloc, shrink/equal-size/zero
matrix cases, arbitrary element side effects, throwing destructors, Windows/MSVC,
strict Neuro C++20, full characterization, extraction, a package version or
consumer migration. Existing normal-storage and row-resize evidence retains its
own scope. Independent fixture and exception-contract review remains required.
Rollback removes only this three-file local slice after preserving run evidence;
legacy inputs, consumer owners and prior characterization remain unchanged.
