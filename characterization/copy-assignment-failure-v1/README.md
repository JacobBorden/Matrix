# Legacy matrix copy assignment exceptions

This isolated characterization uses the unchanged, hash-verified Neuro master
`87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c` and Bitmap development
`4df8a9cecdce77c9e0eb98b0b0897a0c473c0710` inputs in
`../multiplication-v1/sources.json`. Each executable includes one legacy definition.
No consumer or legacy implementation is modified.

Run from the Matrix checkout:

```sh
python3 characterization/copy-assignment-failure-v1/run.py --output /absolute/evidence/directory
```

The runner checks every pinned input hash, then builds Neuro C++17 and Bitmap
C++20 separately in GCC Debug, Release with NDEBUG, and ASan+UBSan modes.
Each of six executions must pass exactly 197 active checks (1182 total).
Commands, compiler/platform identity, exit codes and full output are saved.

A nonnumeric throwing Value copies a 2x3 source into a distinct 3x2 destination.
Six default-construction sites and six copy-assignment sites are injected one at
a time; each case requires the exact exception marker and both operation counts.
The source retains shape, all asymmetric sentinels, capacity and storage identity.

* Neuro preserves the original destination shape, values, capacity and address;
  all temporary elements are destroyed. Its matrix copy-and-swap and row
  temporary-storage assignment supply the observed strong guarantee.
* Bitmap changes matrix metadata to the source shape, destroys old destination
  storage and commits rows incrementally. Completed rows contain source values;
  on copy-assignment failure the current allocated row contains the copied prefix
  and default suffix. On default-construction failure its current row reports
  size 3 while its storage allocation never committed. The harness reads only
  that row's metadata and **never indexes its missing storage**. Later rows
  remain empty. Live counts distinguish partial rows from cleaned temporaries.

Both implementations successfully retry with faults disabled, returning the
destination reference and an independent full copy. Scope destruction restores
zero live elements. Armed self-assignment performs zero element operations,
returns the same matrix and preserves its state. These are snapshot observations;
Bitmap's partial/inconsistent post-failure state needs an explicit error/lifetime
contract or repaired adapter before a shared implementation is accepted.

This qualifies only selected element default-constructor/copy-assignment failures
with a throw-before-mutation Value. It does not qualify actual Pixel, bad_alloc,
throwing destructors, arbitrary element side effects, matrix resize, other shapes,
Windows/MSVC, strict Neuro C++20, a versioned shared package or downstream use.
No undefined invalid-row access is used as a compatibility requirement. A future
repair must be checked against these diagnostic observations and a reviewed safe
contract. Independent fixture/guarantee review remains required.

Rollback removes only this three-file uncommitted slice after preserving its
run evidence; existing characterization and source inputs remain intact.
