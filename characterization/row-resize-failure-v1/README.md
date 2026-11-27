# Legacy row resize: element-failure characterization

Run from the canonical Matrix workspace:

```sh
python3 characterization/row-resize-failure-v1/run.py --output /absolute/evidence/directory
```

The runner verifies every immutable input in `../multiplication-v1/sources.json`.
Neuro master `87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c` is compiled as C++17;
Bitmap development `4df8a9cecdce77c9e0eb98b0b0897a0c473c0710` as C++20.
Each header runs in a separate executable. Neuro's strict C++20 constructor-syntax
limitation remains open. No legacy or consumer source is changed.

A nonnumeric Value throws a distinct marker before incrementing its live count
at a selected default-construction site or before copying at a selected assignment
site. A three-element row contains asymmetric 17, -31, 59 sentinels. Grow to five
injects at all five new-element constructors and all three retained-element copy
assignments; shrink to two injects at both constructors and both assignments.
Twelve fault cases each check nine predicates: exact propagated marker, selected
site reached, unchanged original shape/values, capacity, storage identity,
cleanup of temporary elements, successful retry with prefix/default values,
retry live count and destruction balance. A final global balance check makes
109 active checks per execution, including Release with NDEBUG. Each runtime
must produce exactly 109 passes and zero failures. Six Debug/Release/ASan+UBSan
executions give 654 checks. A nonzero build/test or missing execution fails the
runner; logs and exact commands are retained.

Both pinned implementations prepare temporary row storage before committing
resize. These fixtures characterize preservation and cleanup for the selected
Value exceptions. They do not qualify matrix-wide resize, copy assignment,
allocator/bad_alloc failure, throwing destructors, equal-size resize, actual
Bitmap Pixel, arbitrary element side effects or supported-platform behavior.
They do not promise reference validity after successful resize. Original storage
identity is checked only after failure, before retry. No arithmetic API is
instantiated for Value. Unified error/lifetime contracts still require independent
review, followed by full characterization and package qualification. Linux GCC
results do not establish Windows/MSVC or exact-candidate CI acceptance. No shared
package version or downstream migration is accepted by this slice.
