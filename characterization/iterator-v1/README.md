# Legacy iterator characterization

Separate executables use the seven hash-verified inputs in ../multiplication-v1/sources.json (Neuro 87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c; Bitmap 4df8a9cecdce77c9e0eb98b0b0897a0c473c0710).

Acceptance: asymmetric 2x3 row-major golden traversal for int/float/double; row and element prefix/postfix increment/decrement, prefix return identity, storage aliases and default empty range equality, with checks active under NDEBUG. Neuro lacks const Matrix begin/end; Bitmap exposes mutable rows through const Matrix iterators. Both expose mutable elements through const row iterators. This records legacy behavior, not an approved const contract.

Run `python3 characterization/iterator-v1/run.py --output /absolute/fresh/output` from this worktree. Neuro uses C++17; Bitmap C++20. The runner verifies pins and runs Debug, Release and ASan+UBSan separately. It saves exact commands, logs and counts. Expected counts: Neuro 129, Bitmap 96 per execution.

No column iterator traversal is attempted: independently allocated rows do not establish contiguous scalar storage. Full standard iterator conformance, invalidation, degenerate shapes beyond default empty, Pixel, Windows/MSVC, extraction and package/consumer acceptance remain open. No legacy input, kernel, package header or consumer is modified.

November 8 increment: row-element addition/subtraction (including negative offsets), compound return identity, signed distances, all ordering operators, const iterator signed indexing/storage mutation and std::distance, each for int/float/double. Arithmetic and comparisons use only a single nonempty allocated row and its one-past endpoint. No cross-row arithmetic or empty-pointer subtraction is performed. These operations do not establish full standard iterator conformance.

November 9 increment: Neuro matrix-row arithmetic within the allocated two-row array, signed distances, negative offsets, compound identities, ordering, std::distance and aliases. Bitmap MatrixIterator has no arithmetic operators or iterator traits; its existing traversal remains covered. Neuro MatrixIterator indexing is not exercised: its const member calls a non-const dereference. No full iterator conformance or indexing acceptance is claimed.
