# Legacy elementwise addition and subtraction

Run `python3 characterization/elementwise-v1/run.py --output /absolute/fresh/output` from this task worktree. Each pinned definition compiles separately, Neuro in C++17 and Bitmap in C++20. All seven existing input hashes are checked before compilation.

For each of int, float and double, 22 active checks cover independently tabulated asymmetric 2x3 matrix/scalar addition and subtraction, const operands, independent output storage, compound return identity and mutation, self aliasing, incompatible nonempty shapes and empty/nonempty operands. Values are small exactly representable integers, avoiding overflow and tolerance ambiguity. Debug, Release/NDEBUG and ASan+UBSan give six executions and 396 checks. These tests qualify this finite fixture domain, not all floating-point values, allocation failures or scalar references aliasing matrix elements.

Neuro rejects mismatched empty/nonempty dimensions for binary addition and subtraction. Bitmap binary operations treat empty operands as additive identities (or negate the nonempty right operand for empty-minus-nonempty). Both compound operators reject these mismatches and retain the original left operand. Exception messages are not qualified. Static checks record value return for const matrix addition and reference return for compound subtraction.

No legacy source is repaired. Empty-shape combinations beyond the default 0x0 operand, Pixel operations, serialization, Windows/MSVC, extraction and exact-version consumer qualification remain open. Independent review must choose an explicit empty-operand contract or adapter before consolidation.
