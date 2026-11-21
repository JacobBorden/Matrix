# Legacy transpose characterization v1

Acceptance: isolated pinned Neuro C++17 and Bitmap C++20 executables each pass 99 checks in Debug, Release (NDEBUG) and ASan+UBSan, with immutable input hash verification. Const return type is checked at compile time.

Six asymmetric rectangular, row/column vector, scalar and square shapes check exact independent golden mapping, source preservation, double transpose and independent returned storage for int/float/double. Three empty shapes expose a contract difference: Neuro swaps empty dimensions; Bitmap normalizes to 0x0. This requires an explicit adapter decision before extraction.

Run `python3 characterization/transpose-v1/run.py --output <absolute-evidence-directory>` from this worktree. No legacy definitions are linked together. Existing multiplication source manifest pins remain unchanged. This slice does not qualify Pixel, element exceptions, saved models, Windows/MSVC, shared package or consumer migration.
