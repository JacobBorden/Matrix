# Empty shape characterization v1

Run `python3 characterization/empty-shapes-v1/run.py --output /absolute/fresh/output` from this task worktree. The runner verifies the unchanged seven-file legacy manifest and builds separate Neuro C++17 and Bitmap C++20 executables in Debug, Release (NDEBUG), and ASan+UBSan modes. Each executes 81 active checks.

The 0x0, 0x3 and 2x0 fixtures cover int, float and double construction, deep copy, assignment over populated storage, move/source reset, safely observable empty rows, growth default initialization, source independence, resize back to an empty shape and regrowth after data discard. Both constructors retain dimensions even when one dimension is zero. Neuro resize normalizes every zero dimension to 0x0; Bitmap resize retains dimensions. These observed differences need explicit adapter decisions before extraction. No empty element access or invalid iterator arithmetic is performed.

This slice does not qualify negative dimensions, overflow, allocation failure, actual Pixel, serialization, a unified API, Windows/MSVC, extraction or consumer migration. Immutable legacy pins remain in multiplication-v1/sources.json. Full contract and independent acceptance remain open.
