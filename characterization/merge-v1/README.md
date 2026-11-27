# Merge characterization

Run `python3 characterization/merge-v1/run.py --output /absolute/output` from this persistent task worktree. All seven legacy source hashes are verified against the existing multiplication manifest. Neuro C++17 and Bitmap C++20 compile separately in Debug, Release/NDEBUG and ASan+UBSan modes. Checks stay active with NDEBUG: 69 per executable, 414 total.

Integer, float and double fixtures use independently enumerated exact golden values for asymmetric horizontal/vertical concatenation, operand order and self-merges. They check const/value signatures, unchanged inputs, mutation independence in both directions, invalid_argument for incompatible nonempty shapes, source preservation after rejection, default empty identities and matching 0x2/2x0 identities. Both pinned implementations agree on these selected cases.

This bounded slice does not qualify mismatched degenerate dimensions, generic Pixel, allocation/element exceptions, overflow, serialization, OpenMP or Windows/MSVC. In particular, matching empty identities do not establish safety for arbitrary empty-shaped merges. Shared API/error adapters and independent review remain required before extraction. Kernels, source pins and consumers are unchanged.
