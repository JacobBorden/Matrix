# Scalar multiplication characterization

Run `python3 characterization/scalar-multiply-v1/run.py --output /absolute/output` from this task worktree. The runner hash-verifies all seven existing pinned Neuro/Bitmap inputs, compiles them separately and checks Debug, Release/NDEBUG and ASan+UBSan. Neuro uses C++17; Bitmap uses C++20.

Checks remain active with NDEBUG: int/float/double rectangular golden values for zero, identity and negative scalars; const source preservation; independent results; compound reference identity; aliased scalar behavior; 0x0, 0x3 and 2x0 shapes; exactly representable fractional products for float/double. Neuro's left scalar overload is also tested. Expected counts are Neuro 64 and Bitmap 61 per execution (375 total).

For initial values [-4,2,7;0,-3,9], `m *= m[0][0]` yields [16,-8,-28;0,12,-36] in Neuro (scalar copied) and [16,32,112;0,-48,144] in Bitmap (scalar reference observes the first update). Both retain empty dimensions. These are observations, not an accepted shared alias contract. Values avoid integer overflow; nonfinite, overflow, generic Pixel, OpenMP and other platform behavior are outside this slice. No kernels, consumers or pins change.
