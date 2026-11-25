# Bitmap Pixel storage characterization

This fixture includes Bitmap's public header from revision
`4df8a9cecdce77c9e0eb98b0b0897a0c473c0710` and instantiates its
`Matrix::Matrix<Pixel>` in a Bitmap-only executable. The runner verifies the seven pinned
source inputs byte for byte against [sources.json](sources.json). It links the
original bitmap and bitmap-file translation units alongside the fixture.

The fixture checks all four BGRA channels, including zero and 255. It verifies
zero initialization, mutable aliases, independent const reads, checked bounds,
independent row storage on copy, and row storage transfer on move. It also
checks that a moved-from matrix can be reused, that growing or shrinking
preserves the overlapping pixels, and that the fill overload initializes newly
exposed pixels. The storage checks account for 21 of the 28 checks in each
build, including in Release with NDEBUG.

The conversion fixture uses asymmetric 2×3 images and three sentinel padding
bytes per 24-bit row. Positive and negative heights both produce matrix rows in
physical buffer order. This means two input buffers for the same logical image
yield opposite matrix row orders. The conversion skips padding, sets alpha to
255 for 24-bit input, and yields an empty matrix when even the last padding byte
is missing. Export writes a positive-height, uncompressed 32-bit header and
exact row-ordered BGRA bytes. Reading those bytes back preserves all four
channels, including non-opaque alpha values.

These are observations of the pinned in-memory conversion routines. On-disk
BMP I/O, tensor export orientation, SIMD paths, Windows/MSVC, and an extracted
Matrix package require separate validation.

```sh
python3 characterization/pixel-v1/run.py --output /tmp/matrix-pixel
```

The runner saves the exact commands, toolchain and platform details, exit codes,
test counts, and diagnostics in the output directory. The November 24, 2026
local run used GCC 15.2.1 (SUSE build 20251006), Linux x86-64 and libstdc++ in
C++20 mode. Debug, Release and ASan+UBSan each passed 21 checks with zero
failures or runtime diagnostics. Matrix base was
`1500576c382036085bcd5b267c815e9b9117d1ba`.

The November 25, 2026 run used the same Matrix base and toolchain. Debug,
Release and ASan+UBSan each passed 28 checks, including seven conversion checks,
with zero failures or runtime diagnostics. The runner's local results.json
has SHA-256 `9f5f3bd5b652ee27eeca3c1bfbbb08d565c0a3c71c9c0498e855fd7fb44e2917`.
