# Bitmap Pixel storage characterization

This fixture includes Bitmap's public header from revision
`4df8a9cecdce77c9e0eb98b0b0897a0c473c0710` and instantiates its
`Matrix::Matrix<Pixel>` in a Bitmap-only executable. The runner verifies the
public bitmap header, its bitmap-file header dependency, and the existing Bitmap
Matrix header byte for byte against [sources.json](sources.json).

The fixture checks all four BGRA channels, including zero and 255. It verifies
zero initialization, mutable aliases, independent const reads, checked bounds,
independent row storage on copy, and row storage transfer on move. It also
checks that a moved-from matrix can be reused, that growing or shrinking
preserves the overlapping pixels, and that the fill overload initializes newly
exposed pixels. Each build runs 21 checks, including in Release with NDEBUG.

```sh
python3 characterization/pixel-v1/run.py --output /tmp/matrix-pixel
```

The runner saves the exact commands, toolchain and platform details, exit codes,
test counts, and diagnostics in the output directory. The November 24, 2026
local run used GCC 15.2.1 (SUSE build 20251006), Linux x86-64 and libstdc++ in
C++20 mode. Debug, Release and ASan+UBSan each passed 21 checks with zero
failures or runtime diagnostics. Matrix base was
`1500576c382036085bcd5b267c815e9b9117d1ba`.

The exercise covers ordinary storage operations with the public Pixel type.
Bitmap file conversion, serialization, failure injection, Windows/MSVC and
compatibility with an extracted Matrix package require separate validation.
