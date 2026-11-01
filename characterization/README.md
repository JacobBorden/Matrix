# Pinned legacy Matrix characterization

These eleven completed slices record the behavior of two separately compiled legacy implementations before a shared package is extracted. They do not change the Matrix package header, either consumer repository, numerical kernels, or storage layout.

The immutable inputs and their SHA-256 hashes are listed in [multiplication-v1/sources.json](multiplication-v1/sources.json): Neuro-Net-DLL at `87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c` and Bitmap at `4df8a9cecdce77c9e0eb98b0b0897a0c473c0710`. The seven source/test support files are preserved byte-for-byte. Every runner verifies these hashes before compiling, and each implementation has its own executable.

Run a slice from this repository with a fresh absolute output directory:

```sh
python3 characterization/multiplication-v1/run.py --output /tmp/matrix-multiplication
```

The other slice names are `storage`, `api`, `determinant`, `division`, `inverse`, `bounds`, `matrix-at`, `row-resize-failure`, `copy-assignment-failure`, and `matrix-resize-failure`. The focused CI workflow runs all eleven independently. Runners retain exact compiler commands, test counts, exit codes, and diagnostics in their output directories. Expected compiler failures, rejected mutation controls, and qualified integer-division/singular-matrix diagnostics remain explicit negative controls; they are not repaired behavior or passing feature acceptance.

Neuro is characterized in C++17 because its pinned constructor syntax is rejected by strict C++20. Bitmap uses C++20. Debug, Release, and sanitizer results have the limits described by each slice; direct `matrix.at` coverage is Bitmap-only. Historical API evidence contains 402 passing checks across six executions plus two rejected mutation controls.

Scoped branch publication and this draft PR were authorized on November 1, 2026. Older slice notes describe historical local runs. Their `work/` paths identify local evidence excluded from the repository; use the runner commands to produce fresh evidence. Publication does not qualify a shared API, extracted package, package version, actual Pixel type, Windows/MSVC, or downstream migration. Those technical acceptance requirements remain open.
