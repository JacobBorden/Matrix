# Legacy split and Sigmoid characterization v1

Run from the canonical Matrix checkout:

```sh
python3 characterization/api-v1/run.py --output /absolute/Matrix/work/RUN/validation
```

The runner verifies the complete unchanged multiplication-v1 source manifest and compiles Neuro and Bitmap into separate executables. Neuro remains C++17 because strict C++20 constructor syntax was rejected in the independently reviewed multiplication slice; Bitmap uses C++20. GCC Debug, Release with NDEBUG, and ASan+UBSan each run 67 checks: 13 split checks for float/double/int and four Sigmoid checks each for float/double, plus ten checker rejection regressions for each floating type. Checks are active in Release.

Splits use a 4x6 asymmetric coordinate fixture, with vertical/horizontal values and copy independence, default halves, zero/nondividing counts and odd default shapes. Default 0x0 splitting returns no parts in Neuro; explicit three-way splitting returns three 0x0 matrices. Bitmap rejects both empty calls. Both reject a zero split count even for empty input.

Sigmoid uses six independent precomputed mathematical references from -2 through 3 with absolute tolerance 1e-6. Compile-time return-type assertions and mutation/alias probes distinguish Neuro's const value-returning copy from Bitmap's mutating reference. Both empty Sigmoid contracts are checked. The checker explicitly rejects nonfinite values. Negative regressions inject NaN at each of the six positions and into the entire matrix, positive/negative infinity, and a finite out-of-tolerance value for float and double.

These are observed legacy differences for reviewer/adapter decisions, not an accepted unified API. No integer Sigmoid, noncanonical zero shapes, signed dimension/index conversion, unchecked Bitmap column indexing, exception/allocation failures, Pixel, determinants or serialization are qualified here. No extraction, package acceptance, downstream migration or supported-platform CI is established. Consumer owners retain their work.
