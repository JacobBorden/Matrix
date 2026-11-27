# Legacy scalar division characterization v1

Run from the canonical Matrix checkout:

```sh
python3 characterization/division-v1/run.py --output /absolute/Matrix/work/RUN/validation
```

The runner verifies every immutable multiplication-v1 source hash. Separate executables include Neuro master 87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c (C++17) or Bitmap development 4df8a9cecdce77c9e0eb98b0b0897a0c473c0710 (C++20). GCC Debug, Release with NDEBUG, and ASan+UBSan execute 71 checks per binary, 426 total. Checks remain active in Release.

A signed asymmetric 2x3 fixture covers float/double/int copy shape/values, const input preservation, copy independence, assignment return alias and values, negative divisors, and default-empty behavior. Integer quotient expectations follow truncation toward zero. Floating checks cover signed zero and signed 1e-12 rejection, preservation after errors, and nextafter values below/at/above T(1e-9) for both operators. Neuro compares against T(1e-9); Bitmap promotes to double 1e-9. On this GCC platform float(1e-9) is rounded below the double threshold: Neuro accepts it, Bitmap rejects it. Boundary observations require a reviewed adapter decision.

Both implementations reject floating near-zero divisors before checking emptiness. Bitmap rejects integer zero, including empty matrices. Neuro's T(1e-9) becomes integer zero, so its `<` check misses zero: empty operations silently succeed, while nonempty integer `/` and `/=` reach undefined division. Four separate O0 UBSan-only diagnostic child processes examine nonempty integer zero. Neuro must exit nonzero with `runtime error: division by zero`; Bitmap must throw runtime_error, exit zero and emit INTEGER_ZERO_REJECTED without stderr. The runner counts these as qualified observations, not successful arithmetic or a required compatibility contract. No unsanitized integer-zero execution occurs. The runtime divisor is read from volatile storage into a normal int to support Bitmap's const-reference API.

The initial harness failed Bitmap compilation because it passed a volatile int to const int&; failed logs and first results remain in the run's validation directory. The corrected source is validated separately in validation-final. Legacy sources were never patched.

Review must decide the explicit zero-error contract; preserve ordinary value/copy/alias behavior while repairing the integer threshold during an authorized extraction/adapter task. This slice does not qualify INT_MIN abs/INT_MIN divided by -1, NaN/infinity inputs, nondefault zero shapes, inverse, actual Pixel, allocation failures, saved-model serialization, Windows/MSVC, remote CI, strict Neuro C++20 or a shared package version. Neuro/Bitmap owners retain migration. This is bounded local characterization pending independent review.
