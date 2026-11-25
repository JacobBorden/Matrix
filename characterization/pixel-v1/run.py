#!/usr/bin/env python3
"""Run the pinned Bitmap Pixel fixture in three isolated executables."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess

SLICE = Path(__file__).resolve().parent
ROOT = SLICE.parents[1]
LEGACY = SLICE.parent / "multiplication-v1" / "legacy" / "bitmap"
EXPECTED_CHECKS = 28


def run(output, results, label, command):
    env = {**os.environ, "ASAN_OPTIONS": "detect_leaks=1:halt_on_error=1",
           "UBSAN_OPTIONS": "halt_on_error=1"}
    log = output / f"{label}.log"
    try:
        completed = subprocess.run(
            command, cwd=ROOT, capture_output=True, timeout=90, env=env
        )
        log.write_bytes(completed.stdout + completed.stderr)
        result = dict(label=label, command=command, exit_code=completed.returncode,
                      log=str(log))
        if label.endswith("-tests"):
            counts = re.fullmatch(rb"RESULT passed=(\d+) failed=(\d+)\n", completed.stdout)
            result.update(passed=int(counts[1]) if counts else 0,
                          failed=int(counts[2]) if counts else None,
                          stderr_bytes=len(completed.stderr))
    except subprocess.TimeoutExpired as exc:
        log.write_bytes((exc.stdout or b"") + (exc.stderr or b"") + b"\nTimed out\n")
        result = dict(label=label, command=command, exit_code=124, log=str(log))
    results.append(result)
    print(label, result["exit_code"], result.get("passed", ""), flush=True)
    return result["exit_code"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    output = parser.parse_args().output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    sources = json.loads((SLICE / "sources.json").read_text())
    for path, expected in sources["sha256"].items():
        if hashlib.sha256((LEGACY / path).read_bytes()).hexdigest() != expected:
            raise SystemExit(f"Source hash mismatch: {path}")

    results = []
    run(output, results, "compiler", ["g++", "--version"])
    run(output, results, "platform", ["uname", "-a"])
    modes = {
        "debug": ["-O0", "-g"],
        "release": ["-O2", "-DNDEBUG"],
        "sanitized": ["-O1", "-g", "-fsanitize=address,undefined",
                      "-fno-sanitize-recover=all", "-fno-omit-frame-pointer"],
    }
    for mode, flags in modes.items():
        executable = output / f"bitmap-{mode}"
        command = ["g++", "-std=c++20", "-Wall", "-Wextra", *flags,
                   "-I", str(LEGACY / "src"), str(SLICE / "pixel.cpp"),
                   str(LEGACY / "src/bitmap/bitmap.cpp"),
                   str(LEGACY / "src/bitmapfile/bitmap_file.cpp"), "-o", str(executable)]
        if run(output, results, f"bitmap-{mode}-build", command) == 0:
            run(output, results, f"bitmap-{mode}-tests", [str(executable)])

    evidence = dict(sources=sources, expected_checks=EXPECTED_CHECKS, results=results)
    (output / "results.json").write_text(json.dumps(evidence, indent=2) + "\n")
    tests = [result for result in results if "passed" in result]
    return int(
        any(result["exit_code"] != 0 for result in results)
        or len(tests) != len(modes)
        or any((result["passed"], result["failed"], result["stderr_bytes"])
               != (EXPECTED_CHECKS, 0, 0) for result in tests)
    )


if __name__ == "__main__":
    raise SystemExit(main())
