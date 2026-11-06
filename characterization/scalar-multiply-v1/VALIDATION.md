# November 4, 2026 validation

Fetched canonical origin/master: `1500576c382036085bcd5b267c815e9b9117d1ba`.
Persistent branch `codex/matrix-legacy-characterization`, parent `7fdca64b211c072cad15f6096210a511b145a07a`; base already an ancestor, local and remote starting tips identical. The commit containing this note identifies the tested candidate. Canonical checkout edits and task work/ preserved.

Command from this persistent worktree:

```
python3 characterization/scalar-multiply-v1/run.py --output "$PWD/work/20261104-scalar-multiply-final"
```

Exit 0, GCC (SUSE Linux) 15.2.1 20251006. Six isolated Debug, Release/NDEBUG and ASan+UBSan executions: Neuro C++17 64 checks each, Bitmap C++20 61 each; 375 passes, zero failed checks or sanitizer diagnostics. All seven pinned hashes verified. Exact commands, compiler/platform, logs and counts are retained in excluded `work/20261104-scalar-multiply-final/results.json` and neighboring logs. Other slices unchanged; no unchanged tests rerun. Prior published head had 28 successful hosted checks. Candidate CI is checked after publication.

Historical initial run at work/20261104-scalar-multiply exited 1 because expected totals were mistakenly 75/72; all six executables exited 0 at 66/63. Removed two placeholder integer-fraction checks, corrected expected totals, then reran the exact final fixture successfully. Initial evidence remains preserved.

Adds one CI slice. Scoped whitespace check passes. Independent fixture and scalar-alias adapter review remains next: Bitmap's reference scalar changes during aliased compound multiplication, Neuro copies the scalar. This does not accept a shared contract, actual Pixel, serialization, Windows/MSVC, extraction, package version or downstream migration. No consumer changes, merge or release.
