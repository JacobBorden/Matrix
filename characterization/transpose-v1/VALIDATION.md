# November 2, 2026 validation

Canonical fetched origin/master: 1500576c382036085bcd5b267c815e9b9117d1ba.
Persistent branch: codex/matrix-legacy-characterization; preceding head 34193ae62d2c588842bf2d244d941d137015739d. Canonical base is already an ancestor; no integration conflicts.

Command: `python3 characterization/transpose-v1/run.py --output "$PWD/work/20261102-transpose"`, exit 0. GCC (SUSE Linux) 15.2.1 20251006, Linux/libstdc++; Neuro C++17 and Bitmap C++20. Six separate Debug/Release(NDEBUG)/ASan+UBSan executions each passed 99 checks, total 594, zero failures or sanitizer diagnostics. Exact compiler commands, pinned source hashes and logs remain in local work/20261102-transpose/results.json. All seven input hashes verified. Scoped git diff --check passed. No production/consumer source changed; unchanged slices retain prior evidence. The commit containing this note identifies the validated source tree.

Empty transpose dimensions differ and require independent adapter review. Actual Pixel, failures during transpose, serialization, Windows/MSVC, full characterization, extraction, exact-package and downstream qualification remain open. Remote current-head CI is reported separately and is not implied by local success.
