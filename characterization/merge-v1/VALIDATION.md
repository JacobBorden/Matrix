# November 5, 2026 bounded merge increment

Canonical origin/master fetched from JacobBorden/Matrix: `1500576c382036085bcd5b267c815e9b9117d1ba`, already an ancestor of task parent `a876a7cc8954cddd6184e14360d4a99c8b13348d`. Persistent branch: `codex/matrix-legacy-characterization`; existing draft PR #1 targets master. Primary checkout edits and task work/ were preserved.

Toolchain: g++ (SUSE Linux) 15.2.1 20251006, Linux/libstdc++; Neuro C++17 and Bitmap C++20. Exact legacy revisions remain Neuro `87ee5b2e65fe0a7bb4d1cab1de8feda64b07832c` and Bitmap `4df8a9cecdce77c9e0eb98b0b0897a0c473c0710`; runner verifies all seven source file hashes.

Validation command from this worktree: `python3 characterization/merge-v1/run.py --output "$PWD/work/20261105-merge-final"`. Exit 0: six executions × 69 checks = 414 passes, zero failed checks or sanitizer diagnostics. Exact commands and raw logs are recorded in that directory's results.json. Initial work/20261105-merge retains six compilation failures caused by helper template deduction between float/double matrices and integer golden literals. The helper now takes independent integer golden values (all exactly representable), resolving the fixture failure without legacy edits.

Scope and remaining qualification limits are in README.md. Current-head CI is separate from local evidence. Independent review, remaining characterization, strict Neuro C++20, supported-platform CI, extraction, exact package qualification and consumer acceptance remain open. No merge/release or consumer changes.
