# November 1, 2026 validation

Canonical base: `1500576c382036085bcd5b267c815e9b9117d1ba` (fresh origin/master fetch). Persistent branch: `codex/matrix-legacy-characterization`, prior head `d345230f2f7928328765cd74de999fe181437fd4`; base already ancestor, no integration conflicts. Canonical checkout user changes were preserved.

Command from the persistent worktree: `python3 characterization/empty-shapes-v1/run.py --output /mnt/c/Users/jacob/Documents/Codex/2026-10-31/task/four-publication-repair/publication/worktrees/matrix/work/20261101-empty-shapes`, exit 0. GCC (SUSE Linux) 15.2.1 20251006, Linux/libstdc++, Neuro C++17 and Bitmap C++20. Six independent Debug/Release/ASan+UBSan executions each passed 81 checks: 486 total, zero failures or sanitizer diagnostics. Seven pinned input hashes verified by the runner. Full commands, output and counts remain in the local output results.json and logs; generated artifacts are excluded from Git.

The unchanged eleven slices retain prior exact-source validation and green hosted CI at d345230; this increment changes none of their source or runners. CI now adds empty-shapes as a twelfth job. Independent contract review and Windows/MSVC remain open. Next useful step: qualify these empty shape adapter decisions alongside remaining Pixel and serialization characterization before extraction. No package or downstream acceptance is claimed.
