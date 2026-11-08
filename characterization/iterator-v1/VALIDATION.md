# November 7, 2026 validation

Canonical fetch: `git fetch origin master` in the canonical Matrix checkout. Base origin/master: 1500576c382036085bcd5b267c815e9b9117d1ba, already ancestor of persistent codex/matrix-legacy-characterization. Parent task head: a47967d9743cf3a644a80729910ea66339c190b8. Candidate is the commit containing this note; source and workflow were tested before this documentation-only addition.

Working directory: /mnt/c/Users/jacob/Documents/Codex/2026-10-31/task/four-publication-repair/publication/worktrees/matrix

Command: `python3 characterization/iterator-v1/run.py --output "$PWD/work/20261107-iterator"`

Exit 0. Six executions: Neuro 45 checks each, Bitmap 51 each, 288 total passed, zero failed checks or sanitizer diagnostics. GCC SUSE 15.2.1 Linux; Neuro C++17, Bitmap C++20; Debug, Release with NDEBUG, ASan+UBSan. Seven immutable source hashes verified. Detailed commands/toolchain/logs: work/20261107-iterator/results.json (local generated evidence, not committed).

Independent fixture/const adapter review and broader platform/package gates remain open. Prior 32 checks on parent head were successful; candidate remote CI is separate. Primary checkout edits and untracked work preserved. No merge/release or consumer migration.

# November 8, 2026 row-element arithmetic increment

Canonical origin/master fetched in the canonical checkout: 1500576c382036085bcd5b267c815e9b9117d1ba; already ancestor, so no integration merge needed. Persistent branch: codex/matrix-legacy-characterization. Parent: 308eeed2160b696ad410528530c54f3e887b1b00. Candidate is the commit containing this note.

Command from the same persistent worktree recorded above: `python3 characterization/iterator-v1/run.py --output "$PWD/work/20261108-iterator-arithmetic"`. GCC (SUSE Linux) 15.2.1 20251006, Linux/libstdc++; Neuro C++17 and Bitmap C++20. Debug, Release/NDEBUG and ASan+UBSan in isolated executables; seven immutable source hashes checked. Exit 0: 90 Neuro and 96 Bitmap checks per execution (558 total), zero failed checks or sanitizer diagnostics; exact commands, logs and results remain in work/20261108-iterator-arithmetic/results.json.

Scope: signed within-row arithmetic/distance, ordering, indexed aliases, compound identities and std::distance for int/float/double. No kernels, pinned inputs or consumers changed. Column traversal and full standard iterator conformance remain unqualified; Pixel, serialization, supported-platform and exact package/consumer acceptance remain open. Existing const alias behavior is observed, not approved. Primary checkout changes and untracked work preserved. Existing draft PR1 reused; remote CI status is recorded separately after publication.
