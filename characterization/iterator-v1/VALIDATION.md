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

# November 9, 2026 matrix-row arithmetic increment

Canonical fetch origin/master: 1500576c382036085bcd5b267c815e9b9117d1ba, already ancestor of persistent codex/matrix-legacy-characterization. Parent task head: c95e5de214f24471cbb8b6ecdbc148e54ea70b01. Candidate is the commit containing this note; documentation does not change tested inputs.

Command in the persistent worktree above: `python3 characterization/iterator-v1/run.py --output "$PWD/work/20261109-matrix-row-arithmetic"`. GCC SUSE 15.2.1 20251006, Linux/libstdc++; Neuro C++17, Bitmap C++20. Seven source hashes verified; isolated Debug, Release/NDEBUG and ASan+UBSan executions. Exact commands/logs/toolchain are retained in that output directory.

Scope: Neuro arithmetic over its allocated row-object array, signed distances, negative offsets, compound identities, ordering, std::distance and mutation aliases for int/float/double. Bitmap traversal coverage retained; its matrix iterator lacks arithmetic/traits. Neuro matrix iterator indexing is not qualified (const indexing calls non-const dereference in the pinned source). No empty-pointer subtraction or cross-allocation scalar arithmetic. Full iterator conformance, Pixel, serialization, platform and exact package/consumer acceptance remain open. Primary unrelated edits and work/ preserved. Prior-head CI: all 34 checks succeeded; candidate CI is separate. Existing draft PR1 reused, no merge/release.

Result: exit 0, six executions: Neuro 129 checks each and Bitmap 96 each, 675 total passes, zero failures or sanitizer diagnostics.

# November 10, 2026 narrow-shape iterator increment

Canonical origin/master fetched in the canonical checkout: 1500576c382036085bcd5b267c815e9b9117d1ba, already ancestor. Persistent branch codex/matrix-legacy-characterization; parent c3c4a56087edcf954c0571aabfc71d2a5f56f9cd. Candidate is the commit containing this note; tested source hashes below bind the unchanged executable inputs.

Command in the persistent worktree recorded above: `python3 characterization/iterator-v1/run.py --output "$PWD/work/20261110-narrow-iterators"`. Exit 0: six executions, Neuro174/Bitmap141 each, 945 checks passed, zero failures or sanitizer diagnostics. GCC SUSE15.2.1 20251006 Linux/libstdc++; Neuro C++17, Bitmap C++20; Debug, Release/NDEBUG and ASan+UBSan. Seven immutable legacy source hashes verified. Exact commands and logs: work/20261110-narrow-iterators/results.json.

Acceptance for this increment: 1x1, 1x3 and 3x1 golden forward/reverse traversal, endpoint storage aliases, mutation and restoration for int/float/double. No empty decrement or cross-allocation scalar arithmetic. Full iterator conformance/invalidation, Pixel, serialization, Windows/MSVC and exact package/consumer acceptance remain open. Existing workflow includes this slice; prior head has 34 successful checks, candidate CI is separate. Primary checkout changes and unrelated work/ preserved. Draft PR1 reused; no merge/release.

Tested inputs (documentation-only additions followed execution):
- iterator.cpp: `94384c8ec9cb9366cb34d9bd10152f573372125146e2079b0fd4e36cf90c9382`
- run.py: `840b675ae3ca138dc5b4f47abab51b7ccd244c705733e3b0715eb775fa7c412c`
