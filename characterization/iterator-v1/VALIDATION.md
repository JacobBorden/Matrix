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

# November 11, 2026 standard row algorithm increment

Canonical origin/master fetched in the designated canonical checkout: 1500576c382036085bcd5b267c815e9b9117d1ba, already ancestor. Persistent branch codex/matrix-legacy-characterization; parent 511965abaf3a3d991e6a6670208bae0f0f16dfc5. Candidate is the commit containing this note.

Command from /mnt/c/Users/jacob/Documents/Codex/2026-10-31/task/four-publication-repair/publication/worktrees/matrix: `python3 characterization/iterator-v1/run.py --output "$PWD/work/20261111-row-algorithms"`. Exit 0: six separate Debug, Release/NDEBUG and ASan+UBSan executions, Neuro195/Bitmap162 each, 1071 total passes, zero failures or sanitizer diagnostics. GCC SUSE15.2.1 20251006 Linux/libstdc++; Neuro C++17 and Bitmap C++20. All seven pinned legacy hashes verified; exact commands and raw logs saved in that output's results.json.

Acceptance: std::reverse and std::copy on one allocated row, golden values, returned copy endpoints, restoration after two reversals, storage aliases and untouched adjacent row, for int/float/double. All pass. Source review confirms no scalar arithmetic between row allocations. Full standard iterator conformance, Pixel, serialization, supported platforms and exact package/consumer qualification remain open. Primary unrelated edits and work preserved. Prior-head 34 CI checks succeeded; new candidate CI is separate. Existing draft PR1 retained; no merge/release or consumer changes.

Tested input SHA-256 (only documentation added after execution):
- iterator.cpp: `5b8144438c193b12bdb8618958591fced6d2266c74e9da695f0ecf4d2490bbe0`
- run.py: `11b67d047a524b0c67cbff1be93bca8554c4ec582843e5df43123325ea651d8d`

# November 12, 2026 within-row sort increment

Canonical fetch origin/master: 1500576c382036085bcd5b267c815e9b9117d1ba, already ancestor. Persistent branch codex/matrix-legacy-characterization; parent 8f048b0bab17b0a95ab9b0b810a55165c16e9f95. Candidate is the commit containing this note.

Working directory: /mnt/c/Users/jacob/Documents/Codex/2026-10-31/task/four-publication-repair/publication/worktrees/matrix

Command: `python3 characterization/iterator-v1/run.py --output "$PWD/work/20261112-row-sort"`. Exit 0: six isolated executions, Neuro207/Bitmap174 each, 1143 passes, zero failures or sanitizer diagnostics. GCC SUSE15.2.1 20251006 Linux/libstdc++; Neuro C++17 and Bitmap C++20; Debug, Release/NDEBUG, ASan+UBSan. Seven pinned input hashes verified. Exact commands/logs in that output's results.json.

Acceptance: std::sort within a five-element allocated row, ascending/descending independent golden values including duplicate, custom comparator, unchanged storage/endpoints and adjacent row for int/float/double. All passed. No cross-row scalar arithmetic. Full iterator conformance, Pixel, serialization, supported platforms and exact package/consumer acceptance remain open. Primary unrelated changes and work preserved. Existing draft PR1 reused; prior-head34 CI checks successful, candidate CI recorded separately. No merge/release or consumer edits.

Tested input SHA-256 (documentation-only addition after execution):
- iterator.cpp: `f1651b3f19a2289dc22024f8a897dde6c04fc1e2714be2d8f5726629e5ef7599`
- run.py: `2cfbc8bd72382e8d182afbd7f1fbe6eed04abbdee4c89d03c48d5599de3b95fd`

# November 13, 2026 move iterator lifetime increment

Canonical fetched origin/master: 1500576c382036085bcd5b267c815e9b9117d1ba (already ancestor). Persistent branch codex/matrix-legacy-characterization, parent f158a64a79c80f16632f695623ee9bbe06421fbb. Candidate is the commit containing this note.

Acceptance scope: int/float/double saved matrix-row and row-element iterators, including endpoints, retain transferred storage after move construction and move assignment into a nonempty destination. Golden values, storage aliases, mutation through saved iterators, assignment return identity and emptied sources are checked. Replaced-destination iterators are never used. Resize invalidation, full conformance, Pixel, serialization, supported platforms and package/consumer acceptance remain open.

Explicit working directory: /mnt/c/Users/jacob/Documents/Codex/2026-10-31/task/four-publication-repair/publication/worktrees/matrix.
Command: `python3 characterization/iterator-v1/run.py --output "$PWD/work/20261113-move-iterators"` exited 1 after the Bitmap ASan+UBSan compiler exceeded its 90-second timeout. Before that failure, Neuro Debug/Release/ASan+UBSan each passed231/231; Bitmap Debug/Release each passed198/198 (1089 total, no test failures). Those unchanged executions are reused. The runner aborts before writing results.json on timeout; individual completed command logs are retained. No timeout is reported as a pass.

Toolchain: GCC (SUSE Linux)15.2.1 20251006, Linux/libstdc++; Neuro C++17 and Bitmap C++20. Seven legacy source hashes verified by runner, pins unchanged. Focused retry uses the identical Bitmap sanitizer compile flags with a 180-second bound and a fresh executable bitmap-sanitized-retry; exact commands/exits in work/20261113-move-iterators/retry-results.json, raw logs adjacent.

SHA256 iterator.cpp: a81f5d8c19011d91af61efa29b2b0c1c5fbd7e60c5e777e5859d2ea4db640e4a

SHA256 run.py: 180dcabd6756e2996c4ea72a2bc0a30d0d4d7a2607ee05cb06ee3eef994be629

Focused retry completed: build exit0, test exit0,198/198 checks and no sanitizer diagnostics. Combined six executions:1287 passes,zero test failures. Initial runner exit1/compiler timeout remains preserved above. Cached diff whitespace check passed. CI for this candidate must be inspected separately; previous head34 successful checks do not validate this revision.
