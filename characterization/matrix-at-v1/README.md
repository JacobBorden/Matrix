# Valid Bitmap matrix.at characterization

Run `python3 characterization/matrix-at-v1/run.py --output <absolute-evidence-directory>` from the canonical Matrix workspace.

This bounded follow-up closes the direct valid-access gap recorded by the October 16 bounds steward. The runner verifies the unchanged multiplication-v1 source manifest, then builds the pinned Bitmap header alone under C++20 with GCC in Debug, Release (NDEBUG), and ASan+UBSan configurations. Each executes 39 active checks (13 per float/double/int): exact mutable reference and const value return types; all six asymmetric values and mutable storage addresses; independent const snapshots; mutation visibility in both directions; direct assignment; neighbor and shape preservation. Integer-valued fixtures are exactly representable in all tested types. No production or dependency source changes.

Expected contract observation: mutable matrix.at aliases the selected element; const matrix.at returns an independent value. Neuro has no corresponding two-index matrix.at API (covered by bounds-v1); no Neuro executable or unified API acceptance is claimed here. Legacy bounds/errors remain covered by the unchanged bounds-v1 slice. References are used only while storage remains stable; validity across resize, move or destruction is not asserted.

Independent fixture/adapter review remains required. Actual Pixel, failure paths, serialization, full characterization, strict Neuro C++20, Windows/MSVC, exact-candidate CI, extraction, package-version and consumer acceptance remain open. No downstream migration or publication is authorized by these local results.
