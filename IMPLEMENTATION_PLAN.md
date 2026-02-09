# PSMBasedQuantificationTIMs Performance Implementation Plan (P0/P1)

## 1. Objective

Deliver measurable runtime improvements in `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs` without changing scientific correctness.

1. P0 target: at least 35% end-to-end runtime reduction on reference dataset.
2. P1 target: additional 15-30% improvement over P0.
3. Output requirement: numerical equivalence within defined tolerance.

## 1.0 Latest Status Snapshot (Updated 2026-02-08, 17:30)

Completed high-impact code changes now on branch:

1. XIC hot path m/z-first + allocation reductions in 3D profile lookup (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:30`).
2. Isotopic compare path avoids full-spectrum unzip/rezip materialization (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:639`).
3. Streamed SQL RT-index builder is default with legacy fallback (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1000`).
4. New direct RT/intensity extraction path removes `Peak2D` object creation in XIC pipeline (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:65`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:380`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1139`).
5. Runtime tuning is available via CLI parameters (`src/PSMBasedQuantificationTIMs/CLIArgumentParsing.fs:48`, `src/PSMBasedQuantificationTIMs/Program.fs:30`).
6. Environment-variable runtime tuning fallback removed from quantification code path (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`).

Latest full-run benchmark results (same dataset/params, N15 labeling):

| Run | Cache Peak Cap | Setup ms | QuantLoop ms @49300 | Total ms | Max WS MB | Max Heap MB | Notes |
|---|---:|---:|---:|---:|---:|---:|---|
| `full_latest_noCap_noenv_20260208_125529` | unbounded | 61781 | 2845563 | 2909738 | 23826 | 17469 | Fastest validated full run; highest RAM |
| `full_latest_80M_noenv_keepawake_20260208_161940` | 80M | 63513 | 3546352 | 3612353 | 15838 | 13348 | Balanced profile; ~7.99 GB lower max WS than no-cap |
| `full_latest_80M_noenv_20260208_115053` | 80M | 63826 | 3361870 | 3427782 | 16670 | 15344 | Earlier balanced full run (same code family) |
| `full_latest_80M_noenv_20260208_134422` | 80M | 64475 | 9110743 | 9177546 | 16400 | 14746 | Invalid for speed comparison: host sleep/pause artifact between 1100 and 1200 peptides |

Current recommendation from measured data:

1. Performance-first profile: `--runtime-peak-cache-max 700` and no peak-count cap.
2. Balanced profile: `--runtime-peak-cache-max 700 --runtime-peak-cache-max-peaks 80000000`.
3. Measured tradeoff (`no-cap` vs `80M keepawake`): about `1.25x` faster quant loop and `1.24x` faster total runtime, at cost of about `+50.4%` max working set.
4. Avoid `20000000` peak cap for full runs; it repeatedly caused late quant-loop cache-thrash.
5. For long benchmarks, keep host awake to avoid invalid timing data from standby/sleep gaps.

## 1.1 Execution Tracker (Updated 2026-02-07)

| Task | Status | Code Evidence | Validation Status | Next Gate |
|---|---|---|---|---|
| P0-0 Stage timing + read counters | Implemented | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:541`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:588`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1359` | Build passes (`dotnet build ... -c Release`), dataset-level equivalence not yet run | Run baseline/candidate benchmark and output diff |
| P0-1 Memoized reads in `countMatchedMasses` | Implemented | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:602`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:685` | Build passes, runtime delta pending benchmark | Measure `countMatchedMassesMs` and read count delta |
| P0-2 Memoized reads in isotopic compare path | Implemented | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:181`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:497`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:833` | Build passes, output parity pending benchmark | Compare isotopic stage timing and quant equivalence |
| P0-3 Binary closest-MS1 + remove duplicate MS metadata pass | Implemented + extended | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:180`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:260`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:937` | Build passes; `stream_sql` RT index mode benchmarked against legacy (major setup RAM/time improvement) | Keep `stream_sql` default with legacy rollback option |
| P0-4 Low-allocation 3D `initRTProfile` rewrite | Implemented | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:30` | Build passes; bounded probes show additional quant-loop gains, full `.quant` parity still pending | Complete full-run parity once memory gate is addressed |
| P0-5 Direct RT/intensity XIC path (remove Peak2D hot-loop allocations) | Implemented | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:65`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:380`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1139` | Build passes; additional quant-loop gain confirmed (e.g. `latest_60M` at 1200 peptides: `189133 ms` vs pre-change `222176 ms`) | Validate long-run parity + settle default peak-cap profile |
| P1-1 Fragment matching algorithm optimization | Implemented (code), validation pending | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:663`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:685` | Build passes, equivalence/perf run pending | Execute P1 benchmark and compare to P0 baseline |
| P1-2 Bounded peak cache | Implemented (code), full-run validation complete | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:845`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:860`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1679`, `src/PSMBasedQuantificationTIMs/CLIArgumentParsing.fs:45` | Build passes; full-run tradeoff validated (`no-cap` fastest, `80M` reduced memory) | Keep dual-profile guidance: performance=`unbounded`, balanced=`80M` |
| P1-3 Isotopic distribution memoization | Implemented (code), validation pending | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:647`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:833`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1369` | Build passes, cache hit/miss behavior and output parity pending benchmark | Run benchmark and inspect `PerfCaches` metrics + `.quant` diff |
| P1-4 Controlled parallelization | Not started | - | - | Add feature-flagged worker model |

Execution notes:

1. Current implementation progress is code-first; runtime benchmark and output-diff gates are blocked on selecting and running the fixed reference dataset contract in section 4.
2. Multiple successful compilation checks were run after each P0/P1 edit block: `dotnet build src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fsproj -c Release`.
3. TIMs runtime smoke test was attempted with legacy non-TIMS minimal data (`testsOld/ProteomIQon.Tests/data/PSMBasedQuantification`), but this dataset is not a valid benchmark gate.
4. Legacy params file is schema-incompatible for current parser (`Unexpected property 'XicExtraction'`).
5. Current default params run reaches quantification but exits on empty result sequence in CSV writer for this dataset.
6. Remaining performance and correctness gates require a real TIMs reference dataset that yields non-empty quant output.
7. TIMs benchmark dataset provided by user was executed on 2026-02-06 and compared to legacy benchmark log from 2025-02-17 using timestamp milestones:
8. Setup path improved strongly:
9. `Create RetentionTime index`: 1420.742s (old) -> 379.664s (current), ~3.74x faster.
10. `rt_start -> quant_start`: 1861.599s (old) -> 381.244s (current), ~4.88x faster.
11. Quantification loop is currently slower in this run window:
12. `quant_start -> 4400 peptides quantified`: 182.366s (old) vs 1478.821s (current), ~8.11x slower to same marker.
13. Current run did not reach `executing quantification:finished` before external timeout/abort; old run finished at 1532.772s after quant start.
14. Immediate follow-up hypothesis: default bounded peak cache policy (`bounded_lru`, max 1000) and/or fragment matching refactor introduced throughput regression in quant loop on this dataset.
15. Important comparability caveat: old log (`binned_spectra_1.000_log_old.txt`) runs with `PerformLabeledQuantification = Unlabeled`, while current benchmark params are `N15Labeling`; direct quant-loop speed comparison is therefore not strictly apples-to-apples.
16. Quant-loop A/B probes (all stopped after `500 peptides quantified`, same N15 params, same machine) were executed to isolate regressions:
17. `probe_default_newcode_20260206_153406`: `sec_to_100=109.088`, `sec_to_500=136.993`.
18. `probe_legacyFrag_20260206_150828`: `sec_to_100=104.337`, `sec_to_500=128.853`.
19. `probe_isoOff_20260206_151755`: `sec_to_100=95.062`, `sec_to_500=118.621`.
20. `probe_legacyFrag_isoOff_20260206_152558`: `sec_to_100=98.601`, `sec_to_500=121.983`.
21. `probe_default_afterIsoDefaultOff_20260206_154320`: `sec_to_100=108.085`, `sec_to_500=135.265`.
22. `probe_isoOff_rep2_20260206_155238`: `sec_to_100=114.427`, `sec_to_500=140.866`.
23. Interpretation: quant-loop benchmark noise is substantial (likely driven by random wavelet padding and run heterogeneity), but this probe matrix is sufficient to justify keeping new P1 optimizations behind toggles and prioritizing quant-loop-safe defaults.
24. Implemented quant-loop feature toggles for controlled rollout and isolation:
25. `PQ_TIMS_FRAGMENT_MATCH_MODE=legacy` enables legacy fragment matching.
26. `PQ_TIMS_ISO_CACHE_MODE=on` enables isotopic pattern memoization (default is disabled).
27. Existing `PQ_TIMS_PEAK_CACHE_MODE` and `PQ_TIMS_PEAK_CACHE_MAX` remain active for peak-read cache tuning.
28. Memory-safe benchmarking protocol now in use: short probes stopped at 100 or 500 peptides and cache cap kept at <=1000 to avoid RAM pressure.
29. Added in-loop stage progress logging every 100 peptides via `PerfProgress` to allow hotspot diagnosis without full-run completion.
30. Confirmed hotspot at 100 peptides (`probe_perfProgress100_20260206_160241`): `xicExtractionMs=95735` dominated quant loop (`quantLoopMs=104966`) while `countMatchedMassesMs` and isotopic compare were minor.
31. After low-allocation `initRTProfile` update and memory-safe defaults, matching 100-peptide probe (`probe_safe100_cache1000_20260207_101024`) showed improved quant-loop stage timing:
32. `quantLoopMs`: 104966 -> 93559 (~10.9% faster), `xicExtractionMs`: 95735 -> 85723 (~10.5% faster), with stable `readSpectrumPeaksDbReads=7742`.
33. Remaining priority is further XIC-path optimization under strict memory constraints; do not increase cache size beyond validated limits on this workstation.
34. Additional XIC/internal experiments were run and kept/reverted strictly by measured effect:
35. Kept: low-allocation `initGetProcessedXIC` rewrite (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:247`) and all-zero short-circuit in baseline subtraction (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:234`).
36. Reverted: one-pass no-dictionary 3D aggregation variant (regression in probe), `MzSearch`-first IM variant (regression in probe), bounded FIFO cache mode (regression in probe).
37. Added quant-loop progress cache counters for evidence-driven tuning:
38. `PerfProgress` now logs `peakCacheHits`, `peakCacheMisses`, `peakCacheEvictions` (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1307`).
39. Cache-churn evidence at 100 peptides (cache 1000): `peakCacheHits=40004`, `peakCacheMisses=7742`, `peakCacheEvictions=6742` (`probe_safe100_cacheStats_20260207_111723`).
40. Cache size sweep (memory-safe, stop at 100 peptides):
41. cache 600 -> `quantLoopMs=105360` (`probe_safe100_cache600_20260207_121320`) [worse]
42. cache 700 -> `quantLoopMs=84045` and `93205` replicates; default-mode recheck `85510` (`probe_safe100_default700_20260207_120545`) [best overall]
43. cache 1000 -> around `94313` to `95771` typical in current branch
44. cache 1200 -> `96914`; cache 1500 -> `98263` [worse]
45. Decision: default `PQ_TIMS_PEAK_CACHE_MAX` changed from 1000 to 700 for better quant-loop throughput and lower memory footprint (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:597`).
46. Longer-loop confirmation with new default:
47. `probe_safe500_default700_20260207_122157` reached 500 peptides at `quantLoopMs=117370` (~117.37s), with `xicExtractionMs=102861`.
48. Compared to earlier reference probe (`probe_default_newcode_20260206_153406`, 500 peptides in ~136.99s), this is an additional ~14.3% quant-loop improvement at the 500-peptide marker.
49. Investigated feasibility of m/z-bounded binary search in `initRTProfile`; sampled check on TIMs data (`tmp_check_sorted.fsx`) showed `tested=400 unsorted=400` for `ReadSpectrumPeaks(...).Peaks`, so m/z binary-search assumptions are invalid on this dataset.
50. Rejected experiment (regression): sorted-on-cache-miss + m/z-bounded scan (`probe_safe100_sortMzBound_20260207_135940`) increased `quantLoopMs` at 100 peptides to `145540` and `xicExtractionMs` to `138369`; reverted immediately.
51. Added algorithm-level parity check between old dictionary-based closest-m/z aggregation and new two-pass strategy (`tmp_compare_algo.fsx`): `mismatches=0` across 20,000 randomized spectra/windows.
52. Kept optimization in `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:30`: replace per-spectrum dictionary aggregation with two-pass closest-m/z + sum for same m/z key (allocation-reducing, semantics-equivalent under existing strict comparisons).
53. Kept optimization in `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:30`: reorder hot filter checks to m/z-first, then ion mobility, to skip `IonMobility.Value` for most peaks outside narrow m/z windows.
54. Kept optimization in `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:30`: hoist `mzLow/mzHigh/mzLock` and `imLow/imHigh/imLock` out of inner loops; avoid repeated property reads in hotspot loops.
55. Kept optimization in `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:52`: create `Peak2D` directly instead of `Peak1D` + `RtIndexEntry.AsPeak2D` in the tight profile loop.
56. Same-session benchmark progression (all memory-safe probes, timeout after progress markers, same params/cache defaults):
57. `probe_safe100_peak2dDirect_20260207_142715`: `quantLoopMs=97856`, `xicExtractionMs=90190` at 100 peptides.
58. `probe_safe100_twoPassMz_20260207_144007`: `quantLoopMs=97592`, `xicExtractionMs=88267` at 100 peptides.
59. `probe_safe100_mzFirst_20260207_145143`: `quantLoopMs=95121`, `xicExtractionMs=87363` at 100 peptides.
60. `probe_safe100_mzFirst_hoist_20260207_150254`: `quantLoopMs=93316`, `xicExtractionMs=85308` at 100 peptides; and `quantLoopMs=118735` at 500 peptides.
61. Relative improvement within this probe series: about `-4.6%` quant-loop time at 100 peptides (`97856 -> 93316`) and about `-5.4%` in `xicExtractionMs` (`90190 -> 85308`) after the kept XIC-loop changes.
62. Full-run candidate attempt with current code (`full_candidate_20260207_151619`) was executed to completion gate but did not finish before timeout; latest marker reached `4500 peptides quantified` with `quantLoopMs=2061791` and `xicExtractionMs=1943991`.
63. During this full-run attempt, process working set grew to about `19.6 GB` (`ProteomIQon.PSMBasedQuantificationTIMs.exe`) and was manually stopped to avoid workstation memory pressure; no `.quant` file was emitted.
64. Conclusion from this gate: quant-loop speed improvements are confirmed on bounded probes, but full-run correctness diff and end-to-end completion remain blocked by memory growth; next phases must include memory-containment work (or streaming/segmented output strategy) before full parity validation.
65. Added peak-count bound for LRU cache (feature-flagged): `PQ_TIMS_PEAK_CACHE_MAX_PEAKS` in `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs` (kept).
66. No-cap control after this change (`probe_noCap_control_20260207_215417`) showed no meaningful regression vs previous no-cap baseline at 1000 peptides (`quantLoopMs=234484`, `xicExtractionMs=210632`).
67. Memory-cap probe (`probe_memcap20M_20260207_211439`) showed major RAM reduction through 800 peptides (`workingSetMB` around `3-7 GB`) with comparable speed up to that point.
68. Longer 20M-cap probe (`probe_memcap20M_long_20260207_212556`) showed severe cache thrash from 900 peptides (`peakCacheResidentSpectra` collapsed to `118`, `peakCacheEvictions=27738`, `quantLoopMs=451695` at 900) and is not acceptable as default.
69. 60M-cap probe (`probe_memcap60M_20260207_214024`) reduced memory relative to no-cap (`~10.2 GB` around 1000 peptides) but slowed quant loop materially (`quantLoopMs=288607` at 1000) and is also not acceptable as default.
70. Tested `cache_ms1_only` policy (bypass caching for non-MS1 IDs) in `probe_ms1cacheOnly_noCap_20260207_220707`; effect was neutral-to-slightly-worse and this experiment was reverted.
71. Current recommendation: keep default as `PQ_TIMS_PEAK_CACHE_MAX=700` with unbounded peak count; use `PQ_TIMS_PEAK_CACHE_MAX_PEAKS` only as an emergency memory guard on constrained hosts, expecting potential throughput penalties.
72. Next memory phase should prioritize algorithmic containment (streamed result materialization and reduced in-memory trace retention) over stricter cache caps, to avoid quant-loop regressions.
73. Added explicit runtime CLI parameters (no env var requirement for normal use) and threaded them to `quantifyPeptides` as `RuntimeTuning`:
74. `--runtime-peptide-db-mode`, `--runtime-peak-cache-mode`, `--runtime-peak-cache-max`, `--runtime-peak-cache-max-peaks`, `--runtime-fragment-match-mode`, `--runtime-iso-cache-mode`, `--runtime-rt-index-mode` (`src/PSMBasedQuantificationTIMs/CLIArgumentParsing.fs`, `src/PSMBasedQuantificationTIMs/Program.fs`).
75. Environment variables were removed from runtime tuning in this module; runtime behavior is configured through CLI parameters only.
76. Implemented a non-cache memory optimization for setup: `stream_sql` RT index builder that streams `SpectrumID, Description` rows from mzLite and extracts only `MS:1000511` (ms level) + `MS:1000016` (scan start time), avoiding full `MassSpectrum` object materialization (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`).
77. `stream_sql` vs `legacy` setup comparison on same dataset:
78. `stream_sql`: setup `62043 ms` (`probe_rtidx_streamsql_20260207_230010`) and `62300 ms` (`probe_rtidx_streamsql_rep2_20260207_231956`), with `rows=586385`, `ms1Entries=6863`.
79. `legacy`: setup `361183 ms` (`probe_rtidx_legacy_20260207_230959`).
80. Setup-phase speedup is about `5.8x` and setup RAM at peptide 0 dropped from about `13.9 GB` (`legacy`) to about `0.5 GB` (`stream_sql`).
81. Quant-loop timing under `stream_sql` showed a modest slowdown in these probes (~6-8% at 100-500 peptides vs one legacy run), but total time-to-1000-peptides is still substantially better due setup savings.
82. Decision: keep `stream_sql` as default with automatic fallback to legacy on parse/empty-index failure; keep `--runtime-rt-index-mode legacy` for explicit rollback/testing.
83. CLI-parameter path validated end-to-end (`probe_cli_peakcap20M_20260207_233054`): runtime parameters correctly drove cache/DB/RT-index behavior without requiring environment variables.

## 2. Scope

Primary code scope:

1. `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`
2. `src/PSMBasedQuantification/PSMBasedQuantification.fs` (reference only)
3. `C:\Users\carol\source\repos\MzIO` (analysis and optional P2 follow-up only)

P0 and P1 do not require schema changes in MzIO.

## 3. Execution Constraints

1. Do not change scientific intent of XIC extraction, quantification, isotopic comparison, or filtering.
2. Keep changes isolated and reversible.
3. Make one commit per task block.
4. Run correctness and performance validation after each task block.
5. Keep diagnostics optional and disabled for benchmarks unless explicitly tested.

## 4. Reference Dataset and Benchmark Harness

## 4.1 Dataset contract

Use one fixed benchmark set for all before/after runs:

1. One representative TIMS `.mzlite` file
2. Matching scored PSM file (`.qpsm` / FragPipe-derived table)
3. Matching peptide DB
4. Fixed quantification params JSON

## 4.2 Benchmark command template

```powershell
$instrument = "<path-to-file.mzlite>"
$psm = "<path-to-file.qpsm>"
$db = "<path-to-peptide-db>"
$params = "<path-to-params.json>"
$out = "<path-to-output-dir>"

Measure-Command {
  dotnet run --project src/PSMBasedQuantificationTIMs -- `
    -i $instrument -ii $psm -d $db -p $params -o $out
}
```

## 4.3 Correctness comparison contract

Compare baseline `.quant` vs candidate `.quant` by key:

1. `StringSequence`
2. `GlobalMod`
3. `Charge`
4. `PepSequenceID`
5. `ModSequenceID`

Tolerance:

1. Relative: `1e-6`
2. Absolute: `1e-8`
3. Exact match for non-float columns

## 4.4 Mandatory pre-work (before P0)

1. Create branch: `perf/tims-p0-p1`
2. Capture baseline runtime (3 runs, median)
3. Capture baseline output file for diffing
4. Capture baseline counters (optional but recommended):
   - `dotnet-counters monitor --counters System.Runtime`

## 5. P0 Plan (Low Risk, High Impact)

## 5.1 P0-0: Add lightweight stage timing and read counters

Goal: Make bottleneck reduction measurable per task.

Files:

1. `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`

Steps:

1. Add `Stopwatch` timing around these stages in `quantifyPeptides`:
   - Setup/init
   - `countMatchedMasses`
   - XIC extraction (`getXIC` path)
   - Isotopic comparison
   - Final filtering/write
2. Wrap spectrum read function with a counter:
   - Increment for each actual `ReadSpectrumPeaks` DB call.
3. Emit summary metrics at end of run via logger.

Validation:

1. Build succeeds.
2. Runtime logs include stage durations and read count.
3. Output is unchanged.

Rollback:

1. Remove instrumentation only commit.

Exit criteria:

1. Baseline metrics recorded and repeatable.

## 5.2 P0-1: Move memoized spectrum reader earlier and use it in `countMatchedMasses`

Goal: Remove repeated uncached reads during fragment-support scoring.

Files:

1. `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`

Current references:

1. Memoization currently created at `readSpecPeaksWithMem` near `:705`
2. `countMatchedMasses` uses direct read at `:569`

Steps:

1. Move memoized reader initialization to immediately after reader/transaction initialization.
2. Replace `inReader.ReadSpectrumPeaks psm.PSMId` in `countMatchedMasses` with `readSpecPeaksWithMem psm.PSMId`.
3. Ensure there is exactly one memoized reader definition in function scope.

Validation:

1. Build succeeds.
2. Output `.quant` equivalence check passes.
3. Stage metrics show reduced spectrum-read count and reduced `countMatchedMasses` time.

Rollback:

1. Revert this commit if output mismatch is detected.

Exit criteria:

1. Correctness preserved.
2. Measurable reduction in `countMatchedMasses` runtime.

## 5.3 P0-2: Use memoized reader in isotopic comparison path (`getSpec`)

Goal: Avoid repeated uncached MS1 peak decode in isotopic KL comparison.

Files:

1. `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`

Current references:

1. `getSpec` currently accepts `reader` and calls `reader.ReadSpectrumPeaks` at `:163`
2. Isotopic compare entry point at `initComparePredictedAndMeasuredIsotopicCluster` around `:478`

Steps:

1. Change `getSpec` to accept a function parameter: `string -> Peak1DArray`.
2. Thread memoized reader into `initComparePredictedAndMeasuredIsotopicCluster`.
3. Update compare function closure creation in `quantifyPeptides` to pass memoized reader.
4. Keep functional behavior identical (same filtering and KL math).

Validation:

1. Build succeeds.
2. Output equivalence check passes.
3. Stage metrics show reduced isotopic-comparison read count/time.

Rollback:

1. Revert this isolated commit if any drift.

Exit criteria:

1. Same quant output within tolerance.
2. Fewer peak reads in isotopic branch.

## 5.4 P0-3: Replace linear closest-MS1 lookup with binary search and remove duplicate full MS metadata read

Goal: Cut repeated O(n) closest-RT searches and eliminate duplicate metadata traversal.

Files:

1. `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`

Current references:

1. Linear closest search at `getClosestMs1` (`:156`)
2. RT index build at `:580`
3. Separate full `ReadMassSpectra` pass at `:585-592`

Steps:

1. Replace `getClosestMs1` implementation with binary search over sorted RT array.
2. Switch closest-MS1 representation from `(float * MassSpectrum)[]` to sorted `(float * spectrumId)[]`.
3. Build this array directly from `retTimeIdxed` (already sorted by RT) to avoid second `ReadMassSpectra` pass.
4. Remove the explicit `ReadMassSpectra` + sort block used only for closest lookup.
5. Ensure isotopic comparison still resolves the same spectrum IDs for equivalent RT queries.

Validation:

1. Build succeeds.
2. Output equivalence check passes.
3. Stage metrics show reduced setup time and reduced isotopic comparison CPU.
4. Spot-check closest-ID parity on sampled RT values before/after.

Rollback:

1. Keep old linear helper in temporary local function during verification.
2. Revert if mismatched IDs impact quant output.

Exit criteria:

1. No extra `ReadMassSpectra` pass.
2. Binary closest lookup in place and validated.

## 5.5 P0-4: Rewrite 3D `Query.initRTProfile` to low-allocation imperative implementation

Goal: Remove dominant sequence-based full-scan overhead while preserving current behavior.

Files:

1. `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`

Current reference:

1. `Query.initRTProfile` at `:30-56`

Implementation requirements:

1. Preserve current semantics:
   - filter by IM and mz windows
   - aggregate by exact mz (sum intensities)
   - choose mz group closest to lock mz
   - fallback zero peak when empty
2. Avoid repeated sequence enumeration and `groupBy` allocations.
3. Use `ToArray()` on RT entries and preallocated output array.
4. Prefer direct loops plus compact dictionary for mz aggregation.
5. Keep output type and call sites unchanged.

Detailed steps:

1. Convert entries to array once (`RtIndexEntry.Search(...).ToArray()`).
2. Preallocate `profile` array.
3. For each entry:
   - get peaks from memoized reader
   - iterate once through peaks
   - apply IM and mz predicates
   - accumulate intensity by exact `Mz` in `Dictionary<float,float>`
4. Select closest mz key to lock value and write `Peak2D`.
5. Write zero lock peak when dictionary is empty.
6. Remove prior `Seq.tryHead`, `Seq.groupBy`, `Seq.map`, `Array.ofSeq` path.

Validation:

1. Build succeeds.
2. Output equivalence check passes on full dataset.
3. Runtime reduction in XIC extraction stage is clearly measurable.
4. No regressions in unlabeled and labeled modes.

Rollback:

1. Keep old implementation in feature branch until validation complete.
2. Revert commit if correctness drift appears.

Exit criteria:

1. P0 aggregate speedup >= 35% median on 3-run benchmark.
2. Numeric equivalence accepted.

## 5.6 P0 sign-off checklist

1. All P0 tasks committed separately.
2. Baseline vs P0 performance report captured.
3. Correctness diff report captured.
4. No unresolved regression.

## 6. P1 Plan (Medium Risk, Additional Throughput)

## 6.1 P1-0: Freeze P0 baseline for P1 comparisons

Steps:

1. Tag commit: `perf-p0-stable`.
2. Store P0 benchmark outputs and logs.
3. Use P0 as baseline for all P1 deltas.

## 6.2 P1-1: Optimize `countMatchedMasses` fragment matching algorithm under ppm tolerance

Goal: Reduce O(P x F) inner-loop overhead.

Files:

1. `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`

Current reference:

1. Matching predicate at `:572-575`

Steps:

1. Convert fragment mass container from linked list to sorted array once per peptide group.
2. Implement tolerance-aware lookup helper:
   - binary search insertion point around `peak.Mz`
   - scan local neighbors while delta within ppm tolerance
3. Replace `List.exists` per peak with helper.
4. Validate against current sums exactly/tolerantly.

Validation:

1. Unit-level parity test for helper on synthetic peaks/fragments.
2. Full output equivalence check.
3. Stage timing improvement in `countMatchedMasses`.

Rollback:

1. Revert helper usage if parity fails.

Exit criteria:

1. Correctness maintained.
2. `countMatchedMasses` improves measurably on large groups.

## 6.3 P1-2: Replace unbounded memoization with bounded cache policy

Goal: Control memory growth and GC pressure.

Files:

1. `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`

Steps:

1. Introduce bounded cache wrapper for peak reads (LRU or segmented FIFO with lock).
2. Expose capacity via CLI runtime parameter for tuning, e.g. `--runtime-peak-cache-max`.
3. Keep default conservative value (for example 500-1000 spectra) and allow override.
4. Add cache metrics logging:
   - hits
   - misses
   - evictions
5. Replace current `FSharpAux.Memoization.memoize` call with bounded cache wrapper.

Validation:

1. Output equivalence check.
2. Memory profile: lower peak heap/LOH vs P0 on long runs.
3. Throughput does not regress significantly.

Rollback:

1. Keep fallback path to unbounded memoization behind a runtime CLI mode switch.

Exit criteria:

1. Stable memory behavior on large datasets.
2. No correctness change.

## 6.4 P1-3: Memoize isotopic distribution generation

Goal: Reduce repeated MIDAs computations.

Files:

1. `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`

Current reference:

1. `generateIsotopicDistributionOfFormulaBySum` at `:470`
2. Called via isotopic compare init at `:478-481`

Steps:

1. Add cache keyed by stable peptide identity + charge.
2. Ensure key generation is deterministic and collision-safe for this scope.
3. Replace direct generation call with cached lookup.
4. Add hit/miss counters to logs.

Validation:

1. Output equivalence check.
2. Reduced CPU in isotopic comparison stage.

Rollback:

1. Disable cache path with feature toggle if needed.

Exit criteria:

1. No output drift.
2. Measurable compare-stage reduction.

## 6.5 P1-4: Controlled intra-file parallelization (feature-flagged)

Goal: Improve throughput on multi-core systems without breaking determinism.

Files:

1. `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`
2. Optionally `src/PSMBasedQuantificationTIMs/CLIArgumentParsing.fs` (if exposing CLI knob)

Design constraints:

1. Do not share a single `MzSQL` reader across worker threads.
2. Keep deterministic final result ordering.
3. Keep default mode sequential.

Steps:

1. Add feature toggle and degree-of-parallelism setting via CLI argument(s).
2. Implement worker context factory:
   - one reader per worker
   - one cache per worker (or read-only shared immutable structures where safe)
3. Parallelize grouped peptide processing with deterministic post-sort by original index.
4. Guard chart writing in parallel mode:
   - disable charts automatically in parallel mode, or
   - synchronize chart output paths and writes.
5. Add stress test mode for concurrency correctness.

Validation:

1. Sequential and parallel outputs equivalent within tolerance.
2. No SQLite/thread exceptions under repeated runs.
3. Throughput scales on 2/4/8 cores for representative datasets.

Rollback:

1. Toggle parallelism off by default.
2. Keep original sequential path intact.

Exit criteria:

1. Optional parallel mode is stable.
2. Performance gain validated where enabled.

## 6.6 P1 sign-off checklist

1. P1 tasks are independently toggleable.
2. P1 provides >= 15% additional gain over P0 median on benchmark set.
3. Correctness criteria still pass.

## 7. Execution Order and Timeline

Recommended order:

1. P0-0
2. P0-1
3. P0-2
4. P0-3
5. P0-4
6. P0 sign-off
7. P1-1
8. P1-2
9. P1-3
10. P1-4
11. P1 sign-off

Estimated timeline:

1. P0: 2-4 working days
2. P1: 4-8 working days

## 8. Deliverables

At completion, produce:

1. Updated source with P0/P1 changes.
2. Benchmark report (baseline vs P0 vs P1):
   - total runtime
   - stage timings
   - peak read count
   - cache metrics
3. Correctness comparison report:
   - mismatch counts
   - max absolute/relative deltas
4. Final recommendation on P2 readiness.

## 9. Go/No-Go Criteria

Go to release candidate only if:

1. Runtime targets met.
2. Correctness checks pass.
3. No new runtime instability.
4. Parallel mode remains off by default unless fully validated.


