# PSMBasedQuantificationTIMs Performance Optimization Report (Consolidated & Validated)

## 1. Executive Summary

This report consolidates:

1. Findings from `OPTIMIZATION_REPORT.md`
2. Independent source-level validation against:
   - `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs`
   - `src/PSMBasedQuantification/PSMBasedQuantification.fs`
   - `C:\Users\carol\source\repos\MzIO`

### Bottom line

- The dominant confirmed bottleneck is the 3D custom XIC extraction path in `Query.initRTProfile`, which performs full-spectrum filtering/grouping per RT entry (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:30`).
- Additional confirmed bottlenecks are repeated spectrum reads/decodes in `countMatchedMasses` and isotopic comparison paths (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:569`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:163`).
- There is duplicated run-metadata traversal (full `ReadMassSpectra` done twice through two separate paths) (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:580`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:585`).
- Several claims in `OPTIMIZATION_REPORT.md` are directionally right but overstate impact or risk scientific drift unless carefully constrained.

### Realistic target

- P0 changes (low risk): **35-55% total runtime reduction** on representative workloads.
- P1 changes (moderate): additional **15-30%**.
- P2 architecture changes (MzIO): higher upside, higher risk and effort.

---

## 2. Validation of `OPTIMIZATION_REPORT.md` Findings

### 2.1 Claim validation table

| Claim from `OPTIMIZATION_REPORT.md` | Validation | Evidence |
|---|---|---|
| 3D `initRTProfile` is the primary bottleneck due to `Seq.filter` + `Seq.groupBy` | **Confirmed** | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:30-56` |
| `countMatchedMasses` performs un-memoized peak reads | **Confirmed** | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:569`; memoization created later at `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:705` |
| Isotopic compare path performs un-memoized peak reads | **Confirmed** | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:162-163`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:485` |
| Main loop is sequential and misses parallelism | **Confirmed (opportunity), not always safe to parallelize directly** | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1148-1155`; shared `MzSQL` reader/connection constraints |
| Replacing custom 3D `initRTProfile` with MzIO `RtProfile` is a low-risk drop-in | **Partially validated; risk is understated** | MzIO IM path uses binary search on potentially non-IM-sorted subset: `C:\Users\carol\source\repos\MzIO\src\MzIO.Processing\MzIOLinq.fs:247-251`; comparator at `...MzIOLinq.fs:111-115`; binary search requires sorted order (`C:\Users\carol\source\repos\MzIO\src\MzIO.Processing\BinarySearch.fs:53`) |
| MzIO `ReadSpectrumPeaks` default `ionMobility=false` implies missing IM decode | **Rejected** | Peak metadata comes from DB JSON and decoder uses actual `Peak1DArray` metadata: `C:\Users\carol\source\repos\MzIO\src\MzIO.SQL\MzIOSQL.fs:289`; IM decode branch in decoder: `C:\Users\carol\source\repos\MzIO\src\MzIO\Binary\BinaryDataDecoder.fs:109-115` |
| Memoization likely causes memory pressure | **Plausible hypothesis** | Unbounded memoize call at `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:705`; payloads are full decoded spectra (`C:\Users\carol\source\repos\MzIO\src\MzIO.SQL\MzIOSQL.fs:282-293`) |
| HashSet replacement for fragment matching is a direct win | **Partially valid; naive HashSet is not equivalent with ppm tolerance** | Current tolerance check is distance-based: `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:574`; exact-key HashSet cannot represent range condition directly |
| 5-20x overall speedup estimate | **Not evidence-backed** | No profiling data provided; retained as upper-bound hypothesis only |

### 2.2 Meaningful additions retained from `OPTIMIZATION_REPORT.md`

1. Strong emphasis on `initRTProfile` regression as top priority.
2. Correct callout that `countMatchedMasses` should use cached peak reads.
3. Correct concern that isotopic comparison currently re-reads peaks.
4. Useful operational note: disable diagnostics (`diagCharts`) in production runs when not needed.

### 2.3 Adjustments/corrections applied

1. `RtProfile` drop-in replacement is not classified as low-risk without verifying IM ordering assumptions.
2. Parallelization proposals are downgraded to medium/high risk due shared reader, cache behavior, and output determinism.
3. Impact estimates are normalized to conservative, measurable targets.

---

## 3. Current Architecture and Data Flow (Validated)

## 3.1 End-to-end flow

1. CLI invokes quantification (`src/PSMBasedQuantificationTIMs/Program.fs:16`, `src/PSMBasedQuantificationTIMs/Program.fs:60`).
2. `quantifyPeptides` initializes outputs and lookup DB (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:523`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:545`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:550`).
3. Opens `MzSQL` reader and transaction (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:555-558`).
4. Builds RT index via MzIO (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:580` -> `C:\Users\carol\source\repos\MzIO\src\MzIO.Processing\Query.fs:15` -> `C:\Users\carol\source\repos\MzIO\src\MzIO.Processing\MzIOLinq.fs:179`).
5. Reads mass spectra again and sorts MS1 by scan time (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:585-592`).
6. Reads scored PSMs and computes mz correction (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:597-649`).
7. Initializes XIC extraction and peak detection (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:705-710`).
8. Groups by peptide-ion and quantifies each group (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1126-1155`).
9. Quality filtering and output write (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1157-1193`).

## 3.2 Critical inner-loop data flow

1. `countMatchedMasses` reads MS2 peaks and scores fragment support (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:559-578`).
2. `average` computes weighted RT/IM and calls XIC extractor (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:283-301`).
3. `initGetProcessedXIC` invokes `Query.initRTProfile` for RT/mz/IM window (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:227-255`).
4. Peak detection + fit (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:726`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:733`).
5. Isotopic compare selects closest MS1 and reads spectrum peaks again (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:478-515`).

---

## 4. Critical Path and Hotspots

## 4.1 Confirmed hotspots

1. **3D RT profile extraction regression**
   - Location: `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:30-56`
   - Behavior: full-spectrum sequence filter + grouping + sum per RT entry.
   - Why slow: O(total_peaks_per_spectrum) per entry plus iterator/grouping allocations.

2. **Uncached peak reads in `countMatchedMasses`**
   - Location: `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:569`
   - Behavior: direct `inReader.ReadSpectrumPeaks` for each PSM.
   - Why slow: each call triggers DB query + JSON metadata parse + binary decode (`C:\Users\carol\source\repos\MzIO\src\MzIO.SQL\MzIOSQL.fs:282-293`).

3. **Linear closest-MS1 lookup**
   - Location: `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:156`
   - Behavior: `Array.minBy` over all MS1 entries for each isotopic comparison.
   - Why slow: O(#MS1) repeated many times.

4. **Duplicate full mass-spectrum metadata traversal**
   - Location: `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:580` and `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:585`
   - Behavior: first for RT index build (internally reads mass spectra), then explicit second full read/sort.
   - Why slow: repeated SQL + JSON deserialization.

5. **Isotopic comparison re-reads peaks bypassing memoized function**
   - Location: `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:163`
   - Behavior: direct `reader.ReadSpectrumPeaks`.
   - Why slow: repeat decode on often-overlapping nearby MS1 IDs.

## 4.2 Hypotheses (profile required)

1. Unbounded memoization (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:705`) may produce LOH pressure for large TIMs spectra.
2. Fragment matching complexity (`List.exists` with tolerance) may become visible on large peak/fragment cardinalities (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:572-575`).
3. Isotopic distribution generation may be worth caching (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:470`).

---

## 5. 3D vs 2D Comparison Findings

## 5.1 Confirmed regressions/misses

1. `initRTProfile` implementation changed from loop + m/z binary search in 2D to sequence full scan/grouping in 3D.
   - 2D: `src/PSMBasedQuantification/PSMBasedQuantification.fs:30-41`
   - 3D: `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:30-56`

2. Inherited bottlenecks in both variants:
   - uncached `countMatchedMasses` reads
   - linear closest-MS1 selection
   - duplicate mass-spectrum metadata traversal

3. 3D hard-coded IM extraction window:
   - `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:708` passes `0.05`.
   - This is a tuning opportunity and should likely be parameterized.

## 5.2 Important caveat on "use MzIO RtProfile directly"

`OPTIMIZATION_REPORT.md` suggested direct use of MzIO `RtProfile`.  
Validated caveat:

- MzIO `RtProfile` IM branch applies `IonMobilitySearch` (binary search) after `MzSearch` on an intermediate array (`C:\Users\carol\source\repos\MzIO\src\MzIO.Processing\MzIOLinq.fs:247-251`).
- `BinarySearch.Search` assumes data sorted by the comparison key (`C:\Users\carol\source\repos\MzIO\src\MzIO.Processing\BinarySearch.fs:53`).
- The intermediate array is known sorted by m/z, not proven sorted by ion mobility.

Therefore direct replacement is **not low-risk without explicit ordering verification**.

---

## 6. MzIO Root Function Findings

## 6.1 Confirmed read/decode cost centers

1. Spectrum metadata read path:
   - `ReadMassSpectra` -> SQL `SELECT Description` -> JSON deserialize each row (`C:\Users\carol\source\repos\MzIO\src\MzIO.SQL\MzIOSQL.fs:266-279`, `C:\Users\carol\source\repos\MzIO\src\MzIO.SQL\MzIOSQL.fs:380-382`).

2. Peak read path:
   - `ReadSpectrumPeaks` -> SQL `SELECT PeakArray, PeakData` -> JSON parse `Peak1DArray` metadata -> binary decode (`C:\Users\carol\source\repos\MzIO\src\MzIO.SQL\MzIOSQL.fs:282-293`, `C:\Users\carol\source\repos\MzIO\src\MzIO.SQL\MzIOSQL.fs:390-392`).

3. Decoder allocates peak objects per element:
   - `Array.map2/Array.map3` to construct `Peak1D` instances (`C:\Users\carol\source\repos\MzIO\src\MzIO\Binary\BinaryDataDecoder.fs:89`, `C:\Users\carol\source\repos\MzIO\src\MzIO\Binary\BinaryDataDecoder.fs:117`).

## 6.2 Validated correction from previous report

`ionMobility=false` default in `SelectPeak1DArray` does **not** imply IM data loss in standard read path, because metadata from DB JSON drives decoder behavior (`C:\Users\carol\source\repos\MzIO\src\MzIO.SQL\MzIOSQL.fs:289`).

---

## 7. Consolidated Optimization Proposals

| Priority | Finding Type | Exact Location | Current Behavior | Proposed Change | Expected Impact | Risk to Correctness | Validation Approach |
|---|---|---|---|---|---|---|---|
| P0 | Confirmed | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:30` | Full-spectrum seq filter/groupBy in 3D `initRTProfile` | Rewrite as imperative, low-allocation loop: m/z narrowing first, IM filtering second, preserve current aggregation semantics | High; largest single gain (approx 25-45% overall) | Low-Medium | Golden-output diff on `.quant`; XIC trace equality checks |
| P0 | Confirmed | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:42` | `Seq.tryHead` then second enumeration | Single-pass fold/materialization | Medium | Low | Per-peptide XIC equivalence tests |
| P0 | Confirmed | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:569` and `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:705` | `countMatchedMasses` uses uncached reads | Move memoized peak reader init before `countMatchedMasses`; use cache there | Medium-High in PSM-heavy workloads | Low | Instrument `ReadSpectrumPeaks` call count before/after |
| P0 | Confirmed | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:163` | Isotopic compare uses direct uncached read | Thread memoized read function through `getSpec`/compare path | Medium | Low | Verify identical KL metrics and isotopic arrays |
| P0 | Confirmed | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:156` | Linear closest-MS1 lookup | Replace with binary-search lookup on sorted RT list | Medium (5-15%) | Low | Compare selected MS1 IDs for sampled RTs |
| P0 | Confirmed | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:580`, `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:585` | Reads mass spectra metadata twice | Reuse RT index data to avoid second full metadata pass | Medium | Low | Startup-phase timings + functional equivalence |
| P1 | Hypothesis | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:559` | Fragment match uses nested `List.exists` with tolerance | Replace with sorted array + range-based lookup (not exact HashSet) | Low-Medium | Low-Medium | Compare matched sums and downstream quant differences |
| P1 | Hypothesis | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:705` | Unbounded memoization | Bounded/LRU cache for peak arrays | Throughput stability + lower memory | Low-Medium | GC counters + memory high-watermark |
| P1 | Hypothesis | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:470` | Recompute isotopic distribution repeatedly | Cache by `(formula/sequence, charge)` | Low-Medium | Low | KL and isotopic pattern equality |
| P1 | Hypothesis | `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:1148` | Per-peptide loop sequential | Controlled parallelism with per-worker reader/caches and deterministic merge | Medium-High | Medium-High | Determinism tests + stress runs + perf scale tests |
| P2 | Confirmed | `C:\Users\carol\source\repos\MzIO\src\MzIO.SQL\MzIOSQL.fs:58` | Schema lacks explicit run-level helper columns/indexes for RT build path | Add indexed persisted fields (e.g., scan time, ms level) to avoid heavy JSON parse for index creation | High for all RT-index users | Medium | MzIO regression suite + migration tests |
| P2 | Confirmed | `C:\Users\carol\source\repos\MzIO\src\MzIO.SQL\MzIOSQL.fs:282`, `C:\Users\carol\source\repos\MzIO\src\MzIO\Binary\BinaryDataDecoder.fs:225` | Always decode full spectrum to objects | Add filtered decode APIs returning only required range or flat arrays | High in XIC-heavy workloads | Medium-High | API tests + perf microbench + integration tests |

---

## 8. Phased Implementation Plan (P0/P1/P2)

## 8.1 P0: Fast, low-risk speedups

Tasks:

1. Optimize 3D `initRTProfile` implementation (preserve current output behavior).
2. Move memoized peak reader setup earlier and use in `countMatchedMasses`.
3. Use memoized reader in isotopic comparison path.
4. Replace linear closest-MS1 selection with binary search.
5. Remove duplicate full mass-spectrum metadata read.

Success criteria:

1. End-to-end runtime improvement >= 35% on reference dataset.
2. `ReadSpectrumPeaks` call count materially reduced.
3. `.quant` output numerically equivalent within tolerance.

Rollback/safety:

1. Keep legacy XIC path behind feature flag.
2. Keep pre-change closest-MS1 implementation callable for A/B verification.

## 8.2 P1: Medium refactors

Tasks:

1. Optimize fragment matching data structure under ppm tolerance.
2. Introduce bounded cache policy for peak memoization.
3. Add isotopic distribution memoization.
4. Evaluate safe parallelization strategy with per-worker readers.

Success criteria:

1. Additional 15-30% improvement over P0.
2. Stable memory footprint, no OOM on large runs.
3. Deterministic outputs preserved.

Rollback/safety:

1. Guard each change with independent config toggle.
2. Keep single-thread fallback as default until validated.

## 8.3 P2: MzIO architecture work

Tasks:

1. Introduce schema/API path that avoids expensive per-row JSON decode for RT indexing.
2. Add partial/filter decode APIs to reduce object creation.

Success criteria:

1. Significant additional gains in I/O-heavy runs (>20%).
2. Lower allocation rate and GC pause time.

Rollback/safety:

1. Versioned MzIO changes with backward-compatible fallback path.
2. Migration scripts and verification tests before rollout.

---

## 9. Validation and Benchmark Plan

## 9.1 Correctness

1. Generate baseline `.quant` output on fixed dataset.
2. Compare keyed rows by `(Sequence, GlobalMod, Charge, PepSequenceID, ModSequenceID)`.
3. Tolerances:
   - Relative: `1e-6`
   - Absolute: `1e-8`
4. Exact equality for categorical/non-float fields.
5. Include edge cases: no peaks, failed fits, NaN-heavy groups.

## 9.2 Performance

1. Measure end-to-end runtime and phase timings:
   - `countMatchedMasses`
   - XIC extraction
   - isotopic comparison
   - write/filter
2. Track:
   - `ReadSpectrumPeaks` call count
   - allocations and GC (gen2, LOH trends)
3. Tools:
   - `dotnet-trace`
   - `dotnet-counters`
   - optional BenchmarkDotNet for focused microbenchmarks

## 9.3 Bottleneck ranking confidence (post-validation expectation)

1. High: 3D `initRTProfile`
2. High: repeated peak decode/read paths
3. Medium: metadata double-read + linear closest-MS1
4. Medium: cache/GC pressure
5. Medium-Low: fragment matching micro-optimizations

---

## 10. Risks, Assumptions, Open Questions

## 10.1 Risks

1. Changing 3D XIC logic can alter intensity aggregation if semantics drift.
2. Parallelization can violate determinism and stress SQLite concurrency.
3. Cache bounding may increase I/O misses if undersized.

## 10.2 Assumptions

1. Production runs usually use `.mzlite` inputs in this module (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:555`).
2. TIMs datasets are substantially larger per spectrum than 2D datasets (domain assumption, not code-confirmed).

## 10.3 Open questions

1. Is current `groupBy mz -> sum intensity` in 3D `initRTProfile` scientifically required, or can intensity selection use alternate aggregation?
2. What is typical peptide-ion count per run (for parallel ROI)?
3. Is `diagCharts` enabled in production workloads (can be large overhead)?
4. Should IM window (`0.05`) be parameterized for instrument/method variability (`src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs:708`)?
5. If using MzIO IM binary-search path, can IM ordering guarantees be formally established?

