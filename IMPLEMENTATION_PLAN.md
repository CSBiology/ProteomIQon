# PSMBasedQuantificationTIMs Performance Implementation Plan (P0/P1)

## 1. Objective

Deliver measurable runtime improvements in `src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs` without changing scientific correctness.

1. P0 target: at least 35% end-to-end runtime reduction on reference dataset.
2. P1 target: additional 15-30% improvement over P0.
3. Output requirement: numerical equivalence within defined tolerance.

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
2. Expose capacity via environment variable for tuning, e.g. `PQ_TIMS_PEAK_CACHE_MAX`.
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

1. Keep fallback path to unbounded memoization behind env switch.

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

1. Add feature toggle and degree-of-parallelism setting (env var or CLI arg).
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
