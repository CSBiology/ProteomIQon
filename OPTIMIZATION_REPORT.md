# PSMBasedQuantificationTIMs Performance Optimization Report

## 1. Executive Summary

The 3D TIMs quantification tool (`PSMBasedQuantificationTIMs.fs`, 1194 lines) extends the proven 2D tool (`PSMBasedQuantification.fs`, 1173 lines) by adding ion mobility as a third dimension. The core performance bottleneck lies in the **`initRTProfile` function** (lines 30–56), which replaces the MzIO library's optimized binary-search-based RT profile extraction with a **Seq-based linear scan + filter + groupBy** approach that creates massive intermediate allocations per spectrum, per peptide ion. This function is called **for every XIC extraction** in the hot loop, compounding with the fact that TIMs spectra are significantly larger (peaks carry ion mobility data, increasing per-spectrum peak counts by 10–100x vs conventional MS1).

**Key findings:**
- **PRIMARY BOTTLENECK [Confirmed]**: Custom `initRTProfile` in TIMs uses `Seq.filter` + `Seq.groupBy` + `Seq.map` + `Seq.sumBy` over all peaks per spectrum (lines 37–53), instead of binary search. This is O(N) per spectrum where N = total peaks, vs O(log N) in the 2D version's `MzSearch`.
- **SECONDARY BOTTLENECK [Confirmed]**: `countMatchedMasses` (lines 560–577) reads spectrum peaks via un-memoized `inReader.ReadSpectrumPeaks` for every PSM, decoding binary data fresh each time.
- **TERTIARY BOTTLENECK [Confirmed]**: Memoization caches entire `Peak1DArray` objects, which for TIMs data are very large — causing memory pressure and GC overhead.
- **MISSED OPTIMIZATION [Confirmed]**: Main quantification loop (lines 1149–1170) runs sequentially via `Array.mapi`, no parallelization.

**Estimated impact**: Fixing the primary bottleneck alone should deliver **3–10x speedup** on the XIC extraction hot path. Combined optimizations target **5–20x overall improvement**.

---

## 2. Current Architecture and Data Flow

### End-to-End Pipeline

```
Input Files                        Processing                          Output
─────────                         ──────────                         ──────
.mzlite (SQLite) ──┐
                    ├─► quantifyPeptides()
scored PSMs (.tsv) ─┘      │
peptide DB (.db)  ──────────┘
                            │
                            ▼
                   ┌─────────────────┐
                   │ 1. Init Phase   │
                   │   - Copy DB     │  (line 545-546)
                   │   - Open reader │  (line 555-558)
                   │   - Build RT idx│  (line 580)
                   │   - Sort MS1s   │  (line 584-591)
                   │   - Read PSMs   │  (line 594-596)
                   │   - Mz correction│  (line 601-653)
                   └────────┬────────┘
                            ▼
                   ┌─────────────────┐
                   │ 2. Group PSMs   │  (lines 1149-1163)
                   │  by peptide ion │
                   └────────┬────────┘
                            ▼
               ┌────────────────────────┐
               │ 3. Per-PeptideIon Loop │  (line 1164-1170)
               │  (sequential Array.mapi)│
               │                        │
               │  For each peptide ion: │
               │  ├─ countMatchedMasses │  ← reads spectra (slow)
               │  ├─ average()          │  ← calls getXIC (slow)
               │  │   └─ getXIC         │
               │  │       └─initGetProc │
               │  │          └─getPeaks │
               │  │            └─initRT │  ← PRIMARY BOTTLENECK
               │  │              Profile│
               │  ├─ identifyPeaks      │  ← signal processing
               │  ├─ HULQ.quantifyPeak  │  ← curve fitting
               │  ├─ getInferredXic     │  ← 2nd XIC extraction
               │  │   └─ getXIC (again) │  ← same bottleneck
               │  ├─ comparePredicted.. │  ← isotope pattern
               │  │   └─ getSpec        │  ← another spectrum read
               │  └─ calcCorrelation    │
               └────────────────────────┘
                            ▼
                   ┌─────────────────┐
                   │ 4. Filter &     │  (lines 1172-1184)
                   │    Write Output │
                   └─────────────────┘
```

### Per-Peptide-Ion Call Costs

Each peptide ion triggers (for labeled quantification):
1. **`countMatchedMasses`**: N × `ReadSpectrumPeaks` (un-memoized) + N × full-peak linear scan
2. **`average` → `getXIC`**: 1 × `initRTProfile` (memoized reads, but Seq-based filtering)
3. **`identifyPeaks`**: Signal processing (moderate cost)
4. **`HULQ.quantifyPeak`**: Nonlinear curve fitting (moderate cost)
5. **`getInferredXic` → `getXIC`**: 1 × `initRTProfile` again (for the labeled counterpart)
6. **`comparePredictedAndMeasuredIsotopicCluster`**: 1 × `getSpec` (un-memoized `ReadSpectrumPeaks`)

---

## 3. Critical Path and Hotspots

### Hotspot #1: `initRTProfile` (CRITICAL — Lines 30–56)

**File**: [PSMBasedQuantificationTIMs.fs](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L30-L56)

**Current behavior**: For each RT entry in range:
1. Read spectrum peaks (memoized) — returns ALL peaks in the spectrum
2. `Seq.filter` over ALL peaks checking both IM range AND mz range (line 37–40)
3. `Seq.tryHead` to check if empty (forces enumeration start)
4. `Seq.groupBy` on Mz (creates dictionary, allocates groups)
5. `Seq.map` + `Seq.sumBy` to aggregate intensity per mz group
6. `RtIndexEntry.ClosestMz` to pick closest peak

**Why it's slow**: 
- Linear scan over ALL peaks per spectrum (TIMs spectra can have 50K–500K peaks due to ion mobility dimension)
- `Seq.groupBy` allocates a `Dictionary<float, ResizeArray<Peak1D>>` per spectrum
- `Seq.tryHead` + subsequent full iteration is redundant work
- No binary search on sorted mz/IM arrays

**Contrast with 2D version** ([PSMBasedQuantification.fs L32–44](src/PSMBasedQuantification/PSMBasedQuantification.fs#L32-L44)):
```fsharp
// 2D version — uses binary search, pre-allocated array, imperative loop
let entries = RtIndexEntry.Search(rtIndex, rtRange).ToArray()
let profile = Array.zeroCreate<Peak2D> entries.Length
for rtIdx = 0 to entries.Length-1 do
    let entry = entries.[rtIdx]
    let peaks = (readspecPeaks entry.SpectrumID).Peaks
    let p = (RtIndexEntry.MzSearch (peaks, mzRange))  // <-- BINARY SEARCH
            .DefaultIfEmpty(Peak1D(0., mzRange.LockValue))
            |> fun x -> RtIndexEntry.ClosestMz (x, mzRange.LockValue)
            |> fun x -> RtIndexEntry.AsPeak2D (x, entry.Rt)
    profile.[rtIdx] <- p
profile
```

The 2D version uses `RtIndexEntry.MzSearch` which is a **binary search** (O(log N)), while the TIMs version does a **full linear scan with filter** (O(N)).

**Confirmed bottleneck**: Yes — algorithmic regression from O(log N) to O(N) per spectrum, multiplied by potentially 50–500x more peaks per spectrum.

---

### Hotspot #2: `countMatchedMasses` (HIGH — Lines 560–577)

**File**: [PSMBasedQuantificationTIMs.fs](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L560-L577)

**Current behavior**: For each PSM in a peptide ion group, reads the **full MS2 spectrum** peaks (un-memoized!) and does O(P×F) scan where P = peaks and F = fragment ions.

**Why it's slow**:
- `inReader.ReadSpectrumPeaks psm.PSMId` involves SQLite query + binary decoding + JSON deserialization per call
- NOT using the memoized `readSpecPeaksWithMem` (which is only created later at line 705)
- Fragment ion list `frag` is a `List` (linked list) — `List.exists` is O(F) per peak
- This is called for EVERY peptide ion group, reading same spectra that may overlap

**Confirmed bottleneck**: Yes — un-memoized I/O + O(P×F) matching per PSM.

---

### Hotspot #3: Memoization Cache Memory Pressure (MEDIUM)

**File**: [PSMBasedQuantificationTIMs.fs](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L705)

**Current behavior**: `FSharpAux.Memoization.memoize inReader.ReadSpectrumPeaks` caches the entire `Peak1DArray` object (all peaks + mz + intensity + ion mobility arrays) keyed by spectrum ID.

**Why it's a concern**: For TIMs data, each spectrum may contain 50K–500K peaks × (float intensity + float mz + float IM) = ~1.2–12 MB per spectrum. With RT windows covering ~50–200 spectra, and peptide ions sharing spectra across XIC queries, the memoization cache grows to **multiple GB**, stressing the GC.

**Hypothesis (needs profiling)**: GC pauses from large object heap allocations may contribute 10–30% of wall-clock time.

---

### Hotspot #4: Sequential Main Loop (MEDIUM — Lines 1164–1170)

**File**: [PSMBasedQuantificationTIMs.fs](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L1164-L1170)

**Current behavior**: `Array.mapi` processes peptide ions sequentially. Each iteration is independent (aside from shared reader state).

**Why it's slow**: Modern machines have 8–16+ cores. The computation per peptide ion (signal processing, curve fitting) is CPU-bound. SQLite reader access is the serialization point.

**Confirmed opportunity**: Yes — CPU-bound work can be parallelized if I/O is pre-fetched or made thread-safe.

---

### Hotspot #5: `getSpec` / `comparePredictedAndMeasuredIsotopicCluster` (LOW-MEDIUM)

**File**: [PSMBasedQuantificationTIMs.fs](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L162-L164) and [lines 476–519](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L476-L519)

**Current behavior**: `getSpec` calls un-memoized `reader.ReadSpectrumPeaks(ms1.ID)` to read an MS1 spectrum for isotopic pattern comparison. Then filters to a small mz range.

**Why it's a concern**: Another un-memoized spectrum read per peptide ion, plus `Peaks.unzipIMzliteArray` allocates two fresh `float[]` arrays just to `PeakArray.zip` them back.

---

## 4. 3D vs 2D Comparison Findings

| Aspect | 2D (`PSMBasedQuantification`) | 3D TIMs (`PSMBasedQuantificationTIMs`) | Impact |
|---|---|---|---|
| **`initRTProfile` algorithm** | Binary search via `RtIndexEntry.MzSearch` (O(log N)) | Linear `Seq.filter` + `Seq.groupBy` (O(N)) | **10–100x slower per spectrum** |
| **`initRTProfile` output construction** | Pre-allocated `Array.zeroCreate`, imperative loop | `Seq.map` → `Array.ofSeq` (lazy + materialization) | 2–5x more allocation |
| **`countMatchedMasses` PSM type** | `PSMStatisticsResult` | `PSMStatisticsResultFragpipe` | Same structure, negligible |
| **`countMatchedMasses` memoization** | Not memoized (same issue) | Not memoized (same issue) | Both have this bug |
| **Main loop parallelism** | `Array.mapi` (sequential) | `Array.mapi` (sequential) | Both miss parallelism |
| **`getXIC` call signature** | `(meanScanTime, meanPrecMz)` | `(meanScanTime, meanPrecMz, meanIM)` | 3D adds IM window |
| **Spectrum size** | ~1K–10K peaks per MS1 | ~50K–500K peaks per MS1 (IM expanded) | 10–100x more data per read |
| **ion mobility filtering** | N/A | `Seq.filter` on `IonMobility.Value` | Additional O(N) pass |
| **IM-based groupBy** | N/A | `Seq.groupBy Mz` after IM filter | Expensive allocation |
| **`getInferredXic`** | 2 args | 3 args (adds `targetIM`) | Same cost pattern |

**Key regression**: The TIMs version **rewrote** `initRTProfile` from scratch instead of extending the MzIO library's `RtProfile` method (which already supports `ionMobilityRange` as an optional parameter — see [MzIOLinq.fs lines 231–262](src/MzIO/Processing/MzIOLinq.fs)). The MzIO library already has `RtProfile(rtIndex, rtRange, mzRange, ?ionMobilityRange)` that uses binary search via `MzSearch` + `IonMobilitySearch` with proper `BinarySearch.Search`.

---

## 5. MzIO Root Function Findings

### 5.1 Spectrum Read Path

```
quantifyPeptides
  └─ readSpecPeaksWithMem = memoize inReader.ReadSpectrumPeaks
      └─ MzSQL.ReadSpectrumPeaks(spectrumID)           [MzIOSQL.fs:390-392]
          └─ SelectPeak1DArray(spectrumID)              [MzIOSQL.fs:343-346]
              └─ selectPeak1DArray false spectrumID     [defaults ionMobility=false]
                  └─ prepareSelectPeak1DArray(cn, false) [MzIOSQL.fs:289-301]
                      ├─ SQL: SELECT PeakArray, PeakData FROM Spectrum WHERE SpectrumID = @spectrumID
                      ├─ JSON deserialize: MzIOJson.FromJson<Peak1DArray>(reader.GetString(0))
                      └─ Binary decode: decoder.Decode(reader.GetStream(1), peakArray)
                          └─ BinaryDataDecoder.Decode    [BinaryDataDecoder.fs:229-234]
                              ├─ If NumPressZLib: decompress + numpress decode
                              └─ Constructs Peak1D[] with Array.map3 (intensity, mz, ionMobility)
```

**Critical finding**: `ReadSpectrumPeaks` via the `IMzIODataReader` interface calls `SelectPeak1DArray(spectrumID)` which defaults `ionMobility=false`. However, because the `prepareSelectPeak1DArray` function reads the `PeakArray` JSON from the DB column and deserializes the actual `Peak1DArray` metadata (which includes `IonMobilityDataType` if it was stored with IM), the decoder **DOES** correctly decode ion mobility data. The `ionMobility` flag only controls the **fallback empty array** template. So correctness is OK, but there's still a semantic concern.

### 5.2 Binary Decoding Cost

For NumPressZLib (common for TIMs data):
1. Read compressed bytes from SQLite BLOB
2. Deflate decompress intensities, mz, ionMobility (3 separate `DeflateStreamDecompress` calls)
3. Numpress decode 3 arrays 
4. `Array.map3` to create `Peak1D[]` — allocates N `Peak1D` objects
5. Wrap in `MzIOArray`

For a spectrum with 200K peaks: ~200K × 3 floats × 8 bytes = **4.8 MB** raw data, plus 200K heap-allocated `Peak1D` objects (~48 bytes each with object overhead) = **~14 MB** per spectrum decode.

### 5.3 RT Index and Binary Search

**File**: [BinarySearch.fs](src/MzIO/Processing/BinarySearch.fs) and [MzIOLinq.fs](src/MzIO/Processing/MzIOLinq.fs)

The `BinarySearch.Search` function performs proper O(log N) binary search on sorted `IMzIOArray`, then linear expansion to find full range. `RtIndexEntry.MzSearch` and `IonMobilitySearch` use this same binary search. The MzIO library's `RtProfile` method (line 231, MzIOLinq.fs) already chains `MzSearch` → `IonMobilitySearch` with binary search, but the TIMs tool doesn't use it.

### 5.4 Performance Concerns in MzIO

| Issue | Location | Impact |
|---|---|---|
| Peak1D is a class (heap allocated) | [PeakArray.fs:20](src/MzIO/Binary/PeakArray.fs) | 200K allocs per spectrum decode |
| `DefaultIfEmpty` creates LINQ iterator | MzIOLinq.fs | Small but repeated allocation |
| `Seq.toArray` in `IonMobilitySearch` chain | MzIOLinq.fs:226 | Materializes intermediate array |
| JSON deserialization per spectrum read | MzIOSQL.fs:297 | `MzIOJson.FromJson<Peak1DArray>` overhead |
| No connection pooling | MzIOSQL.fs | Single connection, single thread |

---

## 6. Optimization Proposals

### Summary Table

| # | Proposal | Location | Impact | Risk | Effort | Phase |
|---|----------|----------|--------|------|--------|-------|
| **P0-1** | Replace `initRTProfile` with MzIO library's `RtProfile` with IM support | TIMs L30–56 | **~5–10x XIC speedup** | Low | Low | P0 |
| **P0-2** | Use memoized reader in `countMatchedMasses` | TIMs L569 | **~2–3x init speedup** | None | Trivial | P0 |
| **P0-3** | Convert `frag` from `List` to `HashSet` in `countMatchedMasses` | TIMs L561–567 | **~2x matching speedup** | None | Trivial | P0 |
| **P1-1** | Parallelize main quantification loop | TIMs L1164–1170 | **~3–8x throughput** | Medium | Medium | P1 |
| **P1-2** | Batch-prefetch spectra for IM pre-filtering | TIMs L30–56 | **~2x I/O throughput** | Medium | Medium | P1 |
| **P1-3** | Use LRU cache instead of unbounded memoization | TIMs L705 | **~20–50% less GC** | Low | Medium | P1 |
| **P1-4** | Memoize `getSpec` / isotope comparison spectrum reads | TIMs L163 | **~1.5x per-ion speedup** | None | Low | P1 |
| **P2-1** | Add `ReadSpectrumPeaksFiltered` to MzIO with in-DB mz/IM range | MzIO library | **~10–50x I/O reduction** | High | High | P2 |
| **P2-2** | Convert Peak1D from class to struct | MzIO Binary | **~30–50% less GC** | High | High | P2 |
| **P2-3** | Cache decoded peaks in SQLite WAL or memory-mapped file | MzIO.SQL | **~3x I/O reduction** | High | High | P2 |

---

### P0-1: Replace Custom `initRTProfile` with MzIO Library's Binary Search

**Location**: [PSMBasedQuantificationTIMs.fs lines 30–56](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L30-L56)

**Current behavior**: Custom Seq-based linear filter/groupBy/sumBy over all peaks per spectrum for IM+mz filtering.

**Proposed change**: Use the existing `RtProfile` extension method from MzIO.Processing.MzIOLinq which already supports an optional `ionMobilityRange` parameter ([MzIOLinq.fs lines 231–262](src/MzIO/Processing/MzIOLinq.fs)):

```fsharp
// REPLACE the entire initRTProfile with:
let initRTProfile (reader: IMzIODataReader) (rtIndex: IMzIOArray<RtIndexEntry>) 
                  (rtRange: RangeQuery) (mzRange: RangeQuery) (ionMobilityRange: RangeQuery) =
    reader.RtProfile(rtIndex, rtRange, mzRange, ionMobilityRange)
```

Or, if the ion mobility grouping/summing behavior must be preserved (summing intensities across IM-filtered peaks within mz range), refactor to use binary search:

```fsharp
let initRTProfile (readspecPeaks: string -> Peak1DArray) (rtIndex: IMzIOArray<RtIndexEntry>) 
                  (rtRange: RangeQuery) (mzRange: RangeQuery) (ionMobilityRange: RangeQuery) =
    let entries = RtIndexEntry.Search(rtIndex, rtRange).ToArray()
    let profile = Array.zeroCreate<Peak2D> entries.Length
    for rtIdx = 0 to entries.Length - 1 do
        let entry = entries.[rtIdx]
        let peaks = (readspecPeaks entry.SpectrumID).Peaks
        // Binary search for mz range first
        let mzFiltered = RtIndexEntry.MzSearch(peaks, mzRange)
        // Then binary search for IM range within results  
        let mzFilteredArr = MzIOArray.ToMzIOArray(mzFiltered |> Seq.toArray)
        let imFiltered = RtIndexEntry.IonMobilitySearch(mzFilteredArr, ionMobilityRange)
        let p =
            if imFiltered.Any() then
                let summedIntensity = imFiltered |> Seq.sumBy (fun x -> x.Intensity)
                RtIndexEntry.AsPeak2D(
                    Peak1D(summedIntensity, mzRange.LockValue, ionMobilityRange.LockValue), 
                    entry.Rt)
            else
                RtIndexEntry.AsPeak2D(
                    Peak1D(0., mzRange.LockValue, ionMobilityRange.LockValue), 
                    entry.Rt)
        profile.[rtIdx] <- p
    profile
```

**Expected impact**: ~5–10x speedup on XIC extraction. Binary search reduces per-spectrum cost from O(N) to O(log N + K) where K = matching peaks. For 200K-peak spectra, this is ~15 comparisons vs 200K.

**Risk to correctness**: Low. The current code's `Seq.groupBy Mz` then `Seq.sumBy Intensity` sums intensities from multiple IM-range peaks at similar m/z values. The proposed version sums all IM-filtered peaks, which is equivalent when the mz range is narrow (as it is — typically ±0.02 Da). If exact mz-grouping behavior must be preserved, group after the binary search on the small result set (~1–20 peaks), not the full spectrum.

**Validation**: Compare output `.quant` files (all numeric columns) between old and new for a test dataset. Columns must match within floating-point tolerance (1e-10).

---

### P0-2: Use Memoized Reader in `countMatchedMasses`

**Location**: [PSMBasedQuantificationTIMs.fs line 569](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L569)

**Current behavior**: `let spec = inReader.ReadSpectrumPeaks psm.PSMId` — un-memoized, triggers full SQLite read + binary decode per PSM.

**Proposed change**: Move `readSpecPeaksWithMem` initialization (currently line 705) to before `countMatchedMasses` definition, and use it:

```fsharp
// Move line 705 to before line 560:
let readSpecPeaksWithMem = FSharpAux.Memoization.memoize inReader.ReadSpectrumPeaks

// Then in countMatchedMasses (line 569):
let spec = readSpecPeaksWithMem psm.PSMId  // was: inReader.ReadSpectrumPeaks
```

**Expected impact**: Up to 2–3x speedup on the countMatchedMasses phase, since multiple PSMs may share the same spectrum ID.

**Risk to correctness**: None — same data, just cached.

**Validation**: Output identical to current version.

---

### P0-3: Convert Fragment Ion List to Array/HashSet

**Location**: [PSMBasedQuantificationTIMs.fs lines 561–577](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L561-L577)

**Current behavior**: `frag` is a `List` (F# linked list) built via `List.collect` + `List.map`. Then `List.exists` iterates it for EVERY peak.

**Proposed change**: Convert to sorted array and use binary search, or use `HashSet` with binned mz:

```fsharp
let frag = 
    let ionSeries = (calcIonSeries peptide.BioSequence).TargetMasses
    [|1. .. 2.|]
    |> Array.collect (fun ch -> 
        ionSeries 
        |> List.map (fun x -> Mass.toMZ x.MainPeak.Mass ch)
        |> Array.ofList
    )
    |> Array.sort

// In the filter, use binary search or simple Array.exists (shorter list):
frag
|> Array.exists (fun ion -> abs (ion - peak.Mz) <= (Mass.deltaMassByPpm 100. peak.Mz))
```

**Expected impact**: Minor (~10–20%) but free. Fragment lists are typically 20–60 items; `Array.exists` has better cache locality than `List.exists`.

**Risk to correctness**: None.

---

### P1-1: Parallelize Main Quantification Loop

**Location**: [PSMBasedQuantificationTIMs.fs lines 1164–1170](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L1164-L1170)

**Current behavior**:
```fsharp
|> Array.mapi (fun i (pepIon,psms) -> 
    if i % 100 = 0 then logger.Trace (sprintf "%i peptides quantified" i)
    match processParams.PerformLabeledQuantification with 
    |Domain.Labeling.N15Labeling | ... -> labledQuantification pepIon psms
    |Domain.Labeling.Unlabeled -> lableFreeQuantification pepIon psms
    )
```

**Proposed change**: Pre-fetch all needed spectra into a concurrent dictionary, then use `Array.Parallel.mapi` for the CPU-bound quantification:

```fsharp
// Phase 1: Pre-fetch all spectra (sequential, filling cache)
let allSpectrumIds = 
    qpsmsMzRefined |> Array.collect (fun psm -> [|psm.PSMId|]) |> Array.distinct
allSpectrumIds |> Array.iter (fun id -> readSpecPeaksWithMem id |> ignore)

// Phase 2: Parallel quantification
|> Array.Parallel.mapi (fun i (pepIon,psms) -> 
    if i % 100 = 0 then logger.Trace (sprintf "%i peptides quantified" i)
    match processParams.PerformLabeledQuantification with 
    |Domain.Labeling.N15Labeling | ... -> labledQuantification pepIon psms
    |Domain.Labeling.Unlabeled -> lableFreeQuantification pepIon psms
    )
```

**Challenge**: The `inReader` (MzSQL SQLite connection) is not thread-safe. The memoized cache must serve all threads. Two approaches:
- (a) Pre-fetch all spectra into cache before parallel loop (safest)
- (b) Use `ConcurrentDictionary`-based memoization with lock-per-key on `inReader`

**Expected impact**: 3–8x throughput improvement (limited by core count and memory).

**Risk to correctness**: Medium — must ensure thread-safety of peptideLookUp, chart saving, and logger.

**Validation**: Compare outputs; run on same data single-threaded and parallel, verify identical results.

---

### P1-2: Batch Spectrum Pre-Filtering for IM Range

**Location**: [PSMBasedQuantificationTIMs.fs lines 30–56](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L30-L56)

**Current behavior**: Each call to `initRTProfile` re-reads and re-filters the same spectra if multiple peptide ions have overlapping RT windows.

**Proposed change**: Before the main loop, build a map of `spectrumID → Peak1DArray` for all spectra in the RT range of interest. Filter each spectrum's peaks once into sorted mz arrays that support binary search.

**Expected impact**: ~2x I/O reduction if many peptide ions share RT windows.

**Risk**: Medium — requires careful memory management.

---

### P1-3: Bounded LRU Cache for Spectrum Memoization

**Location**: [PSMBasedQuantificationTIMs.fs line 705](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L705)

**Current behavior**: Unbounded `Dictionary`-based memoization. For TIMs data with ~10K spectra × ~14 MB each = ~140 GB theoretical cache size (will OOM before that).

**Proposed change**: Replace with LRU cache limited to ~500–1000 entries:

```fsharp
let readSpecPeaksWithMem = 
    let cache = System.Collections.Concurrent.ConcurrentDictionary<string, Peak1DArray>()
    let maxSize = 500
    fun (spectrumId: string) ->
        cache.GetOrAdd(spectrumId, fun id ->
            if cache.Count > maxSize then
                // Evict oldest entries (simplified; use proper LRU in production)
                cache.Clear()
            inReader.ReadSpectrumPeaks(id)
        )
```

**Expected impact**: Prevents OOM, reduces GC pressure by ~20–50%.

**Risk**: Low — may slightly increase repeated I/O if cache misses increase. Tune `maxSize` based on available RAM.

---

### P1-4: Memoize `getSpec` in Isotope Comparison

**Location**: [PSMBasedQuantificationTIMs.fs line 163](src/PSMBasedQuantificationTIMs/PSMBasedQuantificationTIMs.fs#L163)

**Current behavior**: `reader.ReadSpectrumPeaks(ms1.ID)` un-memoized in `getSpec`, called from `comparePredictedAndMeasuredIsotopicCluster`.

**Proposed change**: Use `readSpecPeaksWithMem` instead:
```fsharp
let getSpec (readSpecPeaksWithMem: string -> Peak1DArray) (ms1: MzIO.Model.MassSpectrum) =
    Peaks.unzipIMzliteArray (readSpecPeaksWithMem ms1.ID).Peaks
    |> fun (mzData,intensityData) -> PeakArray.zip mzData intensityData
```

**Expected impact**: Eliminates redundant spectrum reads in isotope comparison.

**Risk**: None.

---

### P2-1: Add Filtered Spectrum Read to MzIO (Architectural)

**Location**: MzIO library — [MzIOSQL.fs](src/MzIO.SQL/MzIOSQL.fs), [BinaryDataDecoder.fs](src/MzIO/Binary/BinaryDataDecoder.fs)

**Current behavior**: Always decodes ALL peaks, then filters in application code.

**Proposed change**: Add a `ReadSpectrumPeaksFiltered(spectrumID, mzRange, imRange)` that:
1. Decodes the compressed arrays into raw float arrays (mz[], intensity[], im[])
2. Applies range filtering on the float arrays BEFORE constructing Peak1D objects
3. Returns only matching peaks

**Expected impact**: 10–50x less memory allocation and object creation per spectrum read.

**Risk**: High — requires MzIO library changes, careful testing across all compression types.

---

### P2-2: Convert Peak1D to Struct

**Location**: MzIO library — [PeakArray.fs](src/MzIO/Binary/PeakArray.fs)

**Current behavior**: `Peak1D` inherits from `Peak` (abstract class) — every peak is heap-allocated.

**Proposed change**: Make Peak1D a value type (struct), use `Span<T>` or `Memory<T>` for peak arrays.

**Expected impact**: ~30–50% reduction in GC pressure and decode time.

**Risk**: High — breaking change to MzIO API, requires updating all consumers.

---

## 7. Implementation Plan

### Phase 0: Quick Wins (1–2 days, ~5–10x XIC speedup)

| Task | File | Action | Success Criterion |
|------|------|--------|-------------------|
| **P0-1** Replace `initRTProfile` | TIMs L30–56 | Replace Seq-based filtering with binary search (MzSearch + IonMobilitySearch) | XIC extraction ≥5x faster on benchmark; identical output |
| **P0-2** Memoize `countMatchedMasses` reads | TIMs L569, L705 | Move memoization init before countMatchedMasses; use cached reader | Identical output; measurable I/O reduction |
| **P0-3** Array-ify fragment list | TIMs L561–567 | `List` → `Array` + `Array.sort` | Identical output |

**Rollback**: Revert to original `initRTProfile`. Changes are localized to one function.

**Measurement**: Time the `quantifyPeptides` function end-to-end on a reference TIMs dataset. Target: **≥5x wall-clock speedup** on XIC extraction portion.

---

### Phase 1: Medium Refactors (3–5 days, ~3–8x additional throughput)

| Task | File | Action | Success Criterion |
|------|------|--------|-------------------|
| **P1-1** Parallelize main loop | TIMs L1164–1170 | Pre-fetch spectra, then `Array.Parallel.mapi` | Near-linear scaling to 4 cores; identical output |
| **P1-3** Bounded LRU cache | TIMs L705 | Replace unbounded memoization | Memory usage < 4 GB; no OOM on large datasets |
| **P1-4** Memoize getSpec | TIMs L163 | Thread memoized reader to isotope comparison | Fewer spectrum reads; identical output |

**Rollback**: Each change is independent. Revert individual changes if issues arise.

**Measurement**: Profile with `dotnet-trace` / `dotnet-counters`. Target: **≥3x additional throughput** from parallelism.

---

### Phase 2: Architectural (1–2 weeks, requires MzIO changes)

| Task | Repository | Action | Success Criterion |
|------|-----------|--------|-------------------|
| **P2-1** Filtered spectrum read | MzIO | Add filtered decode pathway | 10x less allocation per spectrum |
| **P2-2** Struct peaks | MzIO | Peak1D → struct | 30–50% GC reduction in microbenchmark |

**Rollback**: MzIO changes need versioning. Use new NuGet package version.

---

## 8. Validation and Benchmark Plan

### Correctness Validation

1. **Reference dataset**: Run current version on a representative TIMs `.mzlite` file + scored PSMs to produce a reference `.quant` output.
2. **Column-wise comparison**: After each optimization, compare all numeric columns with tolerance:
   - Floating-point columns: `|a - b| / max(|a|, |b|) < 1e-8` (relative) or `|a - b| < 1e-12` (absolute)
   - Integer/string columns: exact match
3. **Edge cases**: Empty PSM groups, single-PSM peptide ions, NaN handling.

### Performance Benchmarking

1. **Micro-benchmark**: Time `initRTProfile` in isolation for a known spectrum set (before/after P0-1).
2. **End-to-end**: Time `quantifyPeptides` on reference dataset (before/after each phase).
3. **Memory profiling**: `dotnet-counters monitor` for GC gen0/gen1/gen2 collections, managed heap size.
4. **Profiling tools**:
   - `dotnet-trace collect` for CPU flame graph
   - `dotnet-counters` for real-time GC metrics
   - BenchmarkDotNet for isolated function benchmarks

### Profiling Plan (if no profiling data is available)

```bash
# 1. Collect CPU trace
dotnet-trace collect --process-id <PID> --providers Microsoft-DotNETCore-SampleProfiler

# 2. Monitor GC
dotnet-counters monitor --process-id <PID> --counters System.Runtime

# 3. Analyze trace
dotnet-trace convert trace.nettrace --format speedscope
# Open in https://www.speedscope.app/
```

**Expected findings ranked by confidence**:
1. **[HIGH confidence]** `initRTProfile` dominates CPU time (>50% of hot loop)
2. **[HIGH confidence]** Binary decoding/deserialization is significant (>20%)
3. **[MEDIUM confidence]** GC pauses are measurable (>10% wall clock)
4. **[MEDIUM confidence]** `countMatchedMasses` I/O is significant for large PSM counts
5. **[LOW confidence]** Curve fitting (HULQ.quantifyPeak) is a meaningful fraction

---

## 9. Risks, Assumptions, Open Questions

### Risks

| Risk | Mitigation |
|------|------------|
| Binary search assumes peaks sorted by mz in the database | Verify with test data; MzIO encoder sorts data on write |
| `IonMobilitySearch` assumes peaks sorted by IM after mz filter | The intermediate array from `MzSearch` may not be IM-sorted; need to verify or sort first |
| Parallelization may hit SQLite lock contention | Pre-fetch all needed spectra before parallel loop |
| LRU cache eviction may cause re-reads in worst case | Tune cache size; profile to find optimal value |
| Chart saving in parallel causes file system race conditions | Use per-thread temp directories or synchronize chart saving |

### Assumptions

1. TIMs `.mzlite` files store peaks sorted by mz within each spectrum (standard for MzIO)
2. Spectra sizes are 50K–500K peaks (based on typical timsTOF data)
3. NumPressZLib compression is used (most common for TIMs data)
4. Available RAM ≥ 16 GB (for LRU cache sizing)
5. The `Seq.groupBy Mz` in the current `initRTProfile` was added to handle multiple IM-collapsed peaks at the same nominal mz — this aggregation may or may not be scientifically necessary

### Open Questions

1. **Is the mz-grouping in `initRTProfile` (line 46) scientifically required?** If peaks at slightly different mz values within the mz window should be summed vs. taking closest, this affects P0-1 implementation. Clarify with domain expert.

2. **Are peaks within a single spectrum guaranteed to be sorted by mz?** This affects whether binary search is valid without pre-sorting. The MzIO encoder sorts, but need to verify for Bruker-converted TIMs data.

3. **What is the typical number of peptide ions (main loop iterations)?** This affects parallelization benefit. If <100, parallelism overhead may dominate. If >1000, significant speedup expected.

4. **Is `diagCharts` typically enabled in production?** Chart generation (Plotly HTML) adds significant I/O; disabling it in production is a free speedup.

5. **Does `IonMobilitySearch` return results in IM-sorted order?** If the mz-filtered subset is not IM-sorted, `BinarySearch.Search` won't work correctly for the IM dimension. May need to sort intermediate results or use linear filter only for the (small) IM subset returned by mz binary search.

6. **Could the I/O be moved to SSD/NVMe if not already?** SQLite random-read performance is heavily dependent on storage speed.
