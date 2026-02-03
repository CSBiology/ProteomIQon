# ProteomIQon Unit Testing Implementation Plan

## Overview
Add comprehensive unit tests using Expecto for all ProteomIQon command-line tools. Each tool will have its own test project that includes source files directly (not as NuGet/project references).

## Identified Tools (19 tools + 1 core library)

### Core Library
0. **ProteomIQon** - Core shared library with Domain, DTO, Json, Logging, Core, Drafo, FDRControl', DTW', PepValueCalculation, etc.

### Core Database & Preprocessing Tools
1. **PeptideDB** - Creates peptide databases from FASTA files
2. **MzMLToMzLite** - Converts MzML to MzLite format
3. **MzMLToMzLiteIonMobility** - MzML to MzLite with ion mobility support
4. **Preprocessing** - Data preprocessing tool (Bruker/Wiff/Thermo readers)

### Peptide & Protein Analysis Tools
5. **PeptideSpectrumMatching** - PSM identification with charge state determination, SEQUEST/Andromeda/X!Tandem scoring
6. **ProteinInference** - Infers proteins from peptides with FDR/Q-value calculation
7. **AddDeducedPeptides** - Adds deduced peptides to quant files
8. **PSMStatistics** - PSM statistical analysis with ML-based scoring, Q-values, PEP values

### Quantification Tools
9. **PSMBasedQuantification** - PSM-based quantification (28 tests, +12 functions needed)
10. **PSMBasedQuantificationTIMs** - PSM quantification for TIMS data with ion mobility
11. **AlignmentBasedQuantification** - Alignment-based quantification (15 tests, +9 functions needed)
12. **LabeledProteinQuantification** - Labeled protein quantification (23 tests)
13. **LabelFreeProteinQuantification** - Label-free protein quantification (16 tests)
14. **JoinQuantPepIonsWithProteins** - Joins quantification data with proteins
15. **RatioLFQ** - Ratio-based label-free quantification with constrained least squares

### Statistical & Quality Control Tools
16. **AlignmentBasedQuantStatistics** - Statistics for alignment-based quant with ML scoring
17. **LabelEfficiencyCalculator** - Calculates N15 label efficiency with MIDA isotope simulation

### Utility Tools
18. **MsFraggerToPSM** - Converts MSFragger to PSM format

### Library
19. **QuantBasedAlignment** - Core alignment library with Spline module (not executable)

---

## Implementation Steps

### Phase 0: Setup & Infrastructure
- [x] Analyze repository structure
- [x] 0.1 Create tests root directory structure
- [ ] 0.2 Update solution file to include test projects
- [ ] 0.3 Create shared test utilities (if needed)

### Phase 0.5: Core Library Tests

#### Tool 0: ProteomIQon (Core Library)
- [x] 0.1 Analyze ProteomIQon source files for testable functions
- [x] 0.2 Create `tests/ProteomIQon.Tests/` directory
- [x] 0.3 Create `ProteomIQon.Tests.fsproj` with direct source includes
- [x] 0.4 Write tests for:
  **Json.fs:**
  - [x] **serialize** function (object to JSON string)
  - [x] **serializeAndWrite** function (serialize and write to file)
  - [x] **deserialize** function (JSON string to typed object)
  - [x] **ReadAndDeserialize** function (read file and deserialize)

  **Core.fs - MzIO module:**
  - [x] **Reader.getMzLiteFiles** function (directory file search)
  - [x] **Reader.getThermoRawFiles** function (directory file search)
  - [x] **Reader.getWiffFiles** function (directory file search)
  - [x] **Reader.getBrukerFiles** function (directory file search)
  - [x] **Reader.getMzMLFiles** function (directory file search)
  - [x] **Reader.getMSFilePaths** function (combined file search)
  - [x] **Reader.getReader** function (file extension pattern matching)
  - [ ] **Reader.getDefaultRunID** function (reader type to run ID)
  - [ ] **Processing.changeScanTimeToMinutes** function (time unit conversion)
  - [ ] **Peaks.unzipIMzliteArray** function (peak array unzipping)
  
  **Core.fs - InputPaths module:**
  - [x] **getRelativePath** function (path combination)
  - [x] **parsePath** function (file/directory path parsing)
  - [x] **parsePaths** function (multiple path parsing)
  
  **Core.fs - Zipping module:**
  - [x] **zipDirectory** function (directory compression)
  - [x] **saveZippedDirectory** function (save compressed data)
  
  **BioFSharp.Mz_Temp.fs - SparsePeakArray' module:**
  - [x] **dot** function (sparse vector dot product)
  - [x] **initMzToBinIdx** function (m/z to bin index conversion)
  - [x] **initBinIdxToMz** function (bin index to m/z conversion)
  - [x] **peaksToNearestBinVector** function (peak binning)
  
  **BioFSharp.Mz_Temp.fs - DTW' module:**
  - [x] **warp** function (dynamic time warping)
  - [x] **bestPath** function (optimal warping path extraction)
  - [x] **warping_Path** function (warping path computation)
  - [x] **distance** function (DTW distance calculation)
  - [x] **align** function (sequence alignment)
  - [x] **align'** function (alignment with source mapping)
  - [x] **zNorm** function (z-score normalization)
  
  **BioFSharp.Mz_Temp.fs - SeqIO' module:**
  - [x] **stringFunction** function (value to string formatting)
  - [x] **csv** function (CSV generation)
  
  **BioFSharp.Mz_Temp.fs - FSharpStats'.Wavelet module:**
  - [x] **weightedMean** function (weighted mean calculation)
  - [ ] **identifyPeaksBy** function (wavelet peak identification)
  - [ ] **identifyPeaks** function (peak identification with parameters)
  - [ ] **toIdentifiedPeak** function (peak conversion)
  
  **BioFSharp.Mz_Temp.fs - Fitting' module:**
  - [x] **LogisticFunction** model (logistic regression)
  - [x] **initialParam** function (initial parameter estimation)
  - [x] **initialParamsOverRange** function (parameter estimation over range)
  - [ ] **estimatedParamsWithRSS** function (parameter estimation with RSS)
  
  **BioFSharp.Mz_Temp.fs - SearchDB' module:**
  - [x] **prepareSelectProteinAccessionByID** function (SQL query)
  - [x] **prepareSelectPepSequenceByPepSequenceID** function (SQL query)
  - [x] **prepareSelectMassByModSequenceAndGlobalMod** function (SQL query)
  - [x] **getProteinPeptideLookUpFromFileBy** function (protein-peptide lookup)
  - [x] **selectProteins** function (protein selection)
  - [x] **getSDBParams** function (search DB params retrieval)
  - [x] **setIndexOnModSequenceAndGlobalMod** function (index creation)
  - [x] **getThreadSafePeptideLookUpFromFileBySequenceAndGMod** function (thread-safe lookup)
  
  **BioFSharp.Mz_Temp.fs - ProteinInference' module:**
  - [x] **createInferredProteinClassItemScored** function (record creation)
  - [x] **createInferredProteinClassItemQValue** function (record creation)
  - [x] **createInferredProteinClassItemOut** function (record creation)
  - [x] **removeModification** function (modification removal from sequence)
  - [x] **proteinGroupToString** function (protein group string conversion)
  - [x] **createProteinModelInfoFromEntry** function (model info creation)
  - [x] **assignTranscriptsToGenes** function (transcript-gene assignment)
  - [x] **createPeptideProteinRelation** function (peptide-protein relation)
  - [x] **createPeptideScoreMap** function (peptide score mapping)
  - [x] **createReverseProteinScores** function (reverse protein scoring)
  - [x] **assignPeptideScores** function (peptide score assignment)
  - [x] **assignDecoyScoreToTargetScore** function (decoy score assignment)
  
  **BioFSharp.Mz_Temp.fs - FDRControl' module:**
  - [x] **getLogisticRegressionFunction** function (logistic regression)
  - [x] **binningFunction** function (score binning)
  - [x] **estimatePi0HG** function (pi0 estimation with hypergeometric)
  - [x] **binProteinsLength** function (protein length binning)
  - [x] **expectedFP** function (expected false positives)
  - [x] **calculateFDRwithMAYU** function (MAYU FDR calculation)
  - [x] **calculateFDRwithDecoyTargetRatio** function (decoy/target FDR)
  - [x] **calculateQValueLogReg** function (Q-value with logistic regression)
  - [x] **calculateQValueStorey** function (Q-value with Storey method)
  - [x] **assignQValueToIPCIS** function (Q-value assignment)
  - [x] **createTargetDecoyHis** function (target/decoy histogram)
  - [x] **calculatePEPValues** function (PEP value calculation)
  - [x] **logitTransformPepValues** function (logit transformation)
  - [x] **initCalculateLin** function (linear PEP calculation)
  
  **BioFSharp.Mz_Temp.fs - Fragmentation' module:**
  - [x] **LadderedTaggedMass** class (mass ladder representation)
  - [x] **createLadderedPeakFamily** function (peak family creation)
  - [x] **ladderAndChargeElement** function (ladder with charge states)
  - [x] **ladderElement** function (ion ladder creation)
  
  **BioFSharp.Mz_Temp.fs - Drafo.Core module:**
  - [x] **Key** class (DynamicObj-based key with equality/hash)
  - [x] **indexWithColumnValues** function (Frame indexing with Key)
  - [x] **readFrame** function (CSV to Frame reading)
  - [x] **readAndIndexFrame** function (read and index Frame)
  - [x] **getColumn** function (typed column extraction)
  - [x] **seriesToFrame** function (Series to Frame conversion)
  - [x] **rowKeyToColumns** function (row keys to columns)
  - [x] **createFilter** function (boolean filter Series creation)
  - [x] **transform** function (Series value transformation)
  - [x] **zip** function (Series zipping with operation)
  - [x] **dropAllKeyColumnsBut** function (keep specified key columns)
  - [x] **dropKeyColumns** function (remove key columns)
  - [x] **groupTransform** function (grouped transformation)
  - [x] **createGroupFilter** function (group-based filter)
  - [x] **aggregate** function (aggregation with filters)
  - [x] **assemble** function (column assembly to Frame)
  - [x] **pivot** function (Frame pivoting)
  
  **pepValueCalculation.fs:**
  - [ ] **initCalculatePEPValueIRLS** function (IRLS-based PEP calculation)
  - [ ] **invlogit** function (inverse logit transformation)
  - [ ] **logit** function (logit transformation)
  - [ ] **splineEval** function (spline evaluation)
  - [ ] **iterativeReweightedLeastSquares** function (IRLS algorithm)
  
  **DeedleExtensionsTemp.fs:**
  - [ ] **distinctRowValues** function (distinct row value selection)
  - [ ] **append** function (Series append)
  - [ ] Private helper functions for Frame operations
- [ ] 0.5 Run and verify tests pass

### Phase 0.6: Test Quality Audit & Fixes (Existing Tests)
- [ ] 0.6.1 Replace placeholder assertions (`Expect.isTrue true`) in `tests/ProteomIQon.Tests/Tests.fs` with real fixtures and verifications for:
  - `Core.MzIO.Processing.changeScanTimeToMinutes`
  - `Core.MzIO.Peaks.unzipIMzliteArray`
  - `FSharpStats'.Wavelet.identifyPeaksBy`, `identifyPeaks`, `toIdentifiedPeak`
  - `Fitting'.NonLinearRegression'.LevenbergMarquardtConstrained'.estimatedParamsWithRSS`
  - `SearchDB'` prepared queries and DB lookups (use in-memory SQLite)
  - `ProteinInference'`, `FDRControl'`, `Fragmentation'` functions currently stubbed
  - `Drafo.Core.indexWithColumnValues`
- [ ] 0.6.2 Fix tautological/non-validating assertions:
  - `tests/PSMBasedQuantification.Tests/Tests.fs` "closest RT" test uses `offset = 10.0 || offset = 10.0`
  - `tests/ProteomIQon.Tests/Tests.fs` `DTW'.bestPath` test uses `List.length result >= 0`
  - `tests/ProteomIQon.Tests/Tests.fs` `Core.MzIO.Reader.getDefaultRunID` does not call the function
- [ ] 0.6.3 Convert Labeled/LabelFree quant tests from local math checks to actual module calls:
  - `LabeledProteinQuantification.performGlobalModAggregation`, `performChargeOrModAggregation`, `performPeptideToProteinAggregation`
  - `LabelFreeProteinQuantification.performChargeOrModAggregation`, `performPeptideToProteinAggregation`
- [ ] 0.6.4 Add edge-case coverage where currently missing:
  - `klDiv` zero/invalid inputs and length mismatch handling
  - `weightedMean` zero total weight and length mismatch
  - `light/heavyQualityFilter` all-NaN medians, zero medians, out-of-range rejection

### Phase 1: Core Tools

#### Tool 1: PeptideDB
- [ ] 1.1 Analyze PeptideDB.fs for testable functions
- [ ] 1.2 Create `tests/PeptideDB.Tests/` directory
- [ ] 1.3 Create `PeptideDB.Tests.fsproj` with direct source includes
- [ ] 1.4 Write tests for:
  - [ ] **SearchDbParams record creation** (validate field mapping from PeptideDBParams)
  - [ ] **SearchDB.connectOrCreateDB** usage patterns (mock DB connection)
  - [ ] **SearchDB'.setIndexOnModSequenceAndGlobalMod** usage
  - [ ] Parameter validation logic (MinMissedCleavages <= MaxMissedCleavages, etc.)
- [ ] 1.5 Run and verify tests pass

#### Tool 2: MzMLToMzLite
- [ ] 2.1 Analyze MzMLToMzLite.fs for testable functions
- [ ] 2.2 Create `tests/MzMLToMzLite.Tests/` directory
- [ ] 2.3 Create `MzMLToMzLite.Tests.fsproj` with direct source includes
- [ ] 2.4 Write tests for:
  - [ ] **initPeakPicking** function logic (PeakPicking DU pattern matching)
  - [ ] **SignalDetection.Padding.createPaddingParameters** usage
  - [ ] **SignalDetection.Wavelet.createWaveletParameters** usage
  - [ ] **SignalDetection.Wavelet.toCentroidWithRicker2D** usage
  - [ ] **insertSprectrum** function (MS level routing logic)
  - [ ] **PeakArray.createPeak1DArray** usage
  - [ ] YThreshold.Fixed vs YThreshold.MinSpectrumIntensity logic
  - [ ] CentroidizationMode pattern matching
- [ ] 2.5 Run and verify tests pass

#### Tool 3: MzMLToMzLiteIonMobility
- [ ] 3.1 Analyze MzMLToMzLiteIonMobility.fs for testable functions
- [ ] 3.2 Create `tests/MzMLToMzLiteIonMobility.Tests/` directory
- [ ] 3.3 Create `MzMLToMzLiteIonMobility.Tests.fsproj` with direct source includes
- [ ] 3.4 Write tests for:
  - [ ] **fixSpectrum** function (fix null Precursors/Scans/Products lists)
  - [ ] **createPeak1DArrayCopy** function (copy peak array with same compression settings)
  - [ ] **toIMPeaks** function (convert peaks with ion mobility data)
  - [ ] **unzipIonMobilityMzliteArray** function (unzip to mz/intensity/ionMobility arrays)
  - [ ] **createPeak1DArray** function (create Peak1DArray from mz/intensity/ionMobility arrays)
  - [ ] **getDefaultRunID** function (get default run ID from reader type)
  - [ ] **initPeakPicking** function (wavelet/profile peak picking logic)
  - [ ] Peak1D record creation with ion mobility
- [ ] 3.5 Run and verify tests pass

#### Tool 4: Preprocessing
- [ ] 4.1 Analyze Preprocessing.fs for testable functions
- [ ] 4.2 Create `tests/Preprocessing.Tests/` directory
- [ ] 4.3 Create `Preprocessing.Tests.fsproj` with direct source includes
- [ ] 4.4 Write tests for:
  - [ ] **getReader** function (file extension pattern matching for reader types)
  - [ ] **getDefaultRunID** function (reader type to run ID mapping)
  - [ ] **initPeakPicking** function (reader/PeakPicking DU pattern matching)
  - [ ] **SignalDetection.Padding.createPaddingParameters** usage
  - [ ] **SignalDetection.Wavelet.createWaveletParameters** usage
  - [ ] **SignalDetection.Wavelet.toCentroidWithRicker2D** usage
  - [ ] **insertSprectrum** function (MS level routing logic)
  - [ ] **PeakArray.createPeak1DArray** usage
  - [ ] File extension detection (.wiff, .d, .mzlite, .raw, .RAW)
  - [ ] BafFileReader manufacturer peak picking branch
- [ ] 4.5 Run and verify tests pass

#### Tool 5: PeptideSpectrumMatching
- [ ] 5.1 Analyze PeptideSpectrumMatching.fs for testable functions
- [ ] 5.2 Create `tests/PeptideSpectrumMatching.Tests/` directory
- [ ] 5.3 Create `PeptideSpectrumMatching.Tests.fsproj` with direct source includes
- [ ] 5.4 Write tests for:
  - [ ] **getPrecursorCharge** function (charge state determination with position metric)
  - [ ] **ChargeState.putativePrecursorChargeStatesBy** usage
  - [ ] **ChargeState.peakPosStdDevBy** usage
  - [ ] **ChargeState.initMzDevOfRndSpec** usage
  - [ ] **ChargeState.empiricalPValueOfSim** usage
  - [ ] **ChargeState.removeSubSetsOfBestHit** usage
  - [ ] **SequestLike.getTheoSpecs** usage
  - [ ] **SequestLike.calcSequestScore** usage
  - [ ] **XScoring.getTheoSpecs** usage
  - [ ] **XScoring.calcAndromedaAndXTandemScore** usage
  - [ ] **Fragmentation.Series.fragmentMasses** usage
  - [ ] Score normalization (NormDeltaBestToRest, NormDeltaNext)
  - [ ] Label assignment (+1 target, -1 decoy)
- [ ] 5.5 Run and verify tests pass

#### Tool 6: ProteinInference
- [ ] 6.1 Analyze ProteinInference.fs for testable functions
- [ ] 6.2 Create `tests/ProteinInference.Tests/` directory
- [ ] 6.3 Create `ProteinInference.Tests.fsproj` with direct source includes
- [ ] 6.4 Write tests for:
  - [ ] **initFDR** function (FDR calculation with MAYU/DecoyTargetRatio/Conservative methods)
  - [ ] **initQValue** function (Q-value calculation with Storey/LogisticRegression/NoQValue methods)
  - [ ] **qValueHitsVisualization** function (histogram/chart generation - test data preparation)
  - [ ] **createClassItemCollection** function (protein-peptide classification map creation)
  - [ ] **readAndInferFile** function (protein inference pipeline)
  - [ ] **ClassInfo record creation** tests (Sequence, Class, Proteins fields)
  - [ ] **FDRControl'.calculateFDRwithMAYU** usage
  - [ ] **FDRControl'.calculateFDRwithDecoyTargetRatio** usage
  - [ ] **FDRControl'.calculateQValueStorey** usage
  - [ ] **FDRControl'.calculateQValueLogReg** usage
  - [ ] **FDRControl'.assignQValueToIPCIS** usage
  - [ ] **ProteinInference'.createPeptideProteinRelation** usage
  - [ ] **ProteinInference'.createPeptideScoreMap** usage
  - [ ] **ProteinInference'.createReverseProteinScores** usage
  - [ ] **ProteinInference'.assignPeptideScores** usage
  - [ ] **ProteinInference'.assignDecoyScoreToTargetScore** usage
  - [ ] **ProteinInference'.createInferredProteinClassItemScored** usage
  - [ ] **PeptideClassification.classify** usage
  - [ ] **PeptideClassification.createLocusSpliceVariantCount** usage
- [ ] 6.5 Run and verify tests pass

#### Tool 7: AddDeducedPeptides
- [ ] 7.1 Analyze AddDeducedPeptides.fs for testable functions
- [ ] 7.2 Create `tests/AddDeducedPeptides.Tests/` directory
- [ ] 7.3 Create `AddDeducedPeptides.Tests.fsproj` with direct source includes
- [ ] 7.4 Write tests for:
  - [ ] **ensureEqualProtDetermination** logic (validate protein inference consistency)
  - [ ] **pepProt expansion** logic (expand proteins by peptides, distinct by pep/prot pairs)
  - [ ] **String.split** peptide sequence parsing (semicolon-separated)
  - [ ] **presentPeptides Set creation** (peptide filtering by presence)
  - [ ] **protInfResultGrouped** grouping logic (group by ProteinGroup)
  - [ ] **protInfResultsCombined** combining logic (new PeptideSequence field creation)
- [ ] 7.5 Run and verify tests pass

#### Tool 8: PSMStatistics
- [ ] 8.1 Analyze PSMStatistics.fs for testable functions
- [ ] 8.2 Create `tests/PSMStatistics.Tests/` directory
- [ ] 8.3 Create `PSMStatistics.Tests.fsproj` with direct source includes
- [ ] 8.4 Write tests for:
  - [ ] **downcastPipeline** function (ML pipeline type conversion)
  - [ ] **restorePSMID** function (PSM ID string restoration from underscores/dashes)
  - [ ] **getFlankingAminoAcids** function (extract N/C terminal flanking amino acids from protein)
  - [ ] **initToPSMToLearn** function (convert PSM to PSMToLearn record with flanking AA, proteins, miss cleavages)
  - [ ] **initProteinAndClvIdxLookUp** function (lookup proteins and cleavage indices by PepSequenceID)
  - [ ] **PSMToLearn record creation** tests
  - [ ] **PSMPrediction record creation** tests
  - [ ] **AppliedModel record creation** tests
  - [ ] **calculateQValueStorey** usage (via BioFSharp.Mz.FDRControl)
  - [ ] **initCalculatePEPValueIRLS** usage (via ProteomIQon.PepValueCalculation)
  - [ ] **ML pipeline training** logic (trainModel function)
  - [ ] **applyModel** function (scoring and Q-value calculation)
  - [ ] **Fixed threshold filtering** logic (Threshold.Fixed branch)
- [ ] 8.5 Run and verify tests pass

### Phase 2: Quantification Tools

#### Tool 9: PSMBasedQuantification
- [x] 9.1 Analyze PSMBasedQuantification.fs for testable functions
- [x] 9.2 Create `tests/PSMBasedQuantification.Tests/` directory
- [x] 9.3 Create `PSMBasedQuantification.Tests.fsproj` with direct source includes
- [ ] 9.4 Write tests for:
  - [x] klDiv function (5 tests)
  - [x] weightedMean function (4 tests)
  - [x] createAveragePSM function (4 tests)
  - [x] getBaseLineCorrectionOffsetAt function (4 tests)
  - [x] getInferredXic function (6 tests)
  - [x] substractBaseLine function (5 tests)
  - [ ] **MISSING**: Query.initRTProfile function (RT profile extraction)
  - [ ] **MISSING**: initSpline function (spline fitting with binning)
  - [ ] **MISSING**: searchRTMinusFittedRtTarget function (RT difference calculation)
  - [ ] **MISSING**: searchRTMinusFittedRtInferred function (RT difference for inferred)
  - [ ] **MISSING**: chooseScanTime function (scan time selection logic)
  - [ ] **MISSING**: optimizeWindowWidth function (window optimization via autocorrelation)
  - [ ] **MISSING**: initGetWindowWidth function (window width initialization)
  - [ ] **MISSING**: calcCorrelation function (peak correlation calculation)
  - [ ] **MISSING**: generateIsotopicDistributionOfFormulaBySum function (isotope distribution)
  - [ ] **MISSING**: lightQualityFilter function (quality filtering)
  - [ ] **MISSING**: heavyQualityFilter function (quality filtering)
  - [ ] **MISSING**: average function (PSM averaging with weighted scan time)
- [ ] 9.5 Run and verify tests pass (**28 tests passing**)

#### Tool 10: PSMBasedQuantificationTIMs
- [ ] 10.1 Analyze PSMBasedQuantificationTIMs.fs for testable functions
- [ ] 10.2 Create `tests/PSMBasedQuantificationTIMs.Tests/` directory
- [ ] 10.3 Create `PSMBasedQuantificationTIMs.Tests.fsproj` with direct source includes
- [ ] 10.4 Write tests for:
  - [ ] **Query.initRTProfile** function (RT profile extraction with ion mobility support)
  - [ ] **createAveragePSM** function (record creation with ion mobility)
  - [ ] **initSpline** function (spline fitting with binning and lambda optimization)
  - [ ] **getBaseLineCorrectionOffsetAt** function (baseline offset calculation)
  - [ ] **getClosestMs1** function (find closest MS1 spectrum by scan time)
  - [ ] **getSpec** function (spectrum retrieval and unzipping)
  - [ ] **lightQualityFilter** function (quality filter for light peptides)
  - [ ] **heavyQualityFilter** function (quality filter for heavy peptides)
  - [ ] **klDiv** function (Kullback-Leibler divergence)
  - [ ] **substractBaseLine** function (ALS baseline subtraction)
  - [ ] **initGetProcessedXIC** function (XIC processing with ion mobility and baseline)
  - [ ] **initGetIsotopicEnvelope** function (isotopic envelope with ion mobility)
  - [ ] **weightedMean** function (weighted mean calculation)
  - [ ] **average** function (PSM averaging with weighted scan time and ion mobility)
- [ ] 10.5 Run and verify tests pass

#### Tool 11: AlignmentBasedQuantification
- [x] 11.1 Analyze AlignmentBasedQuantification.fs for testable functions
- [x] 11.2 Create `tests/AlignmentBasedQuantification.Tests/` directory
- [x] 11.3 Create `AlignmentBasedQuantification.Tests.fsproj` with direct source includes
- [ ] 11.4 Write tests for:
  - [x] klDiv function (5 tests)
  - [x] getBaseLineCorrectionOffsetAt function (4 tests)
  - [x] getInferredXic function (2 tests)
  - [x] lightQualityFilter function (2 tests)
  - [x] heavyQualityFilter function (1 test)
  - [x] scanTimeDifferenceFilter function (1 test)
  - [ ] **MISSING**: Query.initRTProfile function (RT profile extraction with ion mobility)
  - [ ] **MISSING**: initSpline function (spline fitting with binning)
  - [ ] **MISSING**: substractBaseLine function (baseline subtraction with ALS)
  - [ ] **MISSING**: lightScanTimeDifferenceFilter function (scan time difference filter)
  - [ ] **MISSING**: heavyScanTimeDifferenceFilter function (scan time difference filter)
  - [ ] **MISSING**: initGetProcessedXIC function (XIC processing with baseline)
  - [ ] **MISSING**: getRefinedXic function (XIC refinement with DTW alignment)
  - [ ] **MISSING**: getClosestMs1 function (find closest MS1 spectrum)
  - [ ] **MISSING**: getSpec function (spectrum retrieval and unzipping)
- [ ] 11.5 Run and verify tests pass (**15 tests passing**)

#### Tool 12: LabeledProteinQuantification
- [x] 12.1 Analyze LabeledProteinQuantification.fs for testable functions
- [x] 12.2 Create `tests/LabeledProteinQuantification.Tests/` directory
- [x] 12.3 Create `LabeledProteinQuantification.Tests.fsproj` with direct source includes
- [ ] 12.4 Write tests for:
  - [x] tryGetQVal function (5 tests)
  - [x] sem function (6 tests)
  - [x] performGlobalModAggregation function (2 tests)
  - [x] performChargeOrModAggregation function (2 tests)
  - [x] performPeptideToProteinAggregation function (8 tests)
  - [ ] **MISSING**: Seq.cv function usage (coefficient of variation)
  - [ ] **MISSING**: Seq.stDevBy function usage (standard deviation)
  - [ ] **MISSING**: correlation_Light_Heavy_Threshold filter logic
  - [ ] **MISSING**: alignment_QValue filter logic
  - [ ] **MISSING**: modPepFilter function (modification peptide filtering)
  - [ ] **MISSING**: keyCols array construction
  - [ ] **MISSING**: Frame.sliceCols usage
  - [ ] **MISSING**: Frame.mergeAll usage
  - [ ] **MISSING**: labeledQuantification pipeline integration tests
- [ ] 12.5 Run and verify tests pass (**23 tests passing**)

**Note**: Core.* functions (getColumn, zip, transform, createFilter, createGroupFilter, aggregate, assemble, dropKeyColumns, dropAllKeyColumnsBut, indexWithColumnValues, pivot, rowKeyToColumns) are tested in ProteomIQon.Tests (Tool 0)

#### Tool 13: LabelFreeProteinQuantification
- [x] 13.1 Analyze LabelFreeProteinQuantification.fs for testable functions
- [x] 13.2 Create `tests/LabelFreeProteinQuantification.Tests/` directory
- [x] 13.3 Create `LabelFreeProteinQuantification.Tests.fsproj` with direct source includes
- [ ] 13.4 Write tests for:
  - [x] tryGetQVal function (5 tests)
  - [x] sem function (5 tests)
  - [x] performChargeOrModAggregation function (2 tests)
  - [x] performPeptideToProteinAggregation function (4 tests)
  - [ ] **MISSING**: Seq.cv function usage (coefficient of variation)
  - [ ] **MISSING**: Seq.stDevBy function usage (standard deviation)
  - [ ] **MISSING**: alignment_QValue filter logic (Frame.filterRows)
  - [ ] **MISSING**: modPepFilter function (modification peptide filtering)
  - [ ] **MISSING**: keyCols array construction
  - [ ] **MISSING**: Frame.sliceCols usage
  - [ ] **MISSING**: Frame.mergeAll usage
  - [ ] **MISSING**: labelFreeQuantification pipeline integration tests
- [ ] 13.5 Run and verify tests pass (**16 tests passing**)

**Note**: Core.* functions (getColumn, transform, createFilter, createGroupFilter, aggregate, assemble, dropKeyColumns, dropAllKeyColumnsBut, indexWithColumnValues, pivot, rowKeyToColumns) are tested in ProteomIQon.Tests (Tool 0)

#### Tool 14: JoinQuantPepIonsWithProteins
- [ ] 14.1 Analyze JoinQuantPepIonsWithProteins.fs for testable functions
- [ ] 14.2 Create `tests/JoinQuantPepIonsWithProteins.Tests/` directory
- [ ] 14.3 Create `JoinQuantPepIonsWithProteins.Tests.fsproj` with direct source includes
- [ ] 14.4 Write tests for:
  - [ ] **Frame.expandRowsByColumn** logic (expand PeptideSequence by semicolon split)
  - [ ] **Frame.align** logic (align proteins with peptide ions by StringSequence)
  - [ ] **Frame.indexRowsUsing** logic (create composite key from multiple columns)
  - [ ] **Frame.mapRowKeys** logic (create FileName/ProteinGroup/StringSequence/etc. key)
  - [ ] **Frame.mapRows** result record creation (QuantAndProtResult mapping)
  - [ ] **QuantificationResult.tryTargetGetScanTime** filter logic
  - [ ] Semicolon-separated PeptideSequence splitting
  - [ ] QValue column renaming to ProteinGroup_QValue
- [ ] 14.5 Run and verify tests pass

#### Tool 15: RatioLFQ
- [ ] 15.1 Analyze RatioLFQ.fs for testable functions
- [ ] 15.2 Create `tests/RatioLFQ.Tests/` directory
- [ ] 15.3 Create `RatioLFQ.Tests.fsproj` with direct source includes
- [ ] 15.4 Write tests for:
  - [ ] **equalityConstrainedLS** function (constrained least squares solver)
  - [ ] **computeLFQProt** function (LFQ protein computation with ratio matrix)
  - [ ] **readQuantResult** function (read QuantificationResult to Deedle Frame)
  - [ ] **combinedRatio** function (combined ratio calculation with log2 transformation)
  - [ ] **combinedOriginal** function (combined original intensity calculation)
  - [ ] **lfq** function (main LFQ computation pipeline)
  - [ ] Log2 transformation edge cases (zeros, negatives)
  - [ ] Ratio matrix construction from intensities
  - [ ] Least squares with equality constraints
- [ ] 15.5 Run and verify tests pass

### Phase 3: Statistical & QC Tools

#### Tool 16: AlignmentBasedQuantStatistics
- [ ] 16.1 Analyze AlignmentBasedQuantStatistics.fs for testable functions
- [ ] 16.2 Create `tests/AlignmentBasedQuantStatistics.Tests/` directory
- [ ] 16.3 Create `AlignmentBasedQuantStatistics.Tests.fsproj` with direct source includes
- [ ] 16.4 Write tests for:
  - [ ] **downcastPipeline** function (ML pipeline type conversion)
  - [ ] **toPeptideForLearning** function (convert ObjectSeries to PeptideForLearning with DTW distance)
  - [ ] **createDataToScore** function (create scoring data from quant/align files)
  - [ ] **createTrainingsData** function (create ML training data with positive/negative labels)
  - [ ] **PeptideForLearning record creation** tests
  - [ ] **ClassPredition record creation** tests
  - [ ] **DTW distance calculation** (via DTW'.distance)
  - [ ] **zNorm** function usage (via DTW'.zNorm)
  - [ ] **assignScoreAndQValue** function (ML scoring pipeline)
  - [ ] **calculateQValueStorey** usage (Q-value calculation)
- [ ] 16.5 Run and verify tests pass

#### Tool 17: LabelEfficiencyCalculator
- [ ] 17.1 Analyze LabelEfficiencyCalculator.fs for testable functions
- [ ] 17.2 Create `tests/LabelEfficiencyCalculator.Tests/` directory
- [ ] 17.3 Create `LabelEfficiencyCalculator.Tests.fsproj` with direct source includes
- [ ] 17.4 Write tests for:
  - [ ] **label** function (apply N15 labeling efficiency to formula)
  - [ ] **generateIsotopicDistribution** function (MIDA isotope distribution)
  - [ ] **floatArrayOf** function (parse semicolon-separated float string to array)
  - [ ] **getHeavyPattern** function (extract heavy isotope pattern from mz/intensity strings)
  - [ ] **simulateFrom** function (simulate isotope pattern for peptide/charge/efficiency)
  - [ ] **normBySum** function (normalize mz/intensity pairs by sum)
  - [ ] **klDiv** function (Kullback-Leibler divergence calculation)
  - [ ] **compareIsotopicDistributions** function (compare measured vs simulated patterns)
  - [ ] **calcKL** function (calculate KL for extracted pattern vs simulation)
  - [ ] **PeptideIon.create** function (record creation)
  - [ ] **ExtractedIsoPattern.create** function (record creation)
  - [ ] **SimulatedIsoPattern.create** function (record creation)
  - [ ] **Brent optimization** for label efficiency (via FSharp.Stats.Optimization.Brent.minimize)
- [ ] 17.5 Run and verify tests pass

### Phase 4: Utility Tools

#### Tool 18: MsFraggerToPSM
- [ ] 18.1 Analyze MsFraggerToPSM.fs for testable functions
- [ ] 18.2 Create `tests/MsFraggerToPSM.Tests/` directory
- [ ] 18.3 Create `MsFraggerToPSM.Tests.fsproj` with direct source includes
- [ ] 18.4 Write tests for:
  - [ ] **readFragegrPsms** function (read MSFragger tab-separated output to PSM records)
  - [ ] **initModSeqLookup** function (initialize mod sequence lookup from DB)
  - [ ] **prepareSelectModsequenceBySequence** function (SQL query preparation)
  - [ ] **findClosest** function (find closest spectrum by scan time and m/z)
  - [ ] **Modification parsing logic** (N-term acetylation, oxidation M, carbamidomethyl C)
  - [ ] **MSFragger column mapping** (Assigned Modifications, Peptide, Retention, etc.)
  - [ ] PSMStatisticsResultFragpipe record creation tests
- [ ] 18.5 Run and verify tests pass

### Phase 5: Library

#### Tool 19: QuantBasedAlignment (Library)
- [ ] 19.1 Analyze QuantBasedAlignment.fs for testable functions
- [ ] 19.2 Create `tests/QuantBasedAlignment.Tests/` directory
- [ ] 19.3 Create `QuantBasedAlignment.Tests.fsproj` with direct source includes
- [ ] 19.4 Write tests for:
  - [ ] **Spline module functions**:
    - [ ] **ReinschQ type** (indexer, constructor tests)
    - [ ] **multiplyInplace** function (matrix-vector multiplication)
    - [ ] **multiplyInplaceTransposedQ** function (transposed Q multiplication)
    - [ ] **ReinschR type** (indexer, constructor tests)
    - [ ] **QtQpR** function (Q^T * Q + R matrix computation)
    - [ ] **diff** function (array differentiation)
    - [ ] **initFullSparseRpαQtQ** function (sparse matrix initialization)
    - [ ] **solveLeastSquares** function (LAPACK-based least squares)
    - [ ] **fitSplineSparseLapack** function (spline fitting with sparse matrices)
    - [ ] **initPredict** function (prediction function initialization)
    - [ ] **fitSplineSparseLapack'** function (alternative spline fitting)
  - [ ] **Main module functions**:
    - [ ] **getRandomColor** function (random color selection)
    - [ ] **toPeptideIon** function (QuantificationResult to PeptideIon conversion)
    - [ ] **toPeptideForLearning** function (learning record with source/target info)
    - [ ] **getQuantifiedPeptides** function (read and filter quantified peptides)
    - [ ] **createAlignmentFiles** function (alignment file structure creation)
    - [ ] **createAlignmentResult** function (alignment result record creation)
- [ ] 19.5 Run and verify tests pass

### Phase 6: Integration & CI/CD
- [ ] 20.1 Update solution file with all test projects
- [ ] 20.2 Create root test runner (if needed)
- [ ] 20.3 Verify `dotnet test` runs all test projects
- [ ] 20.4 Update CI/CD configuration (if exists)
- [ ] 20.5 Add test coverage reporting (optional)
- [ ] 20.6 Document testing approach in README

---

## Testing Strategy

### For Each Tool:

1. **Identify Testable Functions**
   - Pure functions (no I/O, deterministic)
   - Parsing/validation logic
   - Mathematical calculations
   - Data transformations
   - Algorithm components

2. **Test Categories**
   - **Normal Cases**: Standard inputs, expected workflows
   - **Edge Cases**: Empty inputs, boundary values, extreme values
   - **Error Handling**: Invalid inputs, malformed data, exceptional conditions
   - **Determinism**: Same input produces same output
   - **Invariants**: Properties that should always hold

3. **CLI Argument Parsing**
   - Each tool has `CLIArgumentParsing.fs`
   - Test valid argument combinations
   - Test missing required arguments
   - Test invalid argument values
   - Test help text generation

4. **Minimal Refactoring**
   - Only refactor if absolutely necessary for testability
   - Keep changes minimal and well-documented
   - Preserve existing behavior exactly
   - Document any behavior changes

### Expecto Test Pattern (from existing Sample.fs):

```fsharp
module ToolName.Tests

open Expecto
open ProteomIQon

[<Tests>]
let tests =
  testList "ToolName" [
    testList "FunctionName" [
      testCase "normal case" <| fun _ ->
        let result = ToolName.functionName validInput
        Expect.equal result expectedOutput "should produce correct output"
      
      testCase "edge case: empty input" <| fun _ ->
        let result = ToolName.functionName []
        Expect.isEmpty result "should handle empty input"
      
      testCase "error case: invalid input" <| fun _ ->
        Expect.throws (fun () -> ToolName.functionName invalidInput |> ignore)
          "should throw on invalid input"
    ]
  ]
```

### Project File Template:

```xml
<Project Sdk="Microsoft.NET.Sdk">
  <PropertyGroup>
    <OutputType>Exe</OutputType>
    <TargetFramework>net8.0</TargetFramework>
    <GenerateProgramFile>false</GenerateProgramFile>
  </PropertyGroup>

  <ItemGroup>
    <!-- Include tool source files directly -->
    <Compile Include="..\..\src\ToolName\ToolName.fs" />
    <Compile Include="..\..\src\ToolName\CLIArgumentParsing.fs" />
    <!-- Add other .fs files as needed, in correct order -->
    
    <!-- Test files last -->
    <Compile Include="Tests.fs" />
  </ItemGroup>

  <ItemGroup>
    <!-- Reference ProteomIQon core library -->
    <ProjectReference Include="..\..\src\ProteomIQon\ProteomIQon.fsproj" />
    
    <!-- Expecto test packages -->
    <PackageReference Include="Expecto" Version="9.*" />
    <PackageReference Include="YoloDev.Expecto.TestSdk" Version="0.*" />
    <PackageReference Include="Microsoft.NET.Test.Sdk" Version="17.*" />
    <PackageReference Update="FSharp.Core" Version="6.*" />
    
    <!-- Add other dependencies from tool's .fsproj -->
  </ItemGroup>
</Project>
```

---

## Notes

- **ProteomIQon Core**: Most tools reference `ProteomIQon.fsproj` which contains shared types (Domain, Core, DTO, Drafo, FDRControl', DTW', etc.). The core library has its own test project (Tool 0).
- **Dependencies**: Copy relevant PackageReferences from each tool's .fsproj to its test project.
- **Order Matters**: F# compilation order is strict. Ensure source files are included in the correct dependency order.
- **Shared Utilities**: If multiple tests need similar helpers (e.g., mock data generation), consider creating a shared test utilities module.

---

## Success Criteria

- [ ] All 20 test projects created (19 tools + 1 core library)
- [ ] Each test project includes source files directly (no project/NuGet references to tools)
- [ ] Tests cover core functionality, edge cases, and error handling
- [ ] All tests pass with `dotnet test`
- [ ] Solution file includes all test projects
- [ ] CI/CD runs tests automatically
- [ ] Test coverage is reasonable (aim for >70% on testable functions)
