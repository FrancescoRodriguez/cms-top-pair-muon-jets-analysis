# Top-pair cross-section measurement with CMS Open Data

This repository contains a ROOT-based analysis of CMS 2016 Open Data in the
single-muon plus jets final state. The final goal is to measure the top-quark
pair production cross section and study a differential distribution with an
unfolding procedure.

## Reference publication

This analysis is based on the CMS measurement described in the following publication:

**Reference DOI:** [10.1103/PhysRevD.95.092001](https://doi.org/10.1103/PhysRevD.95.092001)

The reference publication provides the physics methodology for the event selection, top-quark reconstruction, cross-section measurement, and unfolding procedure. The implementation in this repository adapts these methods to the available CMS 2016 Open Data and includes simplifications documented in the sections below.

## Phase 1: NanoAOD inspection

The first macro checks the structure of a remote NanoAOD file and creates a
small set of diagnostic plots without downloading the complete dataset.

```bash
mkdir -p output/phase1
root -l -b -q 'macros/inspectNanoAOD.C("root://.../file.root","sample",20000)'
```

The arguments are the file URL, a label used for output filenames, and the
maximum number of events inspected. The macro writes:

- a text report containing tree sizes and required-branch availability;
- a ROOT file containing the diagnostic histograms;
- a PNG summary with muon, jet, MET, and generator-level information.

The URLs selected for the initial inspection are listed in
[`config/phase1_samples.txt`](config/phase1_samples.txt).

To inspect every configured sample using 1000 events per file:

```bash
./scripts/run_phase1.sh 1000
```

## Input datasets

| Sample | CERN Open Data record | Role |
|---|---:|---|
| SingleMuon Run2016G | 30530 | collision data |
| SingleMuon Run2016H | 30563 | collision data |
| TTToSemiLeptonic | 67993 | signal simulation |
| WJetsToLNu | 69747 | main background |
| DYJetsToLL M-50 | 35671 | secondary background |

## Phase 2: muon selection

The muon selection follows the offline requirements of the reference
analysis as closely as the 2016 NanoAOD content permits:

- `HLT_IsoMu24 || HLT_IsoTkMu24` single-muon trigger;
- at least one muon with `Muon_pt > 30 GeV` and `abs(Muon_eta) < 2.1`;
- exactly one such muon satisfying `Muon_tightId`;
- track relative isolation `Muon_tkRelIso < 0.05`;
- veto of any additional loose muon with `pt > 15 GeV` and `abs(eta) < 2.4`.

Run the Phase-2 selection on one file per configured sample with:

```bash
./scripts/run_phase2.sh 100000
```

The argument is the maximum number of input events read from each file. Text
cut flows, ROOT histograms, and PNG diagnostic figures are written under
`output/phase2/`.

## Phase 3: top-like event selection

The top-like selection extends Phase 2 with the detector level requirements
used in the reference article:

- veto additional loose electrons with `pt > 15 GeV` and `abs(eta) < 2.4`;
- select jets with `pt > 30 GeV`, `abs(eta) < 2.4`, `Jet_jetId >= 2`, and
  `DeltaR(jet, muon) > 0.4`;
- require at least four selected jets;
- require at least two DeepJet-loose b jets (`Jet_btagDeepFlavB > 0.0480`);
- require at least one DeepJet-medium b jet (`Jet_btagDeepFlavB > 0.2489`);
- require at least one of the two highest-b-tag jets and one remaining jet to
  have `pt > 35 GeV`.

The DeepJet thresholds are the UL2016 post-VFP working points appropriate for
Run G/H. Run the selection with:

```bash
./scripts/run_phase3.sh 100000
```

Outputs are written under `output/phase3/`.

## Phase 4: Monte Carlo normalization

For every simulated event the luminosity weight is

```text
eventWeight = genWeight * luminosity * crossSection / sumGenWeights
```

The target luminosity and effective cross sections are stored in
`config/normalization.cfg`.  The current target is the combined 2016 G+H luminosity, 16393 pb^-1.

For W+jets, `sumGenWeights` and the normalization factor are calculated
independently for `LHE_Njets == 0`, `1`, and `2`. Events with more than two LHE
partons are reported but receive no normalization because no cross section was
provided for them.

Process every event in the representative Phase-4 files with:

```bash
./scripts/run_phase4.sh
```

For a quick test limited to 10000 events per sample:

```bash
./scripts/run_phase4.sh 10000
```

Normalization factors must be recomputed whenever the processed event subset
or file list changes. Reports, ROOT histograms, and diagnostic figures are
written under `output/phase4/`.

## Phase 5: hadronic top-quark reconstruction

In every top-like event, the two selected jets with the largest DeepJet scores
are treated as the two b candidates. Every pair of the remaining jets is tested
as the hadronic W decay, together with both possible assignments of the
hadronic b candidate. The selected permutation minimizes

```text
chi2 = ((m_jj - 80.4 GeV) / 10 GeV)^2
     + ((m_jjb - 172.5 GeV) / 20 GeV)^2.
```

The 10 and 20 GeV denominators are detector resolution parameters used to rank
permutations, not the physical W and top widths. No chi2 cut is applied yet.

Run the reconstruction with:

```bash
./scripts/run_phase5.sh 100000
```

Mass, transverse-momentum, chi2, and mass-correlation outputs are written under
`output/phase5/`.

## Phase 6: background subtraction and inclusive measurement

First retrieve the complete CERN Open Data file lists and the official CMS
2016 certified-luminosity JSON:

```bash
./scripts/fetch_file_lists.sh
```

The measurement uses

```text
N_signal = N_data - N_W+jets - N_DY
sigma_tt = N_signal / (luminosity * BR_semileptonic * efficiency_tt).
```

Here `BR_semileptonic = 0.439` is entered directly from the W decay branching
fractions. It is not calculated from a reference ttbar cross section. The
macro reports both the semileptonic cross section before this correction and
the inclusive ttbar cross section after it. The theoretical inclusive value
is printed only as a final comparison.

Run a non-physical technical test on at most 10000 events per sample with:

```bash
./scripts/run_phase6.sh 10000
```

Run the physical full-dataset calculation with:

```bash
./scripts/run_phase6.sh
```

For long remote productions, the recommended reliable command is:

```bash
./scripts/run_reliable_measurement.sh
```

It processes one CERN file at a time and saves a small checkpoint under
`output/checkpoints/`. Every failed file is attempted three times. Running the
same command again skips completed files, so a network timeout does not restart
the analysis from the beginning. When all checkpoints exist, ROOT combines
them and writes the same final Phase-6 report.

For the faster exam-project configuration, use all collision data but only the
first 10 files from each Monte Carlo sample:

```bash
./scripts/run_reliable_measurement.sh 10
```

The efficiency and background normalization denominators are recomputed from
the same MC subsets, so no manual rescaling by the number of files is needed.
The tradeoff is a larger Monte Carlo statistical uncertainty.

Partial mode deliberately does not quote a cross section, because a small
collision-data subset cannot be assigned the full G+H luminosity. It also does
not compare the partial data yield with a Monte Carlo prediction normalized to
the full luminosity, since that ratio would have no physical meaning. The current
background model contains W+jets and DY only; single-top and multijet samples
were not supplied, so the eventual result remains a simplified measurement.

## Unfolding: hadronic-top pT unfolding

`macros/unfoldTopPt.C` constructs a five bin response matrix between generated
and reconstructed hadronic-top transverse momentum. It implements four
iterations of Bayesian (D'Agostini) unfolding, subtracts W+jets and DY, and
performs a closure test using separate halves of the TT simulation.

Use all data files and 10 files from each Monte Carlo sample with:

```bash
root -l -b -q 'macros/unfoldTopPt.C(-1,10,4)'
```

For a quick technical test with two data files and two MC files:

```bash
root -l -b -q 'macros/unfoldTopPt.C(2,2,4,5000)'
```

Outputs are saved under `output/unfolding/`. The reported uncertainties are
approximate statistical uncertainties; this remains a simplified
unfolding rather than the complete procedure used by the publication.

For the full-data production, use the reliable file-by-file workflow:

```bash
./scripts/run_reliable_unfolding.sh 10
```

`processUnfoldingFile.C` writes one small ROOT checkpoint per input file under
`output/unfolding_checkpoints/`. `combineUnfolding.C` combines all 152 data
checkpoints and the requested MC subset, then performs background subtraction,
closure, unfolding, and the differential cross-section calculation. Rerunning
the script skips completed files.
