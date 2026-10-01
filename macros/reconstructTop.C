#include "topLikeSelection.C"

#include <Math/Vector4D.h>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace TopReconstruction {

using ROOT::Math::PtEtaPhiMVector;
using ROOT::VecOps::RVec;

constexpr double kWBosonMass = 80.4;
constexpr double kTopQuarkMass = 172.5;

// Resolution parameters used only to rank the possible jet permutations.
// They are detector-scale resolutions, not the physical W/top decay widths.
constexpr double kWMassResolution = 10.0;
constexpr double kTopMassResolution = 20.0;

enum ResultIndex {
  kValid = 0,
  kWMass,
  kTopMass,
  kTopPt,
  kTopEta,
  kChi2,
  kHadronicBIndex,
  kLeptonicBIndex,
  kFirstWJetIndex,
  kSecondWJetIndex
};

RVec<double> ReconstructHadronicTop(const RVec<float> &jetPt,
                                    const RVec<float> &jetEta,
                                    const RVec<float> &jetPhi,
                                    const RVec<float> &jetMass,
                                    const RVec<float> &btag,
                                    const RVec<bool> &goodMask) {
  RVec<double> result(10, -1.0);

  std::vector<std::size_t> goodIndices;
  for (std::size_t index = 0; index < goodMask.size(); ++index)
    if (goodMask[index]) goodIndices.push_back(index);
  if (goodIndices.size() < 4) return result;

  std::sort(goodIndices.begin(), goodIndices.end(),
            [&btag](std::size_t first, std::size_t second) {
              return btag[first] > btag[second];
            });
  const std::vector<std::size_t> bCandidates = {
      goodIndices[0], goodIndices[1]};
  const std::vector<std::size_t> wCandidates(
      goodIndices.begin() + 2, goodIndices.end());
  if (wCandidates.size() < 2) return result;

  double bestChi2 = std::numeric_limits<double>::infinity();
  for (std::size_t hadronicBPosition = 0;
       hadronicBPosition < bCandidates.size(); ++hadronicBPosition) {
    const std::size_t hadronicB = bCandidates[hadronicBPosition];
    const std::size_t leptonicB = bCandidates[1 - hadronicBPosition];
    const PtEtaPhiMVector bVector(jetPt[hadronicB], jetEta[hadronicB],
                                 jetPhi[hadronicB], jetMass[hadronicB]);

    for (std::size_t first = 0; first < wCandidates.size(); ++first) {
      for (std::size_t second = first + 1;
           second < wCandidates.size(); ++second) {
        const std::size_t firstWJet = wCandidates[first];
        const std::size_t secondWJet = wCandidates[second];
        const PtEtaPhiMVector firstVector(
            jetPt[firstWJet], jetEta[firstWJet],
            jetPhi[firstWJet], jetMass[firstWJet]);
        const PtEtaPhiMVector secondVector(
            jetPt[secondWJet], jetEta[secondWJet],
            jetPhi[secondWJet], jetMass[secondWJet]);
        const auto wVector = firstVector + secondVector;
        const auto topVector = wVector + bVector;

        const double wPull =
            (wVector.M() - kWBosonMass) / kWMassResolution;
        const double topPull =
            (topVector.M() - kTopQuarkMass) / kTopMassResolution;
        const double chi2 = wPull * wPull + topPull * topPull;
        if (chi2 >= bestChi2) continue;

        bestChi2 = chi2;
        result[kValid] = 1.0;
        result[kWMass] = wVector.M();
        result[kTopMass] = topVector.M();
        result[kTopPt] = topVector.Pt();
        result[kTopEta] = topVector.Eta();
        result[kChi2] = chi2;
        result[kHadronicBIndex] = hadronicB;
        result[kLeptonicBIndex] = leptonicB;
        result[kFirstWJetIndex] = firstWJet;
        result[kSecondWJetIndex] = secondWJet;
      }
    }
  }
  return result;
}

void StyleHistogram(TH1 *histogram, const char *xTitle) {
  histogram->SetLineColor(kAzure + 2);
  histogram->SetFillColorAlpha(kAzure - 9, 0.60);
  histogram->SetLineWidth(2);
  histogram->SetStats(false);
  histogram->GetXaxis()->SetTitle(xTitle);
  histogram->GetYaxis()->SetTitle("Events");
}

}  // namespace TopReconstruction

int reconstructTop(const char *fileUrl,
                   const char *sampleLabel = "sample",
                   Long64_t maxEvents = 100000) {
  using namespace TopLikeSelection;
  using namespace TopReconstruction;

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(true);

  if (!fileUrl || std::string(fileUrl).empty()) {
    std::cerr << "ERROR: the input file URL is empty.\n";
    return 1;
  }

  const std::string label = sampleLabel ? sampleLabel : "sample";
  const std::string outputDirectory = "output/phase5";
  gSystem->mkdir(outputDirectory.c_str(), true);

  std::cout << "Opening " << fileUrl << "\n";
  ROOT::RDataFrame dataFrame("Events", fileUrl);
  ROOT::RDF::RNode input = dataFrame;
  if (maxEvents > 0) input = input.Range(maxEvents);

  // Reuse exactly the Phase-3 selection definitions before reconstructing the
  // hadronically decaying top quark.
  auto triggered = input.Filter("HLT_IsoMu24 || HLT_IsoTkMu24");
  auto withMuonMasks = triggered
      .Define("KinematicMuonMask", KinematicMuonMask, {"Muon_pt", "Muon_eta"})
      .Define("nKinematicMuon", "Sum(KinematicMuonMask)")
      .Define("TightMuonMask", TightMuonMask,
              {"Muon_pt", "Muon_eta", "Muon_tightId"})
      .Define("nTightMuon", "Sum(TightMuonMask)")
      .Define("SelectedMuonMask", SelectedMuonMask,
              {"Muon_pt", "Muon_eta", "Muon_tightId", "Muon_tkRelIso"})
      .Define("nSelectedMuon", "Sum(SelectedMuonMask)")
      .Define("LooseMuonMask", LooseMuonMask,
              {"Muon_pt", "Muon_eta", "Muon_looseId"})
      .Define("nLooseMuon", "Sum(LooseMuonMask)");
  auto singleMuon = withMuonMasks
      .Filter("nKinematicMuon >= 1")
      .Filter("nTightMuon == 1")
      .Filter("nSelectedMuon == 1")
      .Filter("nLooseMuon == 1")
      .Define("SelectedMuon_eta", "Muon_eta[SelectedMuonMask][0]")
      .Define("SelectedMuon_phi", "Muon_phi[SelectedMuonMask][0]");
  auto zeroElectrons = singleMuon
      .Define("VetoElectronMask", VetoElectronMask,
              {"Electron_pt", "Electron_eta", "Electron_cutBased"})
      .Define("nVetoElectron", "Sum(VetoElectronMask)")
      .Filter("nVetoElectron == 0");
  auto withJets = zeroElectrons
      .Define("GoodJetMask", GoodJetMask,
              {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId",
               "SelectedMuon_eta", "SelectedMuon_phi"})
      .Define("nGoodJet", "Sum(GoodJetMask)")
      .Define("nLooseBJet",
              "Sum(GoodJetMask && Jet_btagDeepFlavB > 0.0480)")
      .Define("nMediumBJet",
              "Sum(GoodJetMask && Jet_btagDeepFlavB > 0.2489)")
      .Define("PassesJetPt35Topology", PassesJetPt35Topology,
              {"Jet_pt", "Jet_btagDeepFlavB", "GoodJetMask"});
  auto topLike = withJets
      .Filter("nGoodJet >= 4")
      .Filter("nLooseBJet >= 2")
      .Filter("nMediumBJet >= 1")
      .Filter("PassesJetPt35Topology");

  auto reconstructed = topLike
      .Define("TopRecoResult", ReconstructHadronicTop,
              {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass",
               "Jet_btagDeepFlavB", "GoodJetMask"})
      .Filter("TopRecoResult[0] > 0.5")
      .Define("RecoWMass", "TopRecoResult[1]")
      .Define("RecoTopMass", "TopRecoResult[2]")
      .Define("RecoTopPt", "TopRecoResult[3]")
      .Define("RecoTopEta", "TopRecoResult[4]")
      .Define("RecoChi2", "TopRecoResult[5]");

  auto countInput = input.Count();
  auto countTopLike = topLike.Count();
  auto countReconstructed = reconstructed.Count();
  auto meanWMass = reconstructed.Mean("RecoWMass");
  auto meanTopMass = reconstructed.Mean("RecoTopMass");

  auto hWMass = reconstructed.Histo1D(
      {"h_recoWMass", "Hadronic W candidate mass", 50, 30., 180.},
      "RecoWMass");
  auto hTopMass = reconstructed.Histo1D(
      {"h_recoTopMass", "Hadronic top candidate mass", 60, 70., 310.},
      "RecoTopMass");
  auto hTopPt = reconstructed.Histo1D(
      {"h_recoTopPt", "Hadronic top candidate transverse momentum",
       60, 0., 600.}, "RecoTopPt");
  auto hChi2 = reconstructed.Histo1D(
      {"h_recoChi2", "Best jet-permutation chi2", 60, 0., 60.},
      "RecoChi2");
  auto hMassCorrelation = reconstructed.Histo2D(
      {"h_massCorrelation", "Reconstructed W and top masses",
       50, 30., 180., 60, 70., 310.},
      "RecoWMass", "RecoTopMass");

  TopReconstruction::StyleHistogram(hWMass.GetPtr(), "m_{jj} [GeV]");
  TopReconstruction::StyleHistogram(hTopMass.GetPtr(), "m_{jjb} [GeV]");
  TopReconstruction::StyleHistogram(hTopPt.GetPtr(), "p_{T}(t_{had}) [GeV]");
  TopReconstruction::StyleHistogram(hChi2.GetPtr(), "#chi^{2}_{min}");
  hMassCorrelation->SetStats(false);
  hMassCorrelation->GetXaxis()->SetTitle("m_{jj} [GeV]");
  hMassCorrelation->GetYaxis()->SetTitle("m_{jjb} [GeV]");

  const ULong64_t inputEvents = *countInput;
  const ULong64_t topLikeEvents = *countTopLike;
  const ULong64_t reconstructedEvents = *countReconstructed;
  const std::string reportPath = outputDirectory + "/" + label + "_reconstruction.txt";
  std::ofstream report(reportPath);
  report << "Hadronic top-quark reconstruction\n";
  report << "=================================\n";
  report << "Sample: " << label << "\n";
  report << "File: " << fileUrl << "\n";
  report << "Processed input events: " << inputEvents << "\n";
  report << "Top-like events: " << topLikeEvents << "\n";
  report << "Successfully reconstructed events: " << reconstructedEvents << "\n";
  report << "Reconstruction fraction among top-like events: "
         << (topLikeEvents > 0
                 ? static_cast<double>(reconstructedEvents) / topLikeEvents
                 : 0.0) << "\n\n";
  report << "Mass hypotheses\n";
  report << "  mW = " << kWBosonMass << " GeV\n";
  report << "  mt = " << kTopQuarkMass << " GeV\n";
  report << "Ranking resolutions\n";
  report << "  sigmaW = " << kWMassResolution << " GeV\n";
  report << "  sigmaTop = " << kTopMassResolution << " GeV\n\n";
  report << "Mean reconstructed W mass: " << *meanWMass << " GeV\n";
  report << "Mean reconstructed top mass: " << *meanTopMass << " GeV\n";
  report << "W-mass histogram peak bin: "
         << hWMass->GetBinCenter(hWMass->GetMaximumBin()) << " GeV\n";
  report << "Top-mass histogram peak bin: "
         << hTopMass->GetBinCenter(hTopMass->GetMaximumBin()) << " GeV\n";
  report << "No chi2 selection is applied in this phase.\n";

  TCanvas canvas("canvas", "Hadronic top reconstruction", 1300, 900);
  canvas.Divide(2, 2);
  canvas.cd(1); gPad->SetLeftMargin(0.14); hWMass->Draw("HIST");
  canvas.cd(2); gPad->SetLeftMargin(0.14); hTopMass->Draw("HIST");
  canvas.cd(3); gPad->SetLeftMargin(0.14); hTopPt->Draw("HIST");
  canvas.cd(4); gPad->SetLeftMargin(0.14); hChi2->Draw("HIST");
  const std::string figurePath = outputDirectory + "/" + label + "_reconstruction.png";
  canvas.SaveAs(figurePath.c_str());

  TCanvas correlationCanvas("correlationCanvas", "Mass correlation", 700, 600);
  correlationCanvas.SetLeftMargin(0.14);
  correlationCanvas.SetRightMargin(0.15);
  hMassCorrelation->Draw("COLZ");
  const std::string correlationPath =
      outputDirectory + "/" + label + "_mass_correlation.png";
  correlationCanvas.SaveAs(correlationPath.c_str());

  const std::string rootPath = outputDirectory + "/" + label + "_reconstruction.root";
  TFile outputFile(rootPath.c_str(), "RECREATE");
  hWMass->Write();
  hTopMass->Write();
  hTopPt->Write();
  hChi2->Write();
  hMassCorrelation->Write();
  outputFile.Close();

  std::cout << "Reconstructed " << reconstructedEvents << " / "
            << topLikeEvents << " top-like events.\n";
  std::cout << "Report: " << reportPath << "\n";
  std::cout << "Figure: " << figurePath << "\n";
  return 0;
}
