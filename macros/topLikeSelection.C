#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <string>
#include <vector>

namespace TopLikeSelection {

using ROOT::VecOps::RVec;

// Lepton requirements from the reference analysis.
constexpr float kMuonPtMin = 30.0f;
constexpr float kMuonAbsEtaMax = 2.1f;
constexpr float kMuonTrackIsoMax = 0.05f;
constexpr float kVetoLeptonPtMin = 15.0f;
constexpr float kVetoLeptonAbsEtaMax = 2.4f;

// Detector-level jet requirements from Sec. VI of the reference analysis.
constexpr float kJetPtMin = 30.0f;
constexpr float kJetAbsEtaMax = 2.4f;
constexpr float kJetMuonDeltaRMin = 0.4f;
constexpr float kAdditionalJetPtMin = 35.0f;


constexpr float kDeepJetLoose = 0.0480f;
constexpr float kDeepJetMedium = 0.2489f;

RVec<bool> KinematicMuonMask(const RVec<float> &pt,
                             const RVec<float> &eta) {
  return (pt > kMuonPtMin) && (abs(eta) < kMuonAbsEtaMax);
}

RVec<bool> TightMuonMask(const RVec<float> &pt, const RVec<float> &eta,
                         const RVec<bool> &tightId) {
  return KinematicMuonMask(pt, eta) && tightId;
}

RVec<bool> SelectedMuonMask(const RVec<float> &pt, const RVec<float> &eta,
                            const RVec<bool> &tightId,
                            const RVec<float> &trackIso) {
  return TightMuonMask(pt, eta, tightId) && (trackIso < kMuonTrackIsoMax);
}

RVec<bool> LooseMuonMask(const RVec<float> &pt, const RVec<float> &eta,
                         const RVec<bool> &looseId) {
  return (pt > kVetoLeptonPtMin) && (abs(eta) < kVetoLeptonAbsEtaMax) && looseId;
}

RVec<bool> VetoElectronMask(const RVec<float> &pt, const RVec<float> &eta,
                            const RVec<int> &cutBased) {
  return (pt > kVetoLeptonPtMin) && (abs(eta) < kVetoLeptonAbsEtaMax) &&
         (cutBased >= 1);
}

float DeltaPhi(float first, float second) {
  float difference = first - second;
  constexpr float pi = 3.14159265358979323846f;
  while (difference > pi) difference -= 2.0f * pi;
  while (difference <= -pi) difference += 2.0f * pi;
  return difference;
}

RVec<bool> GoodJetMask(const RVec<float> &pt, const RVec<float> &eta,
                       const RVec<float> &phi, const RVec<int> &jetId,
                       float muonEta, float muonPhi) {
  RVec<bool> mask(pt.size(), false);
  for (std::size_t index = 0; index < pt.size(); ++index) {
    const float deltaEta = eta[index] - muonEta;
    const float deltaPhi = DeltaPhi(phi[index], muonPhi);
    const float deltaR2 = deltaEta * deltaEta + deltaPhi * deltaPhi;
    mask[index] = pt[index] > kJetPtMin &&
                  std::abs(eta[index]) < kJetAbsEtaMax &&
                  jetId[index] >= 2 &&
                  deltaR2 > kJetMuonDeltaRMin * kJetMuonDeltaRMin;
  }
  return mask;
}

// Reproduces the extra pT condition in Sec. VI: at least one of the two jets
// with highest b discriminant and at least one of the remaining jets must have
// pT > 35 GeV.
bool PassesJetPt35Topology(const RVec<float> &pt, const RVec<float> &btag,
                          const RVec<bool> &goodMask) {
  std::vector<std::size_t> indices;
  for (std::size_t index = 0; index < goodMask.size(); ++index)
    if (goodMask[index]) indices.push_back(index);
  if (indices.size() < 4) return false;

  std::sort(indices.begin(), indices.end(),
            [&btag](std::size_t first, std::size_t second) {
              return btag[first] > btag[second];
            });
  const bool highPtBcandidate =
      pt[indices[0]] > kAdditionalJetPtMin ||
      pt[indices[1]] > kAdditionalJetPtMin;
  bool highPtNonBcandidate = false;
  for (std::size_t position = 2; position < indices.size(); ++position)
    if (pt[indices[position]] > kAdditionalJetPtMin)
      highPtNonBcandidate = true;
  return highPtBcandidate && highPtNonBcandidate;
}

float LeadingSelectedValue(const RVec<float> &values,
                           const RVec<bool> &mask) {
  float leading = -1.0f;
  for (std::size_t index = 0; index < values.size(); ++index)
    if (mask[index]) leading = std::max(leading, values[index]);
  return leading;
}

void StyleHistogram(TH1 *histogram, const char *xTitle) {
  histogram->SetLineColor(kAzure + 2);
  histogram->SetFillColorAlpha(kAzure - 9, 0.60);
  histogram->SetLineWidth(2);
  histogram->SetStats(false);
  histogram->GetXaxis()->SetTitle(xTitle);
  histogram->GetYaxis()->SetTitle("Events");
}

}  // namespace TopLikeSelection

int topLikeSelection(const char *fileUrl,
                     const char *sampleLabel = "sample",
                     Long64_t maxEvents = 100000) {
  using namespace TopLikeSelection;

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(true);

  if (!fileUrl || std::string(fileUrl).empty()) {
    std::cerr << "ERROR: the input file URL is empty.\n";
    return 1;
  }

  const std::string label = sampleLabel ? sampleLabel : "sample";
  const std::string outputDirectory = "output/phase3";
  gSystem->mkdir(outputDirectory.c_str(), true);

  std::cout << "Opening " << fileUrl << "\n";
  ROOT::RDataFrame dataFrame("Events", fileUrl);
  ROOT::RDF::RNode input = dataFrame;
  if (maxEvents > 0) input = input.Range(maxEvents);

  auto triggered = input.Filter("HLT_IsoMu24 || HLT_IsoTkMu24",
                                "single-muon trigger");
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

  auto kinematicMuon = withMuonMasks.Filter("nKinematicMuon >= 1");
  auto tightMuon = kinematicMuon.Filter("nTightMuon == 1");
  auto isolatedMuon = tightMuon.Filter("nSelectedMuon == 1");
  auto singleMuon = isolatedMuon.Filter("nLooseMuon == 1")
      .Define("SelectedMuon_eta", "Muon_eta[SelectedMuonMask][0]")
      .Define("SelectedMuon_phi", "Muon_phi[SelectedMuonMask][0]");

  auto withElectronVeto = singleMuon
      .Define("VetoElectronMask", VetoElectronMask,
              {"Electron_pt", "Electron_eta", "Electron_cutBased"})
      .Define("nVetoElectron", "Sum(VetoElectronMask)");
  auto zeroElectrons = withElectronVeto.Filter("nVetoElectron == 0");

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

  auto fourJets = withJets.Filter("nGoodJet >= 4");
  auto twoLooseBJets = fourJets.Filter("nLooseBJet >= 2");
  auto oneMediumBJet = twoLooseBJets.Filter("nMediumBJet >= 1");
  auto topLike = oneMediumBJet.Filter("PassesJetPt35Topology")
      .Define("GoodJetPt", "Jet_pt[GoodJetMask]")
      .Define("HT", "Sum(GoodJetPt)")
      .Define("LeadingGoodJetPt", LeadingSelectedValue,
              {"Jet_pt", "GoodJetMask"});

  auto countAll = input.Count();
  auto countTrigger = triggered.Count();
  auto countKinematicMuon = kinematicMuon.Count();
  auto countTightMuon = tightMuon.Count();
  auto countIsolatedMuon = isolatedMuon.Count();
  auto countSingleMuon = singleMuon.Count();
  auto countZeroElectrons = zeroElectrons.Count();
  auto countFourJets = fourJets.Count();
  auto countTwoLooseBJets = twoLooseBJets.Count();
  auto countOneMediumBJet = oneMediumBJet.Count();
  auto countTopLike = topLike.Count();

  auto hNJet = withJets.Histo1D(
      {"h_nGoodJet", "Good-jet multiplicity before jet cut", 13, -0.5, 12.5},
      "nGoodJet");
  auto hNBJet = fourJets.Histo1D(
      {"h_nLooseBJet", "Loose b-jet multiplicity after four-jet cut", 7, -0.5, 6.5},
      "nLooseBJet");
  auto hLeadingJetPt = topLike.Histo1D(
      {"h_leadingJetPt", "Leading selected-jet pT", 50, 0., 300.},
      "LeadingGoodJetPt");
  auto hHT = topLike.Histo1D(
      {"h_ht", "Scalar sum of selected-jet pT", 60, 0., 900.}, "HT");

  const std::vector<std::string> cutNames = {
      "All events", "Single-muon trigger", "Muon pT and eta",
      "Exactly one tight muon", "Track isolation", "Extra-muon veto",
      "Electron veto", "At least 4 good jets", "At least 2 loose b jets",
      "At least 1 medium b jet", "35 GeV topology"};
  const std::vector<ULong64_t> cutCounts = {
      *countAll, *countTrigger, *countKinematicMuon, *countTightMuon,
      *countIsolatedMuon, *countSingleMuon, *countZeroElectrons,
      *countFourJets, *countTwoLooseBJets, *countOneMediumBJet, *countTopLike};

  const std::string reportPath = outputDirectory + "/" + label + "_cutflow.txt";
  std::ofstream report(reportPath);
  report << "Top-like event-selection cut flow\n";
  report << "=================================\n";
  report << "Sample: " << label << "\n";
  report << "File: " << fileUrl << "\n";
  report << "Maximum input events: " << maxEvents << "\n\n";
  report << "Jet selection: pT > " << kJetPtMin << " GeV, |eta| < "
         << kJetAbsEtaMax << ", Jet_jetId >= 2, DeltaR(jet,mu) > "
         << kJetMuonDeltaRMin << "\n";
  report << "DeepJet loose WP: " << kDeepJetLoose << "\n";
  report << "DeepJet medium WP: " << kDeepJetMedium << "\n";
  report << "W2J convention for later normalization: LHE_Njets == 2\n\n";
  report << std::left << std::setw(30) << "Cut" << std::right
         << std::setw(14) << "Events" << std::setw(16) << "Step eff."
         << std::setw(16) << "Total eff." << "\n";
  report << std::string(76, '-') << "\n";
  for (std::size_t index = 0; index < cutCounts.size(); ++index) {
    const double stepEfficiency = index == 0
        ? 1.0
        : (cutCounts[index - 1] == 0
            ? 0.0
            : static_cast<double>(cutCounts[index]) / cutCounts[index - 1]);
    const double totalEfficiency = cutCounts[0] == 0
        ? 0.0 : static_cast<double>(cutCounts[index]) / cutCounts[0];
    report << std::left << std::setw(30) << cutNames[index] << std::right
           << std::setw(14) << cutCounts[index]
           << std::setw(15) << std::fixed << std::setprecision(4)
           << stepEfficiency << std::setw(16) << totalEfficiency << "\n";
  }

  auto hCutflow = std::make_unique<TH1D>(
      "h_cutflow", "Top-like event-selection cut flow",
      cutCounts.size(), 0., cutCounts.size());
  hCutflow->SetDirectory(nullptr);
  for (std::size_t index = 0; index < cutCounts.size(); ++index) {
    hCutflow->SetBinContent(index + 1, cutCounts[index]);
    hCutflow->GetXaxis()->SetBinLabel(index + 1, cutNames[index].c_str());
  }
  StyleHistogram(hCutflow.get(), "Selection step");
  hCutflow->LabelsOption("v", "X");
  StyleHistogram(hNJet.GetPtr(), "N_{jets}");
  StyleHistogram(hNBJet.GetPtr(), "N_{b jets} (loose)");
  StyleHistogram(hLeadingJetPt.GetPtr(), "Leading jet p_{T} [GeV]");
  StyleHistogram(hHT.GetPtr(), "H_{T} [GeV]");

  TCanvas canvas("canvas", "Top-like selection", 1300, 900);
  canvas.Divide(2, 2);
  canvas.cd(1); gPad->SetLeftMargin(0.14); hNJet->Draw("HIST");
  canvas.cd(2); gPad->SetLeftMargin(0.14); hNBJet->Draw("HIST");
  canvas.cd(3); gPad->SetLeftMargin(0.14); hHT->Draw("HIST");
  canvas.cd(4);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.40);
  gPad->SetLogy();
  hCutflow->SetMinimum(0.5);
  hCutflow->Draw("HIST");

  const std::string figurePath = outputDirectory + "/" + label + "_toplike.png";
  canvas.SaveAs(figurePath.c_str());

  const std::string rootPath = outputDirectory + "/" + label + "_toplike.root";
  TFile outputFile(rootPath.c_str(), "RECREATE");
  hNJet->Write();
  hNBJet->Write();
  hLeadingJetPt->Write();
  hHT->Write();
  hCutflow->Write();
  outputFile.Close();

  std::cout << "Selected " << cutCounts.back() << " / " << cutCounts.front()
            << " top-like events.\n";
  std::cout << "Cut flow: " << reportPath << "\n";
  std::cout << "Figure: " << figurePath << "\n";
  std::cout << "Histograms: " << rootPath << "\n";
  return 0;
}
