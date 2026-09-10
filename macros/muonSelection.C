#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace MuonSelection {

using ROOT::VecOps::RVec;

// Offline requirements inspired by Secs. V and VI of Phys. Rev. D 95,
// 092001. Muon_tkRelIso is the NanoAOD quantity closest to the track-only
// relative isolation used in the paper.
constexpr float kMuonPtMin = 30.0f;
constexpr float kMuonAbsEtaMax = 2.1f;
constexpr float kMuonTrackIsoMax = 0.05f;
constexpr float kVetoMuonPtMin = 15.0f;
constexpr float kVetoMuonAbsEtaMax = 2.4f;

RVec<bool> KinematicMask(const RVec<float> &pt, const RVec<float> &eta) {
  return (pt > kMuonPtMin) && (abs(eta) < kMuonAbsEtaMax);
}

RVec<bool> TightMask(const RVec<float> &pt, const RVec<float> &eta,
                     const RVec<bool> &tightId) {
  return KinematicMask(pt, eta) && tightId;
}

RVec<bool> SelectedMask(const RVec<float> &pt, const RVec<float> &eta,
                        const RVec<bool> &tightId,
                        const RVec<float> &trackIso) {
  return TightMask(pt, eta, tightId) && (trackIso < kMuonTrackIsoMax);
}

RVec<bool> LooseVetoMask(const RVec<float> &pt, const RVec<float> &eta,
                         const RVec<bool> &looseId) {
  return (pt > kVetoMuonPtMin) && (abs(eta) < kVetoMuonAbsEtaMax) && looseId;
}

void StyleHistogram(TH1 *histogram, int color, const char *xTitle) {
  histogram->SetLineColor(color);
  histogram->SetLineWidth(2);
  histogram->SetStats(false);
  histogram->GetXaxis()->SetTitle(xTitle);
  histogram->GetYaxis()->SetTitle("Entries");
}

}  // namespace MuonSelection

int muonSelection(const char *fileUrl,
                  const char *sampleLabel = "sample",
                  Long64_t maxEvents = 100000) {
  using namespace MuonSelection;

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(true);

  if (!fileUrl || std::string(fileUrl).empty()) {
    std::cerr << "ERROR: the input file URL is empty.\n";
    return 1;
  }

  const std::string label = sampleLabel ? sampleLabel : "sample";
  const std::string outputDirectory = "output/phase2";
  gSystem->mkdir(outputDirectory.c_str(), true);

  std::cout << "Opening " << fileUrl << "\n";
  ROOT::RDataFrame dataFrame("Events", fileUrl);
  ROOT::RDF::RNode input = dataFrame;
  if (maxEvents > 0) input = input.Range(maxEvents);

  // Trigger decision stored for every event by the CMS High-Level Trigger.
  auto triggered = input.Filter("HLT_IsoMu24 || HLT_IsoTkMu24",
                                "single-muon trigger");

  auto withMasks = triggered
      .Define("KinematicMuonMask", KinematicMask, {"Muon_pt", "Muon_eta"})
      .Define("nKinematicMuon", "Sum(KinematicMuonMask)")
      .Define("TightMuonMask", TightMask,
              {"Muon_pt", "Muon_eta", "Muon_tightId"})
      .Define("nTightMuon", "Sum(TightMuonMask)")
      .Define("SelectedMuonMask", SelectedMask,
              {"Muon_pt", "Muon_eta", "Muon_tightId", "Muon_tkRelIso"})
      .Define("nSelectedMuon", "Sum(SelectedMuonMask)")
      .Define("LooseVetoMuonMask", LooseVetoMask,
              {"Muon_pt", "Muon_eta", "Muon_looseId"})
      .Define("nLooseVetoMuon", "Sum(LooseVetoMuonMask)");

  auto kinematic = withMasks.Filter("nKinematicMuon >= 1",
                                    "at least one kinematic muon");
  auto tight = kinematic.Filter("nTightMuon == 1", "exactly one tight muon");
  auto isolated = tight.Filter("nSelectedMuon == 1",
                               "track isolation below 5 percent");
  auto selected = isolated.Filter("nLooseVetoMuon == 1",
                                  "additional loose-muon veto")
      .Define("SelectedMuon_pt", "Muon_pt[SelectedMuonMask][0]")
      .Define("SelectedMuon_eta", "Muon_eta[SelectedMuonMask][0]")
      .Define("SelectedMuon_tkRelIso", "Muon_tkRelIso[SelectedMuonMask][0]");

  // Book every action before reading a result so RDataFrame can execute them
  // together in a single event loop.
  auto countAll = input.Count();
  auto countTriggered = triggered.Count();
  auto countKinematic = kinematic.Count();
  auto countTight = tight.Count();
  auto countIsolated = isolated.Count();
  auto countSelected = selected.Count();

  auto hPtBefore = triggered.Histo1D(
      {"h_muonPt_before", "Muon pT after trigger", 60, 0., 180.}, "Muon_pt");
  auto hEtaBefore = triggered.Histo1D(
      {"h_muonEta_before", "Muon eta after trigger", 50, -2.5, 2.5}, "Muon_eta");
  auto hIsoBefore = triggered.Histo1D(
      {"h_muonIso_before", "Muon track isolation after trigger", 50, 0., 0.5},
      "Muon_tkRelIso");
  auto hPtSelected = selected.Histo1D(
      {"h_muonPt_selected", "Selected muon pT", 60, 0., 180.},
      "SelectedMuon_pt");
  auto hEtaSelected = selected.Histo1D(
      {"h_muonEta_selected", "Selected muon eta", 50, -2.5, 2.5},
      "SelectedMuon_eta");
  auto hIsoSelected = selected.Histo1D(
      {"h_muonIso_selected", "Selected muon track isolation", 50, 0., 0.5},
      "SelectedMuon_tkRelIso");

  const std::vector<std::string> cutNames = {
      "All events", "Single-muon trigger", "Muon pT and eta",
      "Exactly one tight muon", "Track isolation", "Extra-muon veto"};
  const std::vector<ULong64_t> cutCounts = {
      *countAll, *countTriggered, *countKinematic,
      *countTight, *countIsolated, *countSelected};

  const std::string reportPath = outputDirectory + "/" + label + "_cutflow.txt";
  std::ofstream report(reportPath);
  report << "Muon-selection cut flow\n";
  report << "=======================\n";
  report << "Sample: " << label << "\n";
  report << "File: " << fileUrl << "\n";
  report << "Maximum input events: " << maxEvents << "\n\n";
  report << "Selection constants\n";
  report << "  trigger: HLT_IsoMu24 OR HLT_IsoTkMu24\n";
  report << "  selected muon pT > " << kMuonPtMin << " GeV\n";
  report << "  selected muon |eta| < " << kMuonAbsEtaMax << "\n";
  report << "  Muon_tightId = true\n";
  report << "  Muon_tkRelIso < " << kMuonTrackIsoMax << "\n";
  report << "  veto muon pT > " << kVetoMuonPtMin << " GeV, |eta| < "
         << kVetoMuonAbsEtaMax << ", Muon_looseId = true\n\n";
  report << std::left << std::setw(30) << "Cut" << std::right
         << std::setw(14) << "Events" << std::setw(16) << "Step eff."
         << std::setw(16) << "Total eff." << "\n";
  report << std::string(76, '-') << "\n";
  for (std::size_t index = 0; index < cutCounts.size(); ++index) {
    const double stepEfficiency = index == 0 || cutCounts[index - 1] == 0
        ? 1.0
        : static_cast<double>(cutCounts[index]) / cutCounts[index - 1];
    const double totalEfficiency = cutCounts[0] == 0
        ? 0.0
        : static_cast<double>(cutCounts[index]) / cutCounts[0];
    report << std::left << std::setw(30) << cutNames[index] << std::right
           << std::setw(14) << cutCounts[index]
           << std::setw(15) << std::fixed << std::setprecision(4)
           << stepEfficiency << std::setw(16) << totalEfficiency << "\n";
  }

  auto hCutflow = std::make_unique<TH1D>(
      "h_cutflow", "Muon-selection cut flow", cutCounts.size(), 0., cutCounts.size());
  hCutflow->SetDirectory(nullptr);
  for (std::size_t index = 0; index < cutCounts.size(); ++index) {
    hCutflow->SetBinContent(index + 1, cutCounts[index]);
    hCutflow->GetXaxis()->SetBinLabel(index + 1, cutNames[index].c_str());
  }
  hCutflow->SetFillColor(kAzure - 9);
  hCutflow->SetLineColor(kAzure + 2);
  hCutflow->SetLineWidth(2);
  hCutflow->SetStats(false);
  hCutflow->GetYaxis()->SetTitle("Events");
  hCutflow->LabelsOption("v", "X");

  StyleHistogram(hPtBefore.GetPtr(), kGray + 2, "Muon p_{T} [GeV]");
  StyleHistogram(hEtaBefore.GetPtr(), kGray + 2, "Muon #eta");
  StyleHistogram(hIsoBefore.GetPtr(), kGray + 2, "Muon track relative isolation");
  StyleHistogram(hPtSelected.GetPtr(), kAzure + 2, "Muon p_{T} [GeV]");
  StyleHistogram(hEtaSelected.GetPtr(), kAzure + 2, "Muon #eta");
  StyleHistogram(hIsoSelected.GetPtr(), kAzure + 2, "Muon track relative isolation");

  TCanvas canvas("canvas", "Muon selection", 1200, 900);
  canvas.Divide(2, 2);
  canvas.cd(1);
  gPad->SetLeftMargin(0.14);
  hPtBefore->Draw("HIST");
  hPtSelected->Draw("HIST SAME");
  TLegend legend(0.57, 0.72, 0.88, 0.88);
  legend.AddEntry(hPtBefore.GetPtr(), "After trigger", "l");
  legend.AddEntry(hPtSelected.GetPtr(), "Final muon selection", "l");
  legend.Draw();
  canvas.cd(2);
  gPad->SetLeftMargin(0.14);
  hEtaBefore->Draw("HIST");
  hEtaSelected->Draw("HIST SAME");
  canvas.cd(3);
  gPad->SetLeftMargin(0.14);
  hIsoBefore->Draw("HIST");
  hIsoSelected->Draw("HIST SAME");
  canvas.cd(4);
  gPad->SetLeftMargin(0.14);
  gPad->SetBottomMargin(0.30);
  gPad->SetLogy();
  hCutflow->SetMinimum(0.5);
  hCutflow->Draw("HIST");

  const std::string figurePath = outputDirectory + "/" + label + "_muons.png";
  canvas.SaveAs(figurePath.c_str());

  const std::string rootPath = outputDirectory + "/" + label + "_muons.root";
  TFile outputFile(rootPath.c_str(), "RECREATE");
  hPtBefore->Write();
  hEtaBefore->Write();
  hIsoBefore->Write();
  hPtSelected->Write();
  hEtaSelected->Write();
  hIsoSelected->Write();
  hCutflow->Write();
  outputFile.Close();

  std::cout << "Selected " << cutCounts.back() << " / " << cutCounts.front()
            << " events.\n";
  std::cout << "Cut flow: " << reportPath << "\n";
  std::cout << "Figure: " << figurePath << "\n";
  std::cout << "Histograms: " << rootPath << "\n";
  return 0;
}

