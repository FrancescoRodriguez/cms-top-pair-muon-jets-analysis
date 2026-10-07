#include "LumiMask.h"
#include "topLikeSelection.C"

#include <ROOT/RDataFrame.hxx>
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
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace InclusiveMeasurement {

std::vector<std::string> ReadFileList(const std::string &path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Cannot open file list: " + path);
  std::vector<std::string> files;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#') continue;
    files.push_back(line);
  }
  if (files.empty()) throw std::runtime_error("Empty file list: " + path);
  return files;
}

std::map<std::string, double> ReadConfiguration(const std::string &path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Cannot open configuration: " + path);
  std::map<std::string, double> values;
  std::string line;
  while (std::getline(input, line)) {
    const auto first = line.find_first_not_of(" \t");
    if (first == std::string::npos || line[first] == '#') continue;
    const auto separator = line.find('=', first);
    if (separator == std::string::npos)
      throw std::runtime_error("Invalid configuration line: " + line);
    values[line.substr(first, separator - first)] =
        std::stod(line.substr(separator + 1));
  }
  return values;
}

double Required(const std::map<std::string, double> &values,
                const std::string &key) {
  const auto found = values.find(key);
  if (found == values.end())
    throw std::runtime_error("Missing configuration key: " + key);
  return found->second;
}

ROOT::RDF::RNode Limit(ROOT::RDF::RNode node, Long64_t maxEvents) {
  return maxEvents > 0 ? node.Range(maxEvents) : node;
}

ROOT::RDF::RNode ApplyTopLikeSelection(ROOT::RDF::RNode input) {
  using namespace TopLikeSelection;
  auto triggered = input.Filter("HLT_IsoMu24 || HLT_IsoTkMu24");
  auto muons = triggered
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
      .Define("nLooseMuon", "Sum(LooseMuonMask)")
      .Filter("nKinematicMuon >= 1")
      .Filter("nTightMuon == 1")
      .Filter("nSelectedMuon == 1")
      .Filter("nLooseMuon == 1")
      .Define("SelectedMuon_eta", "Muon_eta[SelectedMuonMask][0]")
      .Define("SelectedMuon_phi", "Muon_phi[SelectedMuonMask][0]");
  auto leptons = muons
      .Define("VetoElectronMask", VetoElectronMask,
              {"Electron_pt", "Electron_eta", "Electron_cutBased"})
      .Define("nVetoElectron", "Sum(VetoElectronMask)")
      .Filter("nVetoElectron == 0");
  return leptons
      .Define("GoodJetMask", GoodJetMask,
              {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_jetId",
               "SelectedMuon_eta", "SelectedMuon_phi"})
      .Define("nGoodJet", "Sum(GoodJetMask)")
      .Define("nLooseBJet",
              "Sum(GoodJetMask && Jet_btagDeepFlavB > 0.0480)")
      .Define("nMediumBJet",
              "Sum(GoodJetMask && Jet_btagDeepFlavB > 0.2489)")
      .Define("PassesJetPt35Topology", PassesJetPt35Topology,
              {"Jet_pt", "Jet_btagDeepFlavB", "GoodJetMask"})
      .Filter("nGoodJet >= 4")
      .Filter("nLooseBJet >= 2")
      .Filter("nMediumBJet >= 1")
      .Filter("PassesJetPt35Topology");
}

struct WeightedYield {
  double totalSumW = 0.0;
  double selectedSumW = 0.0;
  double selectedSumW2 = 0.0;
  double crossSectionPb = 0.0;

  double Efficiency() const {
    return totalSumW != 0.0 ? selectedSumW / totalSumW : 0.0;
  }
  double Expected(double luminosityPb) const {
    return luminosityPb * crossSectionPb * Efficiency();
  }
  double ExpectedStatError(double luminosityPb) const {
    if (totalSumW == 0.0) return 0.0;
    return luminosityPb * crossSectionPb / totalSumW *
           std::sqrt(selectedSumW2);
  }
};

}  // namespace InclusiveMeasurement

int measureInclusive(Long64_t maxEventsPerSample = -1,
                     const char *fileListDirectory = "config/filelists",
                     const char *configurationPath = "config/normalization.cfg") {
  using namespace InclusiveMeasurement;

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(true);

  const bool partialMode = maxEventsPerSample > 0;
  const std::string directory(fileListDirectory);
  std::map<std::string, double> configuration;
  std::vector<std::string> dataFiles;
  std::vector<std::string> ttFiles;
  std::vector<std::string> wFiles;
  std::vector<std::string> dyFiles;
  std::shared_ptr<LumiMask> lumiMask;
  try {
    configuration = ReadConfiguration(configurationPath);
    auto dataG = ReadFileList(directory + "/data_Run2016G.txt");
    auto dataH = ReadFileList(directory + "/data_Run2016H.txt");
    dataFiles = dataG;
    dataFiles.insert(dataFiles.end(), dataH.begin(), dataH.end());
    ttFiles = ReadFileList(directory + "/ttbar_semileptonic.txt");
    wFiles = ReadFileList(directory + "/wjets.txt");
    dyFiles = ReadFileList(directory + "/dyjets.txt");
    lumiMask = std::make_shared<LumiMask>(
        directory + "/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt");
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }

  const double luminosityPb = Required(configuration, "luminosity_pb");
  const double ttSemileptonicXsecPb =
      Required(configuration, "ttbar_semileptonic_xsec_pb");
  const double semileptonicBranchingFraction =
      Required(configuration, "ttbar_semileptonic_br");
  const double ttInclusiveReferencePb =
      Required(configuration, "ttbar_inclusive_reference_pb");
  const double dyXsecPb = Required(configuration, "dyjets_xsec_pb");
  const double wXsecs[3] = {
      Required(configuration, "w0j_xsec_pb"),
      Required(configuration, "w1j_xsec_pb"),
      Required(configuration, "w2j_xsec_pb")};
  if (semileptonicBranchingFraction <= 0.0 ||
      semileptonicBranchingFraction >= 1.0) {
    std::cerr << "ERROR: ttbar_semileptonic_br must be between zero and one.\n";
    return 1;
  }

  ROOT::RDataFrame dataFrame("Events", dataFiles);
  ROOT::RDataFrame ttFrame("Events", ttFiles);
  ROOT::RDataFrame wFrame("Events", wFiles);
  ROOT::RDataFrame dyFrame("Events", dyFiles);
  std::cout << "Input chains created: " << dataFiles.size() << " data, "
            << ttFiles.size() << " TT, " << wFiles.size() << " W+jets, "
            << dyFiles.size() << " DY files.\n" << std::flush;

  ROOT::RDF::RNode dataInput = Limit(dataFrame, maxEventsPerSample);
  ROOT::RDF::RNode ttInput = Limit(ttFrame, maxEventsPerSample);
  ROOT::RDF::RNode wInput = Limit(wFrame, maxEventsPerSample);
  ROOT::RDF::RNode dyInput = Limit(dyFrame, maxEventsPerSample);

  auto certifiedData = dataInput.Filter(
      [lumiMask](unsigned int run, unsigned int luminosityBlock) {
        return lumiMask->Accept(run, luminosityBlock);
      }, {"run", "luminosityBlock"}, "certified luminosity mask");
  auto selectedData = ApplyTopLikeSelection(certifiedData);
  auto ttWeighted = ttInput.Define("analysisWeight", "double(genWeight)");
  auto wWeighted = wInput.Define("analysisWeight", "double(genWeight)");
  auto dyWeighted = dyInput.Define("analysisWeight", "double(genWeight)");
  auto selectedTT = ApplyTopLikeSelection(ttWeighted);
  auto selectedW = ApplyTopLikeSelection(wWeighted);
  auto selectedDY = ApplyTopLikeSelection(dyWeighted);

  auto dataCount = selectedData.Count();

  auto ttTotalSumW = ttWeighted.Sum<double>("analysisWeight");
  auto ttSelectedSumW = selectedTT.Sum<double>("analysisWeight");
  auto ttWithW2 = selectedTT.Define(
      "genWeight2", "analysisWeight * analysisWeight");
  auto ttSelectedSumW2 = ttWithW2.Sum<double>("genWeight2");

  auto dyTotalSumW = dyWeighted.Sum<double>("analysisWeight");
  auto dySelectedSumW = selectedDY.Sum<double>("analysisWeight");
  auto dyWithW2 = selectedDY.Define(
      "genWeight2", "analysisWeight * analysisWeight");
  auto dySelectedSumW2 = dyWithW2.Sum<double>("genWeight2");

  std::vector<ROOT::RDF::RResultPtr<double>> wTotalSumW;
  std::vector<ROOT::RDF::RResultPtr<double>> wSelectedSumW;
  std::vector<ROOT::RDF::RResultPtr<double>> wSelectedSumW2;
  for (int category = 0; category <= 2; ++category) {
    const std::string selection = "LHE_Njets == " + std::to_string(category);
    wTotalSumW.push_back(
        wWeighted.Filter(selection).Sum<double>("analysisWeight"));
    auto selectedCategory = selectedW.Filter(selection)
        .Define("genWeight2_cat" + std::to_string(category),
                "analysisWeight * analysisWeight");
    wSelectedSumW.push_back(
        selectedCategory.Sum<double>("analysisWeight"));
    wSelectedSumW2.push_back(selectedCategory.Sum<double>(
        "genWeight2_cat" + std::to_string(category)));
  }

  std::cout << "Analysis graphs booked; starting event loops.\n" << std::flush;

  WeightedYield tt{*ttTotalSumW, *ttSelectedSumW, *ttSelectedSumW2,
                   ttSemileptonicXsecPb};
  std::cout << "TT event loop completed.\n" << std::flush;
  WeightedYield dy{*dyTotalSumW, *dySelectedSumW, *dySelectedSumW2,
                   dyXsecPb};
  std::cout << "DY event loop completed.\n" << std::flush;
  std::vector<WeightedYield> wCategories;
  for (int category = 0; category <= 2; ++category)
    wCategories.push_back({*wTotalSumW[category], *wSelectedSumW[category],
                           *wSelectedSumW2[category], wXsecs[category]});
  std::cout << "W+jets event loops completed.\n" << std::flush;

  const double dataYield = static_cast<double>(*dataCount);
  std::cout << "Data event loop completed.\n" << std::flush;
  const double dyYield = dy.Expected(luminosityPb);
  double wYield = 0.0;
  double backgroundVariance = std::pow(dy.ExpectedStatError(luminosityPb), 2);
  for (const auto &category : wCategories) {
    wYield += category.Expected(luminosityPb);
    backgroundVariance +=
        std::pow(category.ExpectedStatError(luminosityPb), 2);
  }
  const double backgroundYield = wYield + dyYield;
  const double signalYield = dataYield - backgroundYield;
  const double expectedTTYield = tt.Expected(luminosityPb);
  const double semileptonicDenominator = luminosityPb * tt.Efficiency();
  const double inclusiveDenominator = semileptonicDenominator *
                                      semileptonicBranchingFraction;
  const double measuredSemileptonicCrossSection = semileptonicDenominator > 0.0
      ? signalYield / semileptonicDenominator
      : std::numeric_limits<double>::quiet_NaN();
  const double measuredInclusiveCrossSection = inclusiveDenominator > 0.0
      ? signalYield / inclusiveDenominator
      : std::numeric_limits<double>::quiet_NaN();
  const double dataStatError = inclusiveDenominator > 0.0
      ? std::sqrt(dataYield) / inclusiveDenominator : 0.0;
  const double mcStatError = inclusiveDenominator > 0.0
      ? std::sqrt(backgroundVariance) / inclusiveDenominator : 0.0;
  const double signalToPredictionRatio = expectedTTYield != 0.0
      ? signalYield / expectedTTYield
      : std::numeric_limits<double>::quiet_NaN();

  const std::string outputDirectory = "output/phase6";
  gSystem->mkdir(outputDirectory.c_str(), true);
  const std::string reportPath = outputDirectory + "/inclusive_measurement.txt";
  std::ofstream report(reportPath);
  report << "Simplified inclusive ttbar cross-section measurement\n";
  report << "====================================================\n";
  report << "Mode: " << (partialMode ? "PARTIAL TECHNICAL TEST" : "FULL DATASET")
         << "\n";
  report << "Maximum events per sample: " << maxEventsPerSample << "\n";
  report << "Certified runs in lumi mask: " << lumiMask->NumberOfRuns() << "\n";
  report << "Target luminosity: " << luminosityPb << " pb^-1\n\n";
  report << std::setprecision(10);
  report << "TT semileptonic weighted efficiency: " << tt.Efficiency() << "\n";
  report << "Semileptonic branching fraction used: "
         << semileptonicBranchingFraction << "\n";
  report << "Selected data events: " << dataYield << "\n";
  report << "Expected W+jets background: " << wYield << "\n";
  report << "Expected DY background: " << dyYield << "\n";
  report << "Total modeled background: " << backgroundYield << "\n";
  report << "Background-subtracted signal: " << signalYield << "\n\n";
  if (partialMode) {
    report << "NO PHYSICAL CROSS SECTION IS QUOTED. A limited data event subset\n"
              "does not correspond to the full 16.146 fb^-1 luminosity. Run with\n"
              "maxEventsPerSample = -1 for the physical measurement.\n"
              "For the same reason, no data-signal / full-luminosity MC ratio\n"
              "is quoted in this technical mode.\n";
  } else {
    report << "Expected TT signal from simulation: " << expectedTTYield << "\n";
    report << "(Data-background)/expected TT: " << signalToPredictionRatio
           << "\n\n";
    report << "Measured semileptonic ttbar cross section: "
           << measuredSemileptonicCrossSection << " pb\n";
    report << "Measured inclusive ttbar cross section: "
           << measuredInclusiveCrossSection << " pb\n";
    report << "Inclusive data statistical uncertainty: " << dataStatError
           << " pb\n";
    report << "Inclusive MC-background statistical uncertainty: "
           << mcStatError << " pb\n";
    report << "Reference ttbar cross section: " << ttInclusiveReferencePb << " pb\n";
  }
  report << "\nCAVEAT: only W+jets and DY backgrounds are modeled. Single-top and\n"
            "multijet backgrounds are absent because no corresponding samples\n"
            "were supplied. The result is therefore a simplified measurement.\n";

  auto hYields = std::make_unique<TH1D>(
      "h_yields", "Selected event yields", 4, 0., 4.);
  hYields->SetDirectory(nullptr);
  hYields->GetXaxis()->SetBinLabel(1, "Data");
  hYields->GetXaxis()->SetBinLabel(2, "W+jets");
  hYields->GetXaxis()->SetBinLabel(3, "DY");
  hYields->GetXaxis()->SetBinLabel(4, "Data-background");
  hYields->SetBinContent(1, dataYield);
  hYields->SetBinContent(2, wYield);
  hYields->SetBinContent(3, dyYield);
  hYields->SetBinContent(4, signalYield);
  hYields->SetFillColor(kAzure - 9);
  hYields->SetLineColor(kAzure + 2);
  hYields->SetLineWidth(2);
  hYields->SetStats(false);
  hYields->GetYaxis()->SetTitle("Events");

  TCanvas canvas("canvas", "Inclusive measurement yields", 850, 650);
  canvas.SetLeftMargin(0.14);
  hYields->Draw("HIST");
  const std::string figurePath = outputDirectory + "/inclusive_yields.png";
  canvas.SaveAs(figurePath.c_str());
  const std::string rootPath = outputDirectory + "/inclusive_measurement.root";
  TFile outputFile(rootPath.c_str(), "RECREATE");
  hYields->Write();
  outputFile.Close();

  std::cout << "Selected data events: " << dataYield << "\n";
  if (!partialMode) {
    std::cout << "Measured semileptonic cross section: "
              << measuredSemileptonicCrossSection << " pb\n";
    std::cout << "Measured inclusive cross section: "
              << measuredInclusiveCrossSection << " pb\n";
  }
  std::cout << "Mode: " << (partialMode ? "partial technical test" : "full dataset")
            << "\n";
  std::cout << "Report: " << reportPath << "\n";
  return 0;
}
