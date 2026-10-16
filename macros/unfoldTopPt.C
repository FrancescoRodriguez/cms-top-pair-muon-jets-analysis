#include "reconstructTop.C"
#include "LumiMask.h"

#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// Simplified one-dimensional unfolding of the hadronic-top pT.  The response
// matrix is built from TTToSemiLeptonic simulation and four iterations of
// D'Agostini's Bayesian method are implemented without external libraries.
namespace TopPtUnfolding {

using ROOT::RDF::RNode;
using ROOT::VecOps::RVec;
constexpr int kNBins = 5;
constexpr double kEdges[kNBins + 1] = {0., 100., 200., 300., 500., 800.};

std::vector<std::string> ReadList(const std::string &path, int maximum = -1) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Cannot open file list: " + path);
  std::vector<std::string> files;
  std::string line;
  while (std::getline(input, line)) {
    if (line.empty() || line.front() == '#') continue;
    files.push_back(line);
    if (maximum > 0 && static_cast<int>(files.size()) >= maximum) break;
  }
  if (files.empty()) throw std::runtime_error("Empty file list: " + path);
  return files;
}

std::map<std::string, double> ReadConfig(const std::string &path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Cannot open configuration: " + path);
  std::map<std::string, double> values;
  std::string line;
  while (std::getline(input, line)) {
    const auto first = line.find_first_not_of(" \t");
    if (first == std::string::npos || line[first] == '#') continue;
    const auto separator = line.find('=', first);
    if (separator == std::string::npos) continue;
    std::string key = line.substr(first, separator - first);
    while (!key.empty() && (key.back() == ' ' || key.back() == '\t')) key.pop_back();
    values[key] = std::stod(line.substr(separator + 1));
  }
  return values;
}

double Get(const std::map<std::string, double> &values,
           const std::string &key) {
  const auto item = values.find(key);
  if (item == values.end()) throw std::runtime_error("Missing key: " + key);
  return item->second;
}

double FoldPt(double pt) { return pt < 0. ? pt : std::min(pt, 799.999); }

bool IsLeptonFromTopW(std::size_t particle, int topIndex,
                      const RVec<int> &pdgId, const RVec<int> &mother) {
  bool crossedW = false;
  int current = static_cast<int>(particle);
  for (int depth = 0; depth < 30; ++depth) {
    if (current < 0 || current >= static_cast<int>(mother.size())) return false;
    current = mother[current];
    if (current < 0 || current >= static_cast<int>(pdgId.size())) return false;
    if (std::abs(pdgId[current]) == 24) crossedW = true;
    if (current == topIndex) return crossedW;
  }
  return false;
}

// Find the last copy top whose W has no e, mu or tau descendant.
double GeneratedHadronicTopPt(const RVec<float> &pt,
                              const RVec<int> &pdgId,
                              const RVec<int> &mother,
                              const RVec<int> &statusFlags) {
  constexpr int lastCopyBit = 1 << 13;
  for (std::size_t top = 0; top < pdgId.size(); ++top) {
    if (std::abs(pdgId[top]) != 6 || !(statusFlags[top] & lastCopyBit)) continue;
    bool leptonic = false;
    for (std::size_t particle = 0; particle < pdgId.size(); ++particle) {
      const int id = std::abs(pdgId[particle]);
      if ((id == 11 || id == 13 || id == 15) &&
          IsLeptonFromTopW(particle, top, pdgId, mother)) {
        leptonic = true;
        break;
      }
    }
    if (!leptonic) return FoldPt(pt[top]);
  }
  return -1.;
}

RNode SelectTopLike(RNode input) {
  using namespace TopLikeSelection;
  auto muons = input.Filter("HLT_IsoMu24 || HLT_IsoTkMu24")
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
      .Filter("nKinematicMuon >= 1 && nTightMuon == 1 && nSelectedMuon == 1 && nLooseMuon == 1")
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
      .Define("nLooseBJet", "Sum(GoodJetMask && Jet_btagDeepFlavB > 0.0480)")
      .Define("nMediumBJet", "Sum(GoodJetMask && Jet_btagDeepFlavB > 0.2489)")
      .Define("PassesJetPt35Topology", PassesJetPt35Topology,
              {"Jet_pt", "Jet_btagDeepFlavB", "GoodJetMask"})
      .Filter("nGoodJet >= 4 && nLooseBJet >= 2 && nMediumBJet >= 1")
      .Filter("PassesJetPt35Topology");
}

RNode AddRecoTop(RNode selected) {
  return selected
      .Define("TopRecoResult", TopReconstruction::ReconstructHadronicTop,
              {"Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass",
               "Jet_btagDeepFlavB", "GoodJetMask"})
      .Filter("TopRecoResult[0] > 0.5")
      .Define("RecoTopPt", [](const RVec<double> &result) {
        return FoldPt(result[TopReconstruction::kTopPt]);
      }, {"TopRecoResult"});
}

RNode AddGenTop(RNode input) {
  return input
      .Define("GenTopPt", GeneratedHadronicTopPt,
              {"GenPart_pt", "GenPart_pdgId", "GenPart_genPartIdxMother",
               "GenPart_statusFlags"})
      .Filter("GenTopPt >= 0.");
}

struct UnfoldedResult {
  std::unique_ptr<TH1D> histogram;
  std::vector<std::vector<double>> kernel;  // [truth][reco]
};

UnfoldedResult BayesianUnfold(const TH1D &measured, const TH2D &migration,
                              const TH1D &truth, int iterations,
                              const std::string &name) {
  const int n = truth.GetNbinsX();
  std::vector<std::vector<double>> response(n, std::vector<double>(n, 0.));
  std::vector<std::vector<double>> kernel(n, std::vector<double>(n, 0.));
  std::vector<double> efficiency(n, 0.), prior(n, 0.), result(n, 0.);
  double priorSum = 0.;
  for (int j = 0; j < n; ++j) {
    const double denominator = truth.GetBinContent(j + 1);
    prior[j] = std::max(0., denominator);
    priorSum += prior[j];
    if (denominator == 0.) continue;
    for (int i = 0; i < n; ++i) {
      // x axis = generated bin j, y axis = reconstructed bin i.
      response[i][j] = migration.GetBinContent(j + 1, i + 1) / denominator;
      efficiency[j] += response[i][j];
    }
  }
  if (priorSum <= 0.) throw std::runtime_error("Empty truth prior");
  for (double &value : prior) value /= priorSum;

  for (int iteration = 0; iteration < iterations; ++iteration) {
    for (auto &row : kernel) std::fill(row.begin(), row.end(), 0.);
    for (int i = 0; i < n; ++i) {
      double normalization = 0.;
      for (int j = 0; j < n; ++j) normalization += response[i][j] * prior[j];
      if (normalization <= 0.) continue;
      for (int j = 0; j < n; ++j)
        if (efficiency[j] > 0.)
          kernel[j][i] = response[i][j] * prior[j] /
                         normalization / efficiency[j];
    }
    std::fill(result.begin(), result.end(), 0.);
    for (int j = 0; j < n; ++j)
      for (int i = 0; i < n; ++i)
        result[j] += kernel[j][i] * measured.GetBinContent(i + 1);
    double sum = 0.;
    for (double value : result) sum += std::max(0., value);
    if (sum > 0.)
      for (int j = 0; j < n; ++j) prior[j] = std::max(0., result[j]) / sum;
  }

  auto output = std::make_unique<TH1D>(name.c_str(), "", n, kEdges);
  output->SetDirectory(nullptr);
  for (int j = 0; j < n; ++j) output->SetBinContent(j + 1, result[j]);
  return {std::move(output), std::move(kernel)};
}

std::unique_ptr<TH2D> NormalizeColumns(const TH2D &input) {
  auto output = std::unique_ptr<TH2D>(
      static_cast<TH2D *>(input.Clone("h_response_normalized")));
  output->SetDirectory(nullptr);
  // Normalize each generated-pT column (x) over all reconstructed bins (y).
  for (int x = 1; x <= output->GetNbinsX(); ++x) {
    double sum = 0.;
    for (int y = 1; y <= output->GetNbinsY(); ++y)
      sum += output->GetBinContent(x, y);
    if (sum == 0.) continue;
    for (int y = 1; y <= output->GetNbinsY(); ++y)
      output->SetBinContent(x, y, output->GetBinContent(x, y) / sum);
  }
  return output;
}

}  // namespace TopPtUnfolding

int unfoldTopPt(int maximumDataFiles = -1, int maximumMCFiles = 10,
                int iterations = 4, Long64_t maximumEventsPerSample = -1) {
  using namespace TopPtUnfolding;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(true);
  if (maximumMCFiles < 2 || iterations < 1) return 1;

  std::vector<std::string> dataFiles, ttFiles, wFiles, dyFiles;
  std::map<std::string, double> config;
  try {
    // A positive data limit is divided approximately equally between G and H.
    const int perEra = maximumDataFiles > 0 ? std::max(1, maximumDataFiles / 2) : -1;
    auto dataG = ReadList("config/filelists/data_Run2016G.txt", perEra);
    auto dataH = ReadList("config/filelists/data_Run2016H.txt", perEra);
    dataFiles = dataG;
    dataFiles.insert(dataFiles.end(), dataH.begin(), dataH.end());
    ttFiles = ReadList("config/filelists/ttbar_semileptonic.txt", maximumMCFiles);
    wFiles = ReadList("config/filelists/wjets.txt", maximumMCFiles);
    dyFiles = ReadList("config/filelists/dyjets.txt", maximumMCFiles);
    config = ReadConfig("config/normalization.cfg");
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 2;
  }

  const std::size_t split = ttFiles.size() / 2;
  std::vector<std::string> trainFiles(ttFiles.begin(), ttFiles.begin() + split);
  std::vector<std::string> testFiles(ttFiles.begin() + split, ttFiles.end());
  ROOT::RDataFrame dataFrame("Events", dataFiles);
  ROOT::RDataFrame ttFrame("Events", ttFiles);
  ROOT::RDataFrame trainFrame("Events", trainFiles);
  ROOT::RDataFrame testFrame("Events", testFiles);
  ROOT::RDataFrame wFrame("Events", wFiles);
  ROOT::RDataFrame dyFrame("Events", dyFiles);

  RNode dataInput = dataFrame;
  RNode ttInput = ttFrame;
  RNode trainInput = trainFrame;
  RNode testInput = testFrame;
  RNode wInput = wFrame;
  RNode dyInput = dyFrame;
  if (maximumEventsPerSample > 0) {
    dataInput = dataInput.Range(maximumEventsPerSample);
    ttInput = ttInput.Range(maximumEventsPerSample);
    trainInput = trainInput.Range(maximumEventsPerSample);
    testInput = testInput.Range(maximumEventsPerSample);
    wInput = wInput.Range(maximumEventsPerSample);
    dyInput = dyInput.Range(maximumEventsPerSample);
  }

  auto mask = std::make_shared<LumiMask>(
      "config/filelists/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt");
  auto dataNode = dataInput.Filter(
      [mask](unsigned int run, unsigned int luminosityBlock) {
        return mask->Accept(run, luminosityBlock);
      }, {"run", "luminosityBlock"});
  auto dataReco = AddRecoTop(SelectTopLike(dataNode));

  auto tt = ttInput.Define("analysisWeight", "double(genWeight)");
  auto ttTruth = AddGenTop(tt);
  auto ttResponse = AddGenTop(AddRecoTop(SelectTopLike(tt)));
  auto train = trainInput.Define("analysisWeight", "double(genWeight)");
  auto trainTruth = AddGenTop(train);
  auto trainResponse = AddGenTop(AddRecoTop(SelectTopLike(train)));
  auto test = testInput.Define("analysisWeight", "double(genWeight)");
  auto testTruth = AddGenTop(test);
  auto testReco = AddGenTop(AddRecoTop(SelectTopLike(test)));
  auto w = wInput.Define("analysisWeight", "double(genWeight)");
  auto wReco = AddRecoTop(SelectTopLike(w));
  auto dy = dyInput.Define("analysisWeight", "double(genWeight)");
  auto dyReco = AddRecoTop(SelectTopLike(dy));

  ROOT::RDF::TH1DModel ptModel("pt", "", kNBins, kEdges);
  ROOT::RDF::TH2DModel matrixModel("matrix", "", kNBins, kEdges, kNBins, kEdges);
  auto hData = dataReco.Histo1D(ptModel, "RecoTopPt");
  auto hTruth = ttTruth.Histo1D(ptModel, "GenTopPt", "analysisWeight");
  auto hMatrix = ttResponse.Histo2D(matrixModel, "GenTopPt", "RecoTopPt", "analysisWeight");
  auto hTrainTruth = trainTruth.Histo1D(ptModel, "GenTopPt", "analysisWeight");
  auto hTrainMatrix = trainResponse.Histo2D(matrixModel, "GenTopPt", "RecoTopPt", "analysisWeight");
  auto hTestTruth = testTruth.Histo1D(ptModel, "GenTopPt", "analysisWeight");
  auto hTestReco = testReco.Histo1D(ptModel, "RecoTopPt", "analysisWeight");
  auto dyTotal = dy.Sum<double>("analysisWeight");
  auto hDY = dyReco.Histo1D(ptModel, "RecoTopPt", "analysisWeight");

  std::vector<ROOT::RDF::RResultPtr<double>> wTotal;
  std::vector<ROOT::RDF::RResultPtr<TH1D>> hW;
  for (int category = 0; category < 3; ++category) {
    const std::string cut = "LHE_Njets == " + std::to_string(category);
    wTotal.push_back(w.Filter(cut).Sum<double>("analysisWeight"));
    hW.push_back(wReco.Filter(cut).Histo1D(ptModel, "RecoTopPt", "analysisWeight"));
  }

  std::cout << "Reading " << dataFiles.size() << " data and " << ttFiles.size()
            << "/" << wFiles.size() << "/" << dyFiles.size()
            << " TT/W/DY files.\n" << std::flush;
  auto data = std::unique_ptr<TH1D>(static_cast<TH1D *>(hData->Clone("h_data_reco")));
  auto truth = std::unique_ptr<TH1D>(static_cast<TH1D *>(hTruth->Clone("h_tt_truth")));
  auto matrix = std::unique_ptr<TH2D>(static_cast<TH2D *>(hMatrix->Clone("h_response_counts")));
  auto trainingTruth = std::unique_ptr<TH1D>(static_cast<TH1D *>(hTrainTruth->Clone("h_training_truth")));
  auto trainingMatrix = std::unique_ptr<TH2D>(static_cast<TH2D *>(hTrainMatrix->Clone("h_training_response")));
  auto closureTruth = std::unique_ptr<TH1D>(static_cast<TH1D *>(hTestTruth->Clone("h_closure_truth")));
  auto closureReco = std::unique_ptr<TH1D>(static_cast<TH1D *>(hTestReco->Clone("h_closure_reco")));
  auto dyHistogram = std::unique_ptr<TH1D>(static_cast<TH1D *>(hDY->Clone("h_dy_reco")));
  for (auto *hist : {data.get(), truth.get(), trainingTruth.get(), closureTruth.get(),
                     closureReco.get(), dyHistogram.get()}) hist->SetDirectory(nullptr);
  matrix->SetDirectory(nullptr);
  trainingMatrix->SetDirectory(nullptr);

  const double luminosity = Get(config, "luminosity_pb");
  const double branching = Get(config, "ttbar_semileptonic_br");
  dyHistogram->Scale(luminosity * Get(config, "dyjets_xsec_pb") / *dyTotal);
  auto wHistogram = std::unique_ptr<TH1D>(static_cast<TH1D *>(data->Clone("h_w_reco")));
  wHistogram->Reset("ICES");
  const double wXsec[3] = {Get(config, "w0j_xsec_pb"), Get(config, "w1j_xsec_pb"),
                           Get(config, "w2j_xsec_pb")};
  for (int category = 0; category < 3; ++category) {
    auto component = std::unique_ptr<TH1D>(static_cast<TH1D *>(hW[category]->Clone()));
    component->SetDirectory(nullptr);
    if (*wTotal[category] != 0.)
      component->Scale(luminosity * wXsec[category] / *wTotal[category]);
    wHistogram->Add(component.get());
  }
  auto background = std::unique_ptr<TH1D>(static_cast<TH1D *>(wHistogram->Clone("h_background")));
  background->Add(dyHistogram.get());
  auto signalReco = std::unique_ptr<TH1D>(static_cast<TH1D *>(data->Clone("h_signal_reco")));
  signalReco->Add(background.get(), -1.);
  int negativeBins = 0;
  for (int bin = 1; bin <= kNBins; ++bin)
    if (signalReco->GetBinContent(bin) < 0.) {
      signalReco->SetBinContent(bin, 0.);
      ++negativeBins;
    }

  UnfoldedResult closure, unfolded;
  try {
    closure = BayesianUnfold(*closureReco, *trainingMatrix, *trainingTruth,
                             iterations, "h_closure_unfolded");
    unfolded = BayesianUnfold(*signalReco, *matrix, *truth,
                              iterations, "h_data_unfolded_events");
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 3;
  }

  // Approximate statistical errors with the final prior held fixed.
  for (int j = 0; j < kNBins; ++j) {
    double variance = 0.;
    for (int i = 0; i < kNBins; ++i) {
      const double term = unfolded.kernel[j][i] * signalReco->GetBinError(i + 1);
      variance += term * term;
    }
    unfolded.histogram->SetBinError(j + 1, std::sqrt(variance));
  }
  auto differential = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(unfolded.histogram->Clone("h_differential_xsec")));
  differential->SetDirectory(nullptr);
  for (int bin = 1; bin <= kNBins; ++bin) {
    const double scale = 1. / (luminosity * branching * differential->GetBinWidth(bin));
    differential->SetBinContent(bin, differential->GetBinContent(bin) * scale);
    differential->SetBinError(bin, differential->GetBinError(bin) * scale);
  }
  differential->SetTitle("Unfolded inclusive t#bar{t} cross section;p_{T}(t_{had}) [GeV];d#sigma/dp_{T} [pb/GeV]");
  auto normalizedMatrix = NormalizeColumns(*matrix);
  normalizedMatrix->SetTitle("Response matrix;generated p_{T} [GeV];reconstructed p_{T} [GeV]");

  const std::string directory = "output/unfolding";
  gSystem->mkdir(directory.c_str(), true);
  std::ofstream report(directory + "/unfolding_report.txt");
  report << std::setprecision(7)
         << "Simplified Bayesian unfolding of hadronic-top pT\n"
         << "================================================\n"
         << "Files data/TT/W/DY: " << dataFiles.size() << "/" << ttFiles.size()
         << "/" << wFiles.size() << "/" << dyFiles.size() << "\n"
         << "Iterations: " << iterations << "\n"
         << "Maximum events per sample: " << maximumEventsPerSample << "\n"
         << "Negative data-background bins set to zero: " << negativeBins << "\n\n"
         << "bin [GeV]       efficiency   closure ratio   dSigma/dpT [pb/GeV]\n";
  for (int bin = 1; bin <= kNBins; ++bin) {
    double selected = 0.;
    for (int reco = 1; reco <= kNBins; ++reco)
      selected += matrix->GetBinContent(bin, reco);
    const double efficiency = truth->GetBinContent(bin) != 0.
        ? selected / truth->GetBinContent(bin) : 0.;
    const double ratio = closureTruth->GetBinContent(bin) != 0.
        ? closure.histogram->GetBinContent(bin) / closureTruth->GetBinContent(bin) : 0.;
    report << kEdges[bin - 1] << "-" << kEdges[bin] << "        "
           << efficiency << "       " << ratio << "       "
           << differential->GetBinContent(bin) << " +- "
           << differential->GetBinError(bin) << "\n";
  }
  report << "\nErrors are approximate statistical errors only.\n";

  TCanvas responseCanvas("responseCanvas", "Response", 760, 650);
  responseCanvas.SetRightMargin(0.15);
  normalizedMatrix->Draw("COLZ TEXT");
  responseCanvas.SaveAs((directory + "/response_matrix.png").c_str());
  closureTruth->SetLineColor(kBlack);
  closureTruth->SetLineWidth(2);
  closure.histogram->SetLineColor(kRed + 1);
  closure.histogram->SetMarkerStyle(20);
  closureTruth->SetTitle("Closure test;p_{T}(t_{had}) [GeV];weighted events");
  TCanvas closureCanvas("closureCanvas", "Closure", 760, 650);
  closureTruth->Draw("HIST");
  closure.histogram->Draw("E1 SAME");
  TLegend legend(0.57, 0.73, 0.88, 0.87);
  legend.AddEntry(closureTruth.get(), "Independent MC truth", "l");
  legend.AddEntry(closure.histogram.get(), "Unfolded MC reco", "lep");
  legend.Draw();
  closureCanvas.SaveAs((directory + "/closure_test.png").c_str());
  TCanvas resultCanvas("resultCanvas", "Unfolded result", 760, 650);
  differential->SetMarkerStyle(20);
  differential->Draw("E1");
  resultCanvas.SaveAs((directory + "/unfolded_differential_xsec.png").c_str());

  TFile output((directory + "/unfolding.root").c_str(), "RECREATE");
  data->Write(); wHistogram->Write(); dyHistogram->Write(); background->Write();
  signalReco->Write(); truth->Write(); matrix->Write(); normalizedMatrix->Write();
  trainingTruth->Write(); trainingMatrix->Write(); closureTruth->Write();
  closureReco->Write(); closure.histogram->Write(); unfolded.histogram->Write();
  differential->Write();
  output.Close();
  std::cout << "Unfolding complete: " << directory << "/unfolding_report.txt\n";
  return 0;
}
