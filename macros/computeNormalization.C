#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Normalization {

struct Category {
  std::string name;
  double crossSectionPb = 0.0;
  ULong64_t eventCount = 0;
  double sumWeights = 0.0;
  double sumWeights2 = 0.0;
  double scaleFactor = 0.0;
};

std::string Trim(const std::string &input) {
  const auto first = std::find_if_not(input.begin(), input.end(),
                                      [](unsigned char c) { return std::isspace(c); });
  const auto last = std::find_if_not(input.rbegin(), input.rend(),
                                     [](unsigned char c) { return std::isspace(c); }).base();
  return first < last ? std::string(first, last) : std::string();
}

std::map<std::string, double> ReadConfiguration(const char *path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error(std::string("Cannot open configuration: ") + path);

  std::map<std::string, double> values;
  std::string line;
  while (std::getline(input, line)) {
    line = Trim(line);
    if (line.empty() || line.front() == '#') continue;
    const auto separator = line.find('=');
    if (separator == std::string::npos)
      throw std::runtime_error("Invalid configuration line: " + line);
    values[Trim(line.substr(0, separator))] =
        std::stod(Trim(line.substr(separator + 1)));
  }
  return values;
}

double Require(const std::map<std::string, double> &configuration,
               const std::string &key) {
  const auto found = configuration.find(key);
  if (found == configuration.end())
    throw std::runtime_error("Missing configuration key: " + key);
  return found->second;
}

double SumRunsLeaf(TTree *runs, const char *branchName) {
  if (!runs || !runs->GetBranch(branchName)) return 0.0;
  TLeaf *leaf = runs->GetLeaf(branchName);
  if (!leaf) return 0.0;
  double sum = 0.0;
  for (Long64_t entry = 0; entry < runs->GetEntries(); ++entry) {
    runs->GetEntry(entry);
    sum += leaf->GetValue();
  }
  return sum;
}

void Finalize(Category &category, double luminosityPb,
              double branchingFactor) {
  if (category.sumWeights == 0.0 || category.crossSectionPb == 0.0) {
    category.scaleFactor = 0.0;
    return;
  }
  category.scaleFactor = luminosityPb * category.crossSectionPb *
                         branchingFactor / category.sumWeights;
}

}  // namespace Normalization

int computeNormalization(const char *fileUrl,
                         const char *sampleLabel,
                         Long64_t maxEvents = -1,
                         const char *configurationPath =
                             "config/normalization.cfg") {
  using namespace Normalization;

  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);

  if (!fileUrl || !sampleLabel) {
    std::cerr << "ERROR: file URL and sample label are required.\n";
    return 1;
  }

  const std::string label(sampleLabel);
  const bool isWJets = label == "wjets";
  if (!isWJets && label != "ttbar_semileptonic" && label != "dyjets") {
    std::cerr << "ERROR: unsupported Monte Carlo label: " << label << "\n";
    return 2;
  }

  std::map<std::string, double> configuration;
  try {
    configuration = ReadConfiguration(configurationPath);
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 3;
  }

  const double luminosityPb = Require(configuration, "luminosity_pb");
  const double branchingFactor = Require(configuration, "branching_factor");

  std::cout << "Opening " << fileUrl << "\n";
  std::unique_ptr<TFile> inputFile(TFile::Open(fileUrl, "READ"));
  if (!inputFile || inputFile->IsZombie()) {
    std::cerr << "ERROR: ROOT could not open the input file.\n";
    return 4;
  }
  auto *events = dynamic_cast<TTree *>(inputFile->Get("Events"));
  auto *runs = dynamic_cast<TTree *>(inputFile->Get("Runs"));
  if (!events || !runs || !events->GetBranch("genWeight")) {
    std::cerr << "ERROR: required Monte Carlo trees or branches are missing.\n";
    return 5;
  }
  if (isWJets && !events->GetBranch("LHE_Njets")) {
    std::cerr << "ERROR: LHE_Njets is required for W+jets categories.\n";
    return 6;
  }

  std::vector<Category> categories;
  if (label == "ttbar_semileptonic") {
    categories.push_back({"ttbar_semileptonic",
                          Require(configuration, "ttbar_semileptonic_xsec_pb")});
  } else if (label == "dyjets") {
    categories.push_back({"dyjets", Require(configuration, "dyjets_xsec_pb")});
  } else {
    categories.push_back({"W0J", Require(configuration, "w0j_xsec_pb")});
    categories.push_back({"W1J", Require(configuration, "w1j_xsec_pb")});
    categories.push_back({"W2J", Require(configuration, "w2j_xsec_pb")});
    categories.push_back({"W3plusJ_unscaled", 0.0});
  }

  TTreeReader reader(events);
  TTreeReaderValue<Float_t> genWeight(reader, "genWeight");
  std::unique_ptr<TTreeReaderValue<UChar_t>> lheNjets;
  if (isWJets)
    lheNjets = std::make_unique<TTreeReaderValue<UChar_t>>(reader, "LHE_Njets");

  const Long64_t availableEntries = events->GetEntries();
  const Long64_t entriesToRead = maxEvents > 0
      ? std::min(maxEvents, availableEntries) : availableEntries;
  Long64_t entry = 0;
  while (entry < entriesToRead && reader.Next()) {
    std::size_t categoryIndex = 0;
    if (isWJets) {
      const unsigned int numberLheJets = static_cast<unsigned int>(**lheNjets);
      categoryIndex = numberLheJets <= 2 ? numberLheJets : 3;
    }
    Category &category = categories.at(categoryIndex);
    const double weight = static_cast<double>(*genWeight);
    ++category.eventCount;
    category.sumWeights += weight;
    category.sumWeights2 += weight * weight;
    ++entry;
  }

  for (auto &category : categories)
    Finalize(category, luminosityPb, branchingFactor);

  const double runsEventCount = SumRunsLeaf(runs, "genEventCount");
  const double runsSumWeights = SumRunsLeaf(runs, "genEventSumw");
  double scannedSumWeights = 0.0;
  for (const auto &category : categories) scannedSumWeights += category.sumWeights;

  const std::string outputDirectory = "output/phase4";
  gSystem->mkdir(outputDirectory.c_str(), true);
  const std::string reportPath = outputDirectory + "/" + label + "_normalization.txt";
  std::ofstream report(reportPath);
  report << "Monte Carlo normalization report\n";
  report << "================================\n";
  report << "Sample: " << label << "\n";
  report << "File: " << fileUrl << "\n";
  report << "Processed events: " << entry << " / " << availableEntries << "\n";
  report << "Target luminosity: " << luminosityPb << " pb^-1 = "
         << luminosityPb / 1000.0 << " fb^-1\n";
  report << "Additional branching factor: " << branchingFactor << "\n\n";
  report << "Event weight: normFactor * genWeight\n";
  report << "normFactor = luminosity * crossSection / sumGenWeights\n\n";
  report << std::left << std::setw(22) << "Category" << std::right
         << std::setw(13) << "Events" << std::setw(16) << "xsec [pb]"
         << std::setw(20) << "sumGenWeights" << std::setw(18) << "N_eff"
         << std::setw(20) << "normFactor" << "\n";
  report << std::string(109, '-') << "\n";
  report << std::setprecision(10);
  for (const auto &category : categories) {
    const double effectiveEvents = category.sumWeights2 > 0.0
        ? category.sumWeights * category.sumWeights / category.sumWeights2 : 0.0;
    report << std::left << std::setw(22) << category.name << std::right
           << std::setw(13) << category.eventCount
           << std::setw(16) << category.crossSectionPb
           << std::setw(20) << category.sumWeights
           << std::setw(18) << effectiveEvents
           << std::setw(20) << category.scaleFactor << "\n";
  }
  report << "\nRuns/genEventCount for complete file: " << runsEventCount << "\n";
  report << "Runs/genEventSumw for complete file: " << runsSumWeights << "\n";
  report << "Scanned sumGenWeights: " << scannedSumWeights << "\n";
  if (entry == availableEntries) {
    const double relativeDifference = runsSumWeights != 0.0
        ? (scannedSumWeights - runsSumWeights) / runsSumWeights : 0.0;
    report << "Relative Events-vs-Runs sumw difference: "
           << relativeDifference << "\n";
  } else {
    report << "Events-vs-Runs comparison skipped: only a subset was scanned.\n";
  }
  report << "\nImportant: these factors normalize the processed event subset to the\n"
            "target luminosity. Recompute them whenever the processed file list\n"
            "or maximum event count changes.\n";

  auto hSumWeights = std::make_unique<TH1D>(
      "h_sumGenWeights", "Generator-weight sums", categories.size(),
      0., categories.size());
  auto hNormFactors = std::make_unique<TH1D>(
      "h_normFactors", "Monte Carlo normalization factors", categories.size(),
      0., categories.size());
  hSumWeights->SetDirectory(nullptr);
  hNormFactors->SetDirectory(nullptr);
  for (std::size_t index = 0; index < categories.size(); ++index) {
    hSumWeights->SetBinContent(index + 1, categories[index].sumWeights);
    hNormFactors->SetBinContent(index + 1, categories[index].scaleFactor);
    hSumWeights->GetXaxis()->SetBinLabel(index + 1, categories[index].name.c_str());
    hNormFactors->GetXaxis()->SetBinLabel(index + 1, categories[index].name.c_str());
  }
  hSumWeights->GetYaxis()->SetTitle("sumGenWeights");
  hNormFactors->GetYaxis()->SetTitle("normFactor");
  hSumWeights->SetFillColor(kAzure - 9);
  hSumWeights->SetLineColor(kAzure + 2);
  hNormFactors->SetFillColor(kOrange - 2);
  hNormFactors->SetLineColor(kOrange + 7);
  hSumWeights->SetStats(false);
  hNormFactors->SetStats(false);

  TCanvas canvas("canvas", "Monte Carlo normalization", 1100, 500);
  canvas.Divide(2, 1);
  canvas.cd(1); gPad->SetLeftMargin(0.18); gPad->SetLogy(); hSumWeights->Draw("HIST");
  canvas.cd(2); gPad->SetLeftMargin(0.18); gPad->SetLogy(); hNormFactors->Draw("HIST");
  const std::string figurePath = outputDirectory + "/" + label + "_normalization.png";
  canvas.SaveAs(figurePath.c_str());

  const std::string rootPath = outputDirectory + "/" + label + "_normalization.root";
  TFile outputFile(rootPath.c_str(), "RECREATE");
  hSumWeights->Write();
  hNormFactors->Write();
  outputFile.Close();

  std::cout << "Processed " << entry << " events.\n";
  std::cout << "Report: " << reportPath << "\n";
  std::cout << "Figure: " << figurePath << "\n";
  return 0;
}

