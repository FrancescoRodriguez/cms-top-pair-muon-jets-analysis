#include "measureInclusive.C"

#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TSystem.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>

namespace CheckpointCombination {

std::map<std::string, std::string> ReadCheckpoint(const std::string &path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Missing checkpoint: " + path);
  std::map<std::string, std::string> values;
  std::string line;
  while (std::getline(input, line)) {
    const auto separator = line.find('=');
    if (separator != std::string::npos)
      values[line.substr(0, separator)] = line.substr(separator + 1);
  }
  if (values["status"] != "ok")
    throw std::runtime_error("Incomplete checkpoint: " + path);
  return values;
}

int CountInputFiles(const std::string &path) {
  std::ifstream input(path);
  if (!input) throw std::runtime_error("Missing file list: " + path);
  int count = 0;
  std::string line;
  while (std::getline(input, line))
    if (!line.empty() && line.front() != '#') ++count;
  return count;
}

std::string CheckpointPath(const std::string &directory,
                           const std::string &sample, int index) {
  char buffer[64];
  std::snprintf(buffer, sizeof(buffer), "/%s_%04d.txt", sample.c_str(), index);
  return directory + buffer;
}

double Number(const std::map<std::string, std::string> &values,
              const std::string &key) {
  const auto found = values.find(key);
  if (found == values.end()) throw std::runtime_error("Missing key: " + key);
  return std::stod(found->second);
}

}  // namespace CheckpointCombination

int combineCheckpoints(const char *checkpointDirectory = "output/checkpoints",
                       int maxTTFiles = -1, int maxWFiles = -1,
                       int maxDYFiles = -1) {
  using namespace CheckpointCombination;
  using namespace InclusiveMeasurement;

  const std::string directory(checkpointDirectory);
  const std::map<std::string, std::string> lists = {
      {"dataG", "config/filelists/data_Run2016G.txt"},
      {"dataH", "config/filelists/data_Run2016H.txt"},
      {"tt", "config/filelists/ttbar_semileptonic.txt"},
      {"w", "config/filelists/wjets.txt"},
      {"dy", "config/filelists/dyjets.txt"}};

  double dataYield = 0.0;
  double ttTotal = 0.0, ttSelected = 0.0, ttSelected2 = 0.0;
  double dyTotal = 0.0, dySelected = 0.0, dySelected2 = 0.0;
  double wTotal[3] = {0.0, 0.0, 0.0};
  double wSelected[3] = {0.0, 0.0, 0.0};
  double wSelected2[3] = {0.0, 0.0, 0.0};

  try {
    for (const auto &[sample, list] : lists) {
      int numberOfFiles = CountInputFiles(list);
      if (sample == "tt" && maxTTFiles > 0)
        numberOfFiles = std::min(numberOfFiles, maxTTFiles);
      if (sample == "w" && maxWFiles > 0)
        numberOfFiles = std::min(numberOfFiles, maxWFiles);
      if (sample == "dy" && maxDYFiles > 0)
        numberOfFiles = std::min(numberOfFiles, maxDYFiles);
      for (int index = 1; index <= numberOfFiles; ++index) {
        const auto values = ReadCheckpoint(CheckpointPath(directory, sample, index));
        if (sample == "dataG" || sample == "dataH") {
          dataYield += Number(values, "selectedData");
        } else if (sample == "tt") {
          ttTotal += Number(values, "totalSumW0");
          ttSelected += Number(values, "selectedSumW0");
          ttSelected2 += Number(values, "selectedSumW2_0");
        } else if (sample == "dy") {
          dyTotal += Number(values, "totalSumW0");
          dySelected += Number(values, "selectedSumW0");
          dySelected2 += Number(values, "selectedSumW2_0");
        } else {
          for (int category = 0; category < 3; ++category) {
            wTotal[category] += Number(values, "totalSumW" + std::to_string(category));
            wSelected[category] += Number(values, "selectedSumW" + std::to_string(category));
            wSelected2[category] += Number(values, "selectedSumW2_" + std::to_string(category));
          }
        }
      }
    }
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 1;
  }

  const auto config = ReadConfiguration("config/normalization.cfg");
  const double luminosity = Required(config, "luminosity_pb");
  const double branching = Required(config, "ttbar_semileptonic_br");
  const double ttXsec = Required(config, "ttbar_semileptonic_xsec_pb");
  const double dyXsec = Required(config, "dyjets_xsec_pb");
  const double reference = Required(config, "ttbar_inclusive_reference_pb");
  const double wXsec[3] = {Required(config, "w0j_xsec_pb"),
                           Required(config, "w1j_xsec_pb"),
                           Required(config, "w2j_xsec_pb")};

  const double efficiency = ttTotal != 0.0 ? ttSelected / ttTotal : 0.0;
  const double expectedTT = luminosity * ttXsec * efficiency;
  const double expectedDY = dyTotal != 0.0
      ? luminosity * dyXsec * dySelected / dyTotal : 0.0;
  double expectedW = 0.0;
  double backgroundVariance = dyTotal != 0.0
      ? std::pow(luminosity * dyXsec / dyTotal, 2) * dySelected2 : 0.0;
  for (int category = 0; category < 3; ++category) {
    if (wTotal[category] == 0.0) continue;
    expectedW += luminosity * wXsec[category] *
                 wSelected[category] / wTotal[category];
    backgroundVariance += std::pow(luminosity * wXsec[category] /
                                    wTotal[category], 2) * wSelected2[category];
  }
  const double background = expectedW + expectedDY;
  const double signal = dataYield - background;
  const double semileptonicXsec = signal / (luminosity * efficiency);
  const double inclusiveXsec = semileptonicXsec / branching;
  const double dataStat = std::sqrt(dataYield) /
                          (luminosity * efficiency * branching);
  const double backgroundStat = std::sqrt(backgroundVariance) /
                                (luminosity * efficiency * branching);

  gSystem->mkdir("output/phase6", true);
  std::ofstream report("output/phase6/inclusive_measurement.txt");
  report << std::setprecision(10);
  report << "Simplified inclusive ttbar cross-section measurement\n";
  report << "====================================================\n";
  const bool mcSubset = maxTTFiles > 0 || maxWFiles > 0 || maxDYFiles > 0;
  report << "Mode: FULL DATA, "
         << (mcSubset ? "MONTE CARLO SUBSET" : "FULL MONTE CARLO") << "\n";
  report << "MC files used: TT=" << (maxTTFiles > 0 ? maxTTFiles : 138)
         << ", W=" << (maxWFiles > 0 ? maxWFiles : 68)
         << ", DY=" << (maxDYFiles > 0 ? maxDYFiles : 61) << "\n";
  report << "Luminosity: " << luminosity << " pb^-1\n";
  report << "TT semileptonic efficiency: " << efficiency << "\n";
  report << "Semileptonic branching fraction: " << branching << "\n\n";
  report << "Selected data: " << dataYield << "\n";
  report << "Expected W+jets: " << expectedW << "\n";
  report << "Expected DY: " << expectedDY << "\n";
  report << "Total background: " << background << "\n";
  report << "Background-subtracted signal: " << signal << "\n";
  report << "Expected TT signal: " << expectedTT << "\n";
  report << "Data signal / expected TT: " << signal / expectedTT << "\n\n";
  report << "Measured semileptonic cross section: " << semileptonicXsec << " pb\n";
  report << "Measured inclusive cross section: " << inclusiveXsec << " pb\n";
  report << "Data statistical uncertainty: " << dataStat << " pb\n";
  report << "Background MC statistical uncertainty: " << backgroundStat << " pb\n";
  report << "Reference inclusive cross section: " << reference << " pb\n\n";
  report << "CAVEAT: only the W+jets and DY backgrounds supplied for this\n"
            "simplified exam project are subtracted.\n";
  report.close();

  auto yields = std::make_unique<TH1D>("h_yields", "Selected event yields", 4, 0, 4);
  yields->GetXaxis()->SetBinLabel(1, "Data");
  yields->GetXaxis()->SetBinLabel(2, "W+jets");
  yields->GetXaxis()->SetBinLabel(3, "DY");
  yields->GetXaxis()->SetBinLabel(4, "Data-bkg");
  yields->SetBinContent(1, dataYield);
  yields->SetBinContent(2, expectedW);
  yields->SetBinContent(3, expectedDY);
  yields->SetBinContent(4, signal);
  TCanvas canvas("canvas", "yields", 850, 650);
  yields->Draw("HIST");
  canvas.SaveAs("output/phase6/inclusive_yields.png");
  TFile rootOutput("output/phase6/inclusive_measurement.root", "RECREATE");
  yields->Write();
  rootOutput.Close();

  std::cout << "Inclusive ttbar cross section: " << inclusiveXsec << " pb\n";
  return 0;
}
