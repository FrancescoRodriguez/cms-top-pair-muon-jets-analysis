#include "measureInclusive.C"

#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>

// Process one remote NanoAOD file and save only the few numbers needed by the
// inclusive measurement. A completed text file is an elementary checkpoint.
int processOneFile(const char *fileUrl, const char *sampleType,
                   const char *outputPath, Long64_t maxEvents = -1) {
  using namespace InclusiveMeasurement;

  const std::string url = fileUrl ? fileUrl : "";
  const std::string type = sampleType ? sampleType : "";
  const std::string output = outputPath ? outputPath : "";
  if (url.empty() || output.empty() ||
      (type != "data" && type != "tt" && type != "w" && type != "dy")) {
    std::cerr << "ERROR: invalid processOneFile arguments.\n";
    return 1;
  }

  ROOT::RDataFrame frame("Events", url);
  ROOT::RDF::RNode input = Limit(frame, maxEvents);

  double selectedData = 0.0;
  double totalSumW[3] = {0.0, 0.0, 0.0};
  double selectedSumW[3] = {0.0, 0.0, 0.0};
  double selectedSumW2[3] = {0.0, 0.0, 0.0};

  if (type == "data") {
    auto mask = std::make_shared<LumiMask>(
        "config/filelists/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt");
    auto certified = input.Filter(
        [mask](unsigned int run, unsigned int luminosityBlock) {
          return mask->Accept(run, luminosityBlock);
        }, {"run", "luminosityBlock"});
    selectedData = static_cast<double>(*ApplyTopLikeSelection(certified).Count());
  } else {
    auto weighted = input.Define("analysisWeight", "double(genWeight)");
    auto selected = ApplyTopLikeSelection(weighted)
        .Define("analysisWeight2", "analysisWeight * analysisWeight");
    if (type == "w") {
      for (int category = 0; category < 3; ++category) {
        const std::string cut = "LHE_Njets == " + std::to_string(category);
        totalSumW[category] =
            *weighted.Filter(cut).Sum<double>("analysisWeight");
        auto selectedCategory = selected.Filter(cut);
        selectedSumW[category] =
            *selectedCategory.Sum<double>("analysisWeight");
        selectedSumW2[category] =
            *selectedCategory.Sum<double>("analysisWeight2");
      }
    } else {
      totalSumW[0] = *weighted.Sum<double>("analysisWeight");
      selectedSumW[0] = *selected.Sum<double>("analysisWeight");
      selectedSumW2[0] = *selected.Sum<double>("analysisWeight2");
    }
  }

  // Write atomically: an interrupted job never leaves a file that looks done.
  const std::string temporary = output + ".tmp";
  std::ofstream report(temporary);
  if (!report) {
    std::cerr << "ERROR: cannot create " << temporary << "\n";
    return 2;
  }
  report << std::setprecision(17);
  report << "status=ok\n";
  report << "sample=" << type << "\n";
  report << "selectedData=" << selectedData << "\n";
  for (int category = 0; category < 3; ++category) {
    report << "totalSumW" << category << "=" << totalSumW[category] << "\n";
    report << "selectedSumW" << category << "=" << selectedSumW[category] << "\n";
    report << "selectedSumW2_" << category << "=" << selectedSumW2[category]
           << "\n";
  }
  report.close();
  if (std::rename(temporary.c_str(), output.c_str()) != 0) {
    std::cerr << "ERROR: cannot finalize " << output << "\n";
    return 3;
  }
  std::cout << "Checkpoint written: " << output << "\n";
  return 0;
}
