#include "unfoldTopPt.C"

#include <TNamed.h>
#include <TParameter.h>

#include <cstdio>
#include <memory>
#include <string>

// Process one NanoAOD file and write a ROOT checkpoint containing only
// the histograms needed by the final unfolding.
int processUnfoldingFile(const char *fileUrl, const char *sampleType,
                         const char *outputPath, Long64_t maximumEvents = -1) {
  using namespace TopPtUnfolding;
  const std::string url = fileUrl ? fileUrl : "";
  const std::string type = sampleType ? sampleType : "";
  const std::string output = outputPath ? outputPath : "";
  if (url.empty() || output.empty() ||
      (type != "data" && type != "tt" && type != "w" && type != "dy")) {
    std::cerr << "ERROR: invalid processUnfoldingFile arguments.\n";
    return 1;
  }

  ROOT::RDataFrame frame("Events", url);
  ROOT::RDF::RNode input = frame;
  if (maximumEvents > 0) input = input.Range(maximumEvents);
  ROOT::RDF::TH1DModel ptModel("pt", "", kNBins, kEdges);
  ROOT::RDF::TH2DModel matrixModel("matrix", "", kNBins, kEdges,
                                   kNBins, kEdges);
  const std::string temporary = output + ".tmp";
  TFile checkpoint(temporary.c_str(), "RECREATE");
  if (checkpoint.IsZombie()) return 2;

  if (type == "data") {
    auto mask = std::make_shared<LumiMask>(
        "config/filelists/Cert_271036-284044_13TeV_Legacy2016_Collisions16_JSON.txt");
    auto certified = input.Filter(
        [mask](unsigned int run, unsigned int luminosityBlock) {
          return mask->Accept(run, luminosityBlock);
        }, {"run", "luminosityBlock"});
    auto reco = AddRecoTop(SelectTopLike(certified));
    auto histogram = reco.Histo1D(ptModel, "RecoTopPt");
    histogram->SetName("h_data_reco");
    histogram->Write();
  } else if (type == "tt") {
    auto weighted = input.Define("analysisWeight", "double(genWeight)");
    auto truth = AddGenTop(weighted);
    auto responseEvents = AddGenTop(AddRecoTop(SelectTopLike(weighted)));
    auto hTruth = truth.Histo1D(ptModel, "GenTopPt", "analysisWeight");
    auto hReco = responseEvents.Histo1D(ptModel, "RecoTopPt", "analysisWeight");
    auto hMatrix = responseEvents.Histo2D(
        matrixModel, "GenTopPt", "RecoTopPt", "analysisWeight");
    hTruth->SetName("h_tt_truth");
    hReco->SetName("h_tt_reco");
    hMatrix->SetName("h_tt_matrix");
    hTruth->Write();
    hReco->Write();
    hMatrix->Write();
  } else if (type == "dy") {
    auto weighted = input.Define("analysisWeight", "double(genWeight)");
    auto reco = AddRecoTop(SelectTopLike(weighted));
    auto total = weighted.Sum<double>("analysisWeight");
    auto histogram = reco.Histo1D(ptModel, "RecoTopPt", "analysisWeight");
    histogram->SetName("h_dy_reco_raw");
    histogram->Write();
    TParameter<double>("dy_total_sumw", *total).Write();
  } else {
    auto weighted = input.Define("analysisWeight", "double(genWeight)");
    auto reco = AddRecoTop(SelectTopLike(weighted));
    for (int category = 0; category < 3; ++category) {
      const std::string cut = "LHE_Njets == " + std::to_string(category);
      auto total = weighted.Filter(cut).Sum<double>("analysisWeight");
      auto histogram = reco.Filter(cut).Histo1D(
          ptModel, "RecoTopPt", "analysisWeight");
      histogram->SetName(("h_w" + std::to_string(category) + "_reco_raw").c_str());
      histogram->Write();
      TParameter<double>(("w" + std::to_string(category) + "_total_sumw").c_str(),
                         *total).Write();
    }
  }

  TNamed("checkpoint_status", "ok").Write();
  checkpoint.Close();
  if (std::rename(temporary.c_str(), output.c_str()) != 0) {
    std::cerr << "ERROR: cannot finalize " << output << "\n";
    return 3;
  }
  std::cout << "Unfolding checkpoint written: " << output << "\n";
  return 0;
}
