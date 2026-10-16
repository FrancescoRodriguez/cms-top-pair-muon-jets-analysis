#include "unfoldTopPt.C"

#include <TNamed.h>
#include <TParameter.h>

#include <cstdio>
#include <memory>
#include <stdexcept>

namespace ReliableUnfolding {

std::string Path(const std::string &directory, const std::string &sample,
                 int index) {
  char buffer[80];
  std::snprintf(buffer, sizeof(buffer), "/%s_%04d.root", sample.c_str(), index);
  return directory + buffer;
}

void Check(TFile &file, const std::string &path) {
  if (file.IsZombie()) throw std::runtime_error("Cannot open " + path);
  auto status = dynamic_cast<TNamed *>(file.Get("checkpoint_status"));
  if (!status || std::string(status->GetTitle()) != "ok")
    throw std::runtime_error("Incomplete checkpoint " + path);
}

template <class Histogram>
void AddHistogram(std::unique_ptr<Histogram> &sum, TFile &file,
                  const char *name) {
  auto input = dynamic_cast<Histogram *>(file.Get(name));
  if (!input) throw std::runtime_error(std::string("Missing histogram ") + name);
  if (!sum) {
    sum.reset(static_cast<Histogram *>(input->Clone(name)));
    sum->SetDirectory(nullptr);
  } else {
    sum->Add(input);
  }
}

double Parameter(TFile &file, const char *name) {
  auto value = dynamic_cast<TParameter<double> *>(file.Get(name));
  if (!value) throw std::runtime_error(std::string("Missing parameter ") + name);
  return value->GetVal();
}

}  // namespace ReliableUnfolding

int combineUnfolding(const char *checkpointDirectory =
                         "output/unfolding_checkpoints",
                     int numberOfTTFiles = 10,
                     int numberOfWFiles = 10,
                     int numberOfDYFiles = 10,
                     int iterations = 4) {
  using namespace ReliableUnfolding;
  using namespace TopPtUnfolding;
  gROOT->SetBatch(true);
  gStyle->SetOptStat(0);
  TH1::SetDefaultSumw2(true);
  if (numberOfTTFiles < 2 || numberOfWFiles < 1 ||
      numberOfDYFiles < 1 || iterations < 1) return 1;

  const std::string directory(checkpointDirectory);
  std::unique_ptr<TH1D> data, truth, closureTruth, closureReco;
  std::unique_ptr<TH2D> matrix, trainingMatrix;
  std::unique_ptr<TH1D> trainingTruth, dyRaw;
  std::unique_ptr<TH1D> wRaw[3];
  double dyTotal = 0.;
  double wTotal[3] = {0., 0., 0.};

  try {
    for (int index = 1; index <= 70; ++index) {
      const auto path = Path(directory, "dataG", index);
      TFile file(path.c_str(), "READ"); Check(file, path);
      AddHistogram(data, file, "h_data_reco");
    }
    for (int index = 1; index <= 82; ++index) {
      const auto path = Path(directory, "dataH", index);
      TFile file(path.c_str(), "READ"); Check(file, path);
      AddHistogram(data, file, "h_data_reco");
    }

    const int trainingFiles = numberOfTTFiles / 2;
    for (int index = 1; index <= numberOfTTFiles; ++index) {
      const auto path = Path(directory, "tt", index);
      TFile file(path.c_str(), "READ"); Check(file, path);
      AddHistogram(truth, file, "h_tt_truth");
      AddHistogram(matrix, file, "h_tt_matrix");
      if (index <= trainingFiles) {
        AddHistogram(trainingTruth, file, "h_tt_truth");
        AddHistogram(trainingMatrix, file, "h_tt_matrix");
      } else {
        AddHistogram(closureTruth, file, "h_tt_truth");
        AddHistogram(closureReco, file, "h_tt_reco");
      }
    }
    for (int index = 1; index <= numberOfWFiles; ++index) {
      const auto path = Path(directory, "w", index);
      TFile file(path.c_str(), "READ"); Check(file, path);
      for (int category = 0; category < 3; ++category) {
        const std::string histogram =
            "h_w" + std::to_string(category) + "_reco_raw";
        const std::string parameter =
            "w" + std::to_string(category) + "_total_sumw";
        AddHistogram(wRaw[category], file, histogram.c_str());
        wTotal[category] += Parameter(file, parameter.c_str());
      }
    }
    for (int index = 1; index <= numberOfDYFiles; ++index) {
      const auto path = Path(directory, "dy", index);
      TFile file(path.c_str(), "READ"); Check(file, path);
      AddHistogram(dyRaw, file, "h_dy_reco_raw");
      dyTotal += Parameter(file, "dy_total_sumw");
    }
  } catch (const std::exception &error) {
    std::cerr << "ERROR: " << error.what() << "\n";
    return 2;
  }

  const auto config = ReadConfig("config/normalization.cfg");
  const double luminosity = Get(config, "luminosity_pb");
  const double branching = Get(config, "ttbar_semileptonic_br");
  auto dy = std::unique_ptr<TH1D>(static_cast<TH1D *>(dyRaw->Clone("h_dy_reco")));
  dy->SetDirectory(nullptr);
  dy->Scale(luminosity * Get(config, "dyjets_xsec_pb") / dyTotal);
  auto w = std::unique_ptr<TH1D>(static_cast<TH1D *>(data->Clone("h_w_reco")));
  w->Reset("ICES");
  const double wXsec[3] = {Get(config, "w0j_xsec_pb"),
                           Get(config, "w1j_xsec_pb"),
                           Get(config, "w2j_xsec_pb")};
  for (int category = 0; category < 3; ++category) {
    auto component = std::unique_ptr<TH1D>(
        static_cast<TH1D *>(wRaw[category]->Clone()));
    component->SetDirectory(nullptr);
    if (wTotal[category] != 0.)
      component->Scale(luminosity * wXsec[category] / wTotal[category]);
    w->Add(component.get());
  }
  auto background = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(w->Clone("h_background")));
  background->Add(dy.get());
  auto signalReco = std::unique_ptr<TH1D>(
      static_cast<TH1D *>(data->Clone("h_signal_reco")));
  signalReco->Add(background.get(), -1.);
  int negativeBins = 0;
  for (int bin = 1; bin <= kNBins; ++bin) {
    if (signalReco->GetBinContent(bin) < 0.) {
      signalReco->SetBinContent(bin, 0.);
      ++negativeBins;
    }
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
    const double scale = 1. /
        (luminosity * branching * differential->GetBinWidth(bin));
    differential->SetBinContent(bin, differential->GetBinContent(bin) * scale);
    differential->SetBinError(bin, differential->GetBinError(bin) * scale);
  }
  differential->SetTitle(
      "Unfolded inclusive t#bar{t} cross section;p_{T}(t_{had}) [GeV];d#sigma/dp_{T} [pb/GeV]");
  auto normalized = NormalizeColumns(*matrix);
  normalized->SetTitle(
      "Response matrix;generated p_{T} [GeV];reconstructed p_{T} [GeV]");

  gSystem->mkdir("output/unfolding", true);
  std::ofstream report("output/unfolding/unfolding_report.txt");
  report << std::setprecision(7)
         << "Reliable Bayesian unfolding of hadronic-top pT\n"
         << "==============================================\n"
         << "Mode: FULL DATA, MONTE CARLO SUBSET\n"
         << "Files data/TT/W/DY: 152/" << numberOfTTFiles << "/"
         << numberOfWFiles << "/" << numberOfDYFiles << "\n"
         << "Iterations: " << iterations << "\n"
         << "Negative data-background bins set to zero: " << negativeBins << "\n\n"
         << "bin [GeV]       efficiency   closure ratio   dSigma/dpT [pb/GeV]\n";
  double integratedCrossSection = 0.;
  for (int bin = 1; bin <= kNBins; ++bin) {
    double selected = 0.;
    for (int reco = 1; reco <= kNBins; ++reco)
      selected += matrix->GetBinContent(bin, reco);
    const double efficiency = truth->GetBinContent(bin) != 0.
        ? selected / truth->GetBinContent(bin) : 0.;
    const double closureRatio = closureTruth->GetBinContent(bin) != 0.
        ? closure.histogram->GetBinContent(bin) /
          closureTruth->GetBinContent(bin) : 0.;
    integratedCrossSection += differential->GetBinContent(bin) *
                              differential->GetBinWidth(bin);
    report << kEdges[bin - 1] << "-" << kEdges[bin] << "        "
           << efficiency << "       " << closureRatio << "       "
           << differential->GetBinContent(bin) << " +- "
           << differential->GetBinError(bin) << "\n";
  }
  report << "\nIntegrated unfolded cross section: "
         << integratedCrossSection << " pb\n"
         << "Errors are approximate statistical errors only.\n";

  TCanvas responseCanvas("responseCanvas", "Response", 760, 650);
  responseCanvas.SetRightMargin(0.15);
  normalized->Draw("COLZ TEXT");
  responseCanvas.SaveAs("output/unfolding/response_matrix.png");
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
  closureCanvas.SaveAs("output/unfolding/closure_test.png");
  TCanvas resultCanvas("resultCanvas", "Unfolded", 760, 650);
  differential->SetMarkerStyle(20);
  differential->Draw("E1");
  resultCanvas.SaveAs("output/unfolding/unfolded_differential_xsec.png");

  TFile output("output/unfolding/unfolding.root", "RECREATE");
  data->Write(); w->Write(); dy->Write(); background->Write(); signalReco->Write();
  truth->Write(); matrix->Write(); normalized->Write(); trainingTruth->Write();
  trainingMatrix->Write(); closureTruth->Write(); closureReco->Write();
  closure.histogram->Write(); unfolded.histogram->Write(); differential->Write();
  output.Close();
  std::cout << "Integrated unfolded cross section: "
            << integratedCrossSection << " pb\n";
  return 0;
}
