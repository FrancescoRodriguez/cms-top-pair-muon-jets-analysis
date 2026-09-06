#include <TBranch.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1D.h>
#include <TLeaf.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>


#include <algorithm>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace {

struct BranchGroup {  
  std::string description;
  std::vector<std::string> names;
};

bool HasBranch(TTree *tree, const std::string &name) { 
  return tree && tree->GetBranch(name.c_str()); 
}

double SumLeaf(TTree *tree, const char *branchName) { 
  if (!tree || !tree->GetBranch(branchName)) return 0.0; 

  TLeaf *leaf = tree->GetLeaf(branchName);  
  if (!leaf) return 0.0; 

  double sum = 0.0; 
  for (Long64_t entry = 0; entry < tree->GetEntries(); ++entry) { 
    tree->GetEntry(entry); 
    sum += leaf->GetValue(); 
  }
  return sum;
}


void PrintBranchTable(TTree *tree, std::ostream &output) { 
  const std::vector<BranchGroup> groups = { 
      {"Event identifiers", {"run", "luminosityBlock", "event"}}, 
      {"Single-muon triggers", {"HLT_IsoMu24", "HLT_IsoTkMu24"}}, 
      {"Muon multiplicity and kinematics",
       {"nMuon", "Muon_pt", "Muon_eta", "Muon_phi", "Muon_mass"}}, 
      {"Muon identification and isolation",
       {"Muon_tightId", "Muon_mediumId", "Muon_pfRelIso04_all"}}, 
      {"Electron veto",
       {"nElectron", "Electron_pt", "Electron_eta", "Electron_cutBased"}}, 
      {"Jet reconstruction",
       {"nJet", "Jet_pt", "Jet_eta", "Jet_phi", "Jet_mass", "Jet_jetId"}}, 
      {"b tagging", {"Jet_btagDeepFlavB", "Jet_btagDeepB"}}, 
      {"Missing transverse momentum", {"MET_pt", "MET_phi"}}, 
      {"Monte Carlo normalization", {"genWeight", "LHE_Njets", "Pileup_nTrueInt"}}, 
      {"Generator-level objects",
       {"nGenPart", "GenPart_pt", "GenPart_eta", "GenPart_phi",
        "GenPart_mass", "GenPart_pdgId", "GenPart_status",
        "GenPart_statusFlags", "GenPart_genPartIdxMother"}}}; 

  output << "\nRequired branch check\n"; 
  output << "---------------------\n";
  for (const auto &group : groups) {        
    output << group.description << ":\n";   
    for (const auto &name : group.names) {   
      output << "  " << std::left << std::setw(30) << name
             << (HasBranch(tree, name) ? "available" : "MISSING") << "\n"; 
    }
  }
}

void StyleHistogram(TH1 *histogram, const char *xTitle) { 
  histogram->SetLineColor(kAzure + 2); 
  histogram->SetFillColorAlpha(kAzure - 9, 0.55); 
  histogram->SetLineWidth(2); 
  histogram->GetXaxis()->SetTitle(xTitle); 
  histogram->GetYaxis()->SetTitle("Entries");
  histogram->SetStats(false); 
}

Long64_t FillFromTree(TTree *tree, TH1 *histogram, const char *expression,
                      Long64_t entriesToRead) {
  if (!tree || !histogram) return 0;
  tree->SetEstimate(std::max<Long64_t>(tree->GetEstimate(), entriesToRead * 20));
  const Long64_t selectedRows =
      tree->Draw(expression, "", "goff", entriesToRead);
  const double *values = tree->GetV1();
  if (selectedRows < 0) return 0;
  for (Long64_t row = 0; row < selectedRows; ++row)
    histogram->Fill(values[row]);
  return selectedRows;
}

}  // end namespace. 

int inspectNanoAOD(const char *fileUrl,
                   const char *sampleLabel = "sample",
                   Long64_t maxEvents = 20000) {
  gROOT->SetBatch(true); 
  gStyle->SetOptStat(0); 
  TH1::SetDefaultSumw2(true); 

  if (!fileUrl || std::string(fileUrl).empty()) { 
    std::cerr << "ERROR: the input file URL is empty.\n"; 
    return 1;
  }

  const std::string label = sampleLabel ? sampleLabel : "sample"; 
  const std::string outputDirectory = "output/phase1"; 
  gSystem->mkdir(outputDirectory.c_str(), true); 

  // open ROOT file
  std::cout << "Opening " << fileUrl << "\n"; 
  std::unique_ptr<TFile> inputFile(TFile::Open(fileUrl, "READ"));  

  if (!inputFile || inputFile->IsZombie()) { 
    std::cerr << "ERROR: ROOT could not open the input file.\n"; 
    return 2;
  }

  auto *events = dynamic_cast<TTree *>(inputFile->Get("Events")); 
  auto *runs = dynamic_cast<TTree *>(inputFile->Get("Runs")); 
  auto *luminosityBlocks = 
      dynamic_cast<TTree *>(inputFile->Get("LuminosityBlocks"));
  if (!events || !runs || !luminosityBlocks) {
    std::cerr << "ERROR: one or more standard NanoAOD trees are missing.\n"; 
    return 3;
  }

  const std::string reportPath = outputDirectory + "/" + label + "_report.txt"; 
  std::ofstream report(reportPath); 
  if (!report) { 
    std::cerr << "ERROR: could not create " << reportPath << "\n";
    return 4;
  }

  report << "NanoAOD inspection report\n"; 
  report << "=========================\n";
  report << "Sample: " << label << "\n"; 
  report << "File: " << fileUrl << "\n\n"; 
  report << "Events entries: " << events->GetEntries() << "\n"; 
  report << "Runs entries: " << runs->GetEntries() << "\n"; 
  report << "LuminosityBlocks entries: " << luminosityBlocks->GetEntries() << "\n"; 
  report << "Events branches: " << events->GetListOfBranches()->GetEntries() << "\n"; 

  PrintBranchTable(events, report); 
  PrintBranchTable(events, std::cout); 

  const bool isSimulation = HasBranch(events, "genWeight"); 
  report << "\nSample classification: "
         << (isSimulation ? "Monte Carlo simulation" : "collision data") << "\n"; 

  if (isSimulation) { // solo se è simulazione
    report << std::setprecision(15);  // scrivi i numeri con alta precisione
    report << "Runs/genEventCount sum: " << SumLeaf(runs, "genEventCount") << "\n"; 
    report << "Runs/genEventSumw sum: " << SumLeaf(runs, "genEventSumw") << "\n"; 
    report << "Runs/genEventSumw2 sum: " << SumLeaf(runs, "genEventSumw2") << "\n"; 
  }

  const Long64_t entriesToRead = 
      std::min(std::max<Long64_t>(maxEvents, 1), events->GetEntries()); 
  report << "Events used for diagnostic plots: " << entriesToRead << "\n";

  const std::string rootOutputPath =
      outputDirectory + "/" + label + "_inspection.root"; 
  TFile outputFile(rootOutputPath.c_str(), "RECREATE"); 

  // Make histograms
  auto hNMuon = std::make_unique<TH1D>("h_nMuon", "Muon multiplicity", 6, -0.5, 5.5); 
  auto hMuonPt = std::make_unique<TH1D>("h_muonPt", "Muon transverse momentum", 50, 0., 200.); 
  auto hNJet = std::make_unique<TH1D>("h_nJet", "Jet multiplicity", 16, -0.5, 15.5); 
  auto hMet = std::make_unique<TH1D>("h_met", "Missing transverse momentum", 50, 0., 250.); 
  auto hLheNjets = std::make_unique<TH1D>("h_lheNjets", "LHE parton multiplicity", 8, -0.5, 7.5); 
  auto hGenWeight = std::make_unique<TH1D>("h_genWeight", "Generator-weight sign", 3, -1.5, 1.5); 

  hNMuon->SetDirectory(nullptr); 
  hMuonPt->SetDirectory(nullptr); 
  hNJet->SetDirectory(nullptr);
  hMet->SetDirectory(nullptr);
  hLheNjets->SetDirectory(nullptr);
  hGenWeight->SetDirectory(nullptr);

  StyleHistogram(hNMuon.get(), "N_{mu}"); 
  StyleHistogram(hMuonPt.get(), "Muon p_{T} [GeV]");
  StyleHistogram(hNJet.get(), "N_{jets}");
  StyleHistogram(hMet.get(), "p_{T}^{miss} [GeV]");
  StyleHistogram(hLheNjets.get(), "LHE_Njets");
  StyleHistogram(hGenWeight.get(), "sign(genWeight)");

  FillFromTree(events, hNMuon.get(), "nMuon", entriesToRead); 
  FillFromTree(events, hMuonPt.get(), "Muon_pt", entriesToRead); 
  FillFromTree(events, hNJet.get(), "nJet", entriesToRead);
  FillFromTree(events, hMet.get(), "MET_pt", entriesToRead);
  if (HasBranch(events, "LHE_Njets"))  
    FillFromTree(events, hLheNjets.get(), "LHE_Njets", entriesToRead);
  if (HasBranch(events, "genWeight"))
    FillFromTree(events, hGenWeight.get(),
                 "(genWeight>=0)-(genWeight<0)", entriesToRead); 

  TCanvas canvas("canvas", "NanoAOD inspection", 1200, 800); 
  canvas.Divide(3, 2); 
  for (int pad = 1; pad <= 6; ++pad) { 
    canvas.cd(pad);
    gPad->SetLeftMargin(0.14); 
    gPad->SetBottomMargin(0.13);
    gPad->SetTopMargin(0.16);
  }
  canvas.cd(1); hNMuon->Draw("HIST"); 
  canvas.cd(2); hMuonPt->Draw("HIST");
  canvas.cd(3); hNJet->Draw("HIST");
  canvas.cd(4); hMet->Draw("HIST");
  canvas.cd(5); hLheNjets->Draw("HIST");
  canvas.cd(6); hGenWeight->Draw("HIST");

  const std::string pngOutputPath =
      outputDirectory + "/" + label + "_inspection.png"; 
  canvas.SaveAs(pngOutputPath.c_str()); 
 
  outputFile.cd(); 
  hNMuon->Write(); 
  hMuonPt->Write();
  hNJet->Write();
  hMet->Write();
  hLheNjets->Write();
  hGenWeight->Write();
  outputFile.Close(); 
  report.close(); 

  std::cout << "Read " << entriesToRead << " events.\n"; 
  std::cout << "Report: " << reportPath << "\n"; 
  std::cout << "Histograms: " << rootOutputPath << "\n";
  std::cout << "Figure: " << pngOutputPath << "\n";
  return 0; // 
}
