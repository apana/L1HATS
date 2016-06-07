#include "TFile.h"
#include "TH1.h"
#include "TLegend.h"
#include "TString.h"
#include "TCanvas.h"
#include "TEfficiency.h"
#include "TStyle.h"

void plotEff (const TString& infile = "Eff.root") {

  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);

  TH1::SetDefaultSumw2();
  
  TFile* f_in = new TFile(infile);

  // ---- efficiency 
  
  TH1F* h_jet_denom = (TH1F*) f_in->Get("RecoJetPt");
  TH1F* h_jet_num = (TH1F*) f_in->Get("RecoJetPtTrg");
  
  TCanvas* cEff = new TCanvas("cEff","cEff");
  cEff->SetGrid(1,1);
  cEff->cd();

  TH2F* h_jet_axis = new TH2F("h_jet_axis",";p_{T} [GeV];Efficiency",64,-0.5,512.,20,0,1);
  h_jet_axis->GetYaxis()->SetTitleOffset(0.98);
  h_jet_axis->Draw();
  
  TEfficiency* h_jet_eff = new TEfficiency(*h_jet_num, *h_jet_denom);
  h_jet_eff->SetLineColor(kRed);
  h_jet_eff->SetMarkerColor(kRed);
  
  h_jet_eff->Draw("pe same");
  
 
  return;
}
