/*  Simple C plus program */
#include <iostream>
#include <sstream>
#include <string>
#include <math.h>
#include <vector>

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include <boost/program_options.hpp>
#include <boost/filesystem.hpp>

#include "L1TriggerDPG/L1Menu/macros/L1Ntuple.h"

#include "DataFormats/L1Trigger/interface/EtSum.h"

#include "L1Trigger/L1TNtuples/interface/L1AnalysisEventDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisL1UpgradeDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoJetDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoElectronDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMuon2DataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoMetFilterDataFormat.h"
#include "L1Trigger/L1TNtuples/interface/L1AnalysisRecoTauDataFormat.h"

enum EtSumType {
  ETT = l1t::EtSum::EtSumType::kTotalEt,
  HTT = l1t::EtSum::EtSumType::kTotalHt,
  ETM = l1t::EtSum::EtSumType::kMissingEt,
  HTM = l1t::EtSum::EtSumType::kMissingHt,
}; // Base on "DataFormats/L1Trigger/interface/EtSum.h"

class RateExample: public L1Ntuple
{
public:
  // constructor
  RateExample(std::string infiles, std::string outfile, float JetThresh){

    SelBx=0; // look at bunch crossing 0
    Nevts=0;
    jetThreshold=JetThresh;

    std::cout << "Writing output to: " << outfile << std::endl;
    outrootfile =  new TFile( outfile.c_str(), "RECREATE");

    SelectTree(true);
    if (infiles.find(".root") != std::string::npos) {
      std::cout << "Reading RootFile: " << infiles << std::endl;
      Open(infiles);
    }else{
      std::cout << "Reading Filelist: " << infiles << std::endl;
      if (! OpenWithList(infiles)) exit(0);
    }
  }
  ~RateExample() {
    writeHistograms();
    outrootfile->Close();
    std::cout << "Program execution terminated" << std::endl;
  }
  void loop(int);
  void bookHistograms();
  void writeHistograms();
  void scaleHistograms(int);
  void integrateRate();

  int GetSumEtIdx( EtSumType );
  double SumETVal( EtSumType );
  bool checkPFJetID(int);

private:
  int SelBx, Nevts;
  double jetThreshold;
  std::map<std::string,TH1F*> hTH1Fs,hRates,hEffs;
  TFile        *outrootfile;
};

int RateExample::GetSumEtIdx(EtSumType type)
{
  for (unsigned int i = 0; i < upgrade_->sumType.size(); ++i )
  {
    if (upgrade_->sumType.at(i) == type && upgrade_->sumBx.at(i) == SelBx)
      return i;
  }
  return -1;
}       // -----  end of function L1AlgoFactory::GetSumEtIdx  -----


bool RateExample::checkPFJetID(int ijet){
  bool jid=false;
  double chf   = recoJet_->chef.at(ijet);
  double nhf   = recoJet_->nhef.at(ijet);
  double phf   = recoJet_->pef.at(ijet);
  double elf   = recoJet_->eef.at(ijet);
  double chm   = recoJet_->chMult.at(ijet);
  int npr   = recoJet_->cMult.at(ijet) + recoJet_->nMult.at(ijet);

  jid  = (npr>1 && phf<0.99 && nhf<0.99 && ((fabs(recoJet_->eta.at(ijet))<=2.4 && elf<0.99 && chf>0 && chm>0) || fabs(recoJet_->eta.at(ijet))>2.4)) ;
  return jid;

}

double RateExample::SumETVal( EtSumType type ) {

  Float_t TheVal = -10;
  int idx = GetSumEtIdx(type);
  assert(upgrade_->sumType.at(idx) == type);
  if(upgrade_->sumBx.at(idx)==SelBx) TheVal =upgrade_->sumEt.at(idx);

  return TheVal;
}

void RateExample::bookHistograms(){

  // generic histograms
  hTH1Fs["NEV"] = new TH1F("NEV","NEV; Number of events processed; ",1,-.5,.5);

  // rate histograms
  hRates["HTT"] = new TH1F("HTT","HTT; HTT (GeV); Integrated Rate [kHz]",512,-.5,511.5);
  hRates["HTM"] = new TH1F("HTM","HTM; HTM (GeV); Integrated Rate [kHz]",150,-.5,149.5);
  hRates["ETT"] = new TH1F("ETT","ETT; ETT (GeV); Integrated Rate [kHz]",150,-.5,149.5);
  hRates["ETM"] = new TH1F("ETM","ETM; ETM (GeV); Integrated Rate [kHz]",512,-.5,511.5);

  hRates["JET"] = new TH1F("JET","JET; JET (GeV); Integrated Rate [kHz]",512,-.5,511.5);

  hRates["EG"] = new TH1F("EG","EG; EG (GeV); Integrated Rate [kHz]",100,-.5,99.5);
  hRates["IsoEG"] = new TH1F("IsoEG","IsoEG; IsoEG (GeV); Integrated Rate [kHz]",100,-.5,99.5);

  // Efficiency histograms
  hEffs["RecoJetPt"] = new TH1F("RecoJetPt","JetPt; JetPt (GeV) ",512,-.5,511.5);
  hEffs["RecoJetPtTrg"] = new TH1F("RecoJetPtTrg","JetPtTrg; JetPt (GeV) ",512,-.5,511.5);
}

void RateExample::scaleHistograms(int nBunches)
{

  double scale = 11246.; // ZB per bunch in Hz

  scale /= Nevts*1000.; // in kHz
  scale *= nBunches;

  for(auto h : hRates)
  {
    h.second->Scale(scale);
  }

  return;
}       // -----  end of function L1Menu2016::WriteHistogram  -----

void RateExample::integrateRate()
{

  for(auto h : hRates)
  {
    int nbins=h.second->GetNbinsX();
    for (int ibin=0; ibin<nbins; ++ibin){
      int ibin1=ibin+1;
      int ibin2=nbins;
      double err;
      double fint=h.second->IntegralAndError(ibin1,ibin2,err);

      h.second->SetBinContent(ibin+1,fint);
      //h.second->SetBinError(ibin+1,err);  // error bars looking weird in integrated plot
      h.second->SetBinError(ibin+1,0.);
    }
  }

  return;
}       // -----  end of function L1Menu2016::WriteHistogram  -----

void RateExample::writeHistograms()
{

  outrootfile->cd();
  for(auto h : hTH1Fs)
  {
    h.second->Write();
  }
  for(auto h : hRates)
  {
    h.second->Write();
  }

  for(auto h : hEffs)
  {
    h.second->Write();
  }

  return;
}       // -----  end of function L1Menu2016::WriteHistogram  -----

void RateExample::loop(int maxEvents){

  int i=0;
  while(true)
  {
    Long64_t ientry = LoadTree(i);
    if (ientry < 0) break;

    GetEntry(i); i++;
    if (maxEvents>0 && i>=maxEvents) break;


    if (i<0)  //first n events
      {
    std::cout << "--------------------- Event "<<i<<" ---------------------"<<std::endl;
    //event_
    std::cout << "L1Tree         : event_->run =\t"<<event_->run<< std::endl;
    std::cout << "L1Tree         : event_->event =\t"<<event_->event<<std::endl;

    std::cout << "\nL1Tree:        : upgrade_->nJets=\t" <<  upgrade_->nJets << std::endl;
    }
    if (i % 200000 == 0)
      std::cout << "Processed " << i << " events." << std::endl;

    hTH1Fs["NEV"]->Fill(0.);

    double TheETT(-10),TheETM(-10),TheHTT(-10),TheHTM(-10);
    TheETT =SumETVal( EtSumType::ETT );
    TheETM =SumETVal( EtSumType::ETM );
    TheHTT =SumETVal( EtSumType::HTT );
    TheHTM =SumETVal( EtSumType::HTM );

    if (TheETT>0.5) hRates["ETT"]->Fill(TheETT);
    if (TheETM>0.5) hRates["ETM"]->Fill(TheETM);
    if (TheHTT>0.5) hRates["HTT"]->Fill(TheHTT);
    if (TheHTM>0.5) hRates["HTM"]->Fill(TheHTM);


    //EG rates
    double EGEt(-10.), isoEGEt(-10.);
    for (UInt_t ue=0; ue < upgrade_->nEGs; ue++){
      Int_t bx = upgrade_->egBx.at(ue);
      if(bx != 0) continue;
      Float_t pt  = upgrade_->egEt.at(ue);
      if (pt>EGEt)EGEt=pt;
      if (upgrade_->egIso.at(ue)){
	if (pt>isoEGEt)isoEGEt=pt;
      }
    }  // end loop over EM objects
    if (EGEt>0.5) hRates["EG"]->Fill(EGEt);
    if (isoEGEt>0.5) hRates["IsoEG"]->Fill(isoEGEt);

    //Jet rates
    double JetEt(-10.);
    for (UInt_t ue=0; ue < upgrade_->nJets; ue++){
      Int_t bx = upgrade_->jetBx.at(ue);
      if(bx != 0) continue;
      Float_t pt  = upgrade_->jetEt.at(ue);
      if (pt>JetEt)JetEt=pt;
    }  // end loop over EM objects
    if (JetEt>0.5) hRates["JET"]->Fill(JetEt);

    // Now simple Jet Efficiency
    if (doRecoJet && recoJet_->nJets > 0 ){
      // check leading jet et
      bool goodJet=checkPFJetID(0);
      if (goodJet){
	double ptL=recoJet_->etCorr.at(0);
	hEffs["RecoJetPt"]->Fill(ptL);
	//std::cout << "recoJet:" << recoJet_->etCorr.at(0) << std::endl;
	if (JetEt>jetThreshold){
	  hEffs["RecoJetPtTrg"]->Fill(ptL);
	}
      }
    }
  }

  Nevts=i;
  std::cout << "\nDone with event loop. Processed: " << Nevts << " events" << std::endl;
}
using std::vector;
int main(int argc, char *argv[])
{

  namespace po = boost::program_options;
  // Declare the supported options.
  boost::program_options::options_description desc("Allowed options");
  const std::string defaultntuple = "r259721_tsgv3.list";
  const std::string defaultoutput = "rates.root";

  desc.add_options()
    ("help,h", "produce help message")
    ("filelist,l",    po::value<std::string>()->default_value(defaultntuple), "set the input ntuple list")
    ("outfilename,o", po::value<std::string>()->default_value(defaultoutput), "set output file name")
    ("maxEvent,n",    po::value<int>()->default_value(-1),                    "run number of events; -1 for all")
    ("nBunches,b",    po::value<int>()->default_value(1),                     "set number of bunches")
    ("jetThreshold",   po::value<float>()->default_value(120),                     "jet threshold for efficiency")
    ;


  po::variables_map vm;
  po::store(po::parse_command_line(argc, argv, desc), vm);
  po::notify(vm);

  if (vm.count("help")) {
    std::cout << desc << std::endl;
    return 1;
  }

  int nevts   =vm["maxEvent"].as<int>();
  int nBunches=vm["nBunches"].as<int>();
  double jetThreshold=vm["jetThreshold"].as<float>();
  std::string inFiles=vm["filelist"].as<std::string>();
  std::string outFile=vm["outfilename"].as<std::string>();

  RateExample myRates(inFiles,outFile,jetThreshold);

  myRates.bookHistograms();
  myRates.loop(nevts);
  myRates.scaleHistograms(nBunches);
  myRates.integrateRate();


}
