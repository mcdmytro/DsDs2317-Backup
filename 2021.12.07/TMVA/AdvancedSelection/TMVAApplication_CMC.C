#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"

using namespace TMVA;

void TMVAApplication_CMC(){

  TMVA::Tools::Instance();
  TMVA::Reader *reader = new TMVA::Reader("!Color:!Silent");

  Float_t dm_ds17, mc_ds17, pst_ds17, c2_ds17, pc2_ds17, pst_dsii, c2_dsii, pc2_dsii, pc2_pi0, p_pi0, c2_pi0, egam1, egam2, hel_phi;

  reader->AddVariable("pst_ds17",  &pst_ds17);
  reader->AddVariable("c2_ds17",   &c2_ds17);
  reader->AddVariable("pc2_ds17",  &pc2_ds17);
  reader->AddVariable("pst_dsii",  &pst_dsii);
  reader->AddVariable("c2_dsii",   &c2_dsii);
  reader->AddVariable("pc2_dsii",  &pc2_dsii);
  reader->AddVariable("pc2_pi0",   &pc2_pi0);
  reader->AddVariable("p_pi0",     &p_pi0);
  reader->AddVariable("c2_pi0",    &c2_pi0);
  reader->AddVariable("egam1",     &egam1);
  reader->AddVariable("egam2",     &egam2);
  reader->AddVariable("hel_phi",   &hel_phi);

  reader->AddSpectator( "dm_ds17", &dm_ds17);
  reader->AddSpectator( "mc_ds17", &mc_ds17);

  reader->BookMVA("MLP classifier", "dataset/weights/MVAnalysis_MLP.weights.xml");

  TChain* chain =  new TChain("Tree");
  chain->Add("CMC_Ds2317_forTMVA.root");

  chain->SetBranchAddress("dm_ds17",  &dm_ds17);
  chain->SetBranchAddress("mc_ds17",  &mc_ds17);
  chain->SetBranchAddress("pst_ds17", &pst_ds17);
  chain->SetBranchAddress("c2_ds17",  &c2_ds17);
  chain->SetBranchAddress("pc2_ds17", &pc2_ds17);
  chain->SetBranchAddress("pst_dsii", &pst_dsii);
  chain->SetBranchAddress("c2_dsii",  &c2_dsii);
  chain->SetBranchAddress("pc2_dsii", &pc2_dsii);
  chain->SetBranchAddress("pc2_pi0",  &pc2_pi0);
  chain->SetBranchAddress("p_pi0",    &p_pi0);
  chain->SetBranchAddress("c2_pi0",   &c2_pi0);
  chain->SetBranchAddress("egam1",    &egam1);
  chain->SetBranchAddress("egam2",    &egam2);
  chain->SetBranchAddress("hel_phi",  &hel_phi);

  TNtuple* out  = new TNtuple("out",  "Out tree", "mvaValue:dm_ds17:mc_ds17");

  for(Long64_t ievt=0; ievt<chain->GetEntries(); ievt++){
      chain->GetEntry(ievt);
      Double_t mvaValue = reader->EvaluateMVA("MLP classifier");
      out->Fill(mvaValue, dm_ds17, mc_ds17);
    }

  TFile *MyFile = new TFile("CMC_outTMVA.root","RECREATE");
  out ->Write();
  MyFile->Close();

  delete reader;

  cout << '\n' << ">>>>> Not bad, not bad ... What's next? <<<<<" << '\n';
  return 0;
}
