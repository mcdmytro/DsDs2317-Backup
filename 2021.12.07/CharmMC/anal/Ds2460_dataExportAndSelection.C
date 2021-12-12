#include "RooRealVar.h"
#include "RooDataSet.h"

using namespace RooFit;

void Ds2460_dataExportAndSelection(){
  gROOT->Reset();

  Float_t
    r0_2317=0.25,
    r0_2460=0.25,
    r1_2317=0.45,
    r1_2460=0.45,
    x0_2317=0.25,
    x0_2460=0.25,
    x1_2317=0.45,
    x1_2460=0.45;

  RooRealVar *mass = new RooRealVar("mass", "", r0_2317, r1_2317);

  TChain chain("h1");
  chain.Add("../Ds2460.root");

  Float_t
    cn_dsn, cn_dsw, m_ds60, mc_ds60, hel_phi, pstg_dss, hel_ds60, pst_ds60, c2_ds60, dgf_ds60, c2_dsn, dgf_dsn, c2_pi0, dgf_pi0, c2_dsst, dgf_dsst, pst_dsn, m_dsst;

  chain.SetBranchAddress("cn_dsn", &cn_dsn);
  chain.SetBranchAddress("cn_dsw", &cn_dsw);
  chain.SetBranchAddress("m_dsst", &m_dsst);
  chain.SetBranchAddress("m_ds60", &m_ds60);
  chain.SetBranchAddress("mc_ds60", &mc_ds60);
  chain.SetBranchAddress("pstg_dss", &pstg_dss);
  chain.SetBranchAddress("hel_ds60", &hel_ds60);
  chain.SetBranchAddress("pst_ds60", &pst_ds60);
  chain.SetBranchAddress("c2_ds60", &c2_ds60);
  chain.SetBranchAddress("dgf_ds60", &dgf_ds60);
  chain.SetBranchAddress("c2_dsn", &c2_dsn);
  chain.SetBranchAddress("dgf_dsn", &dgf_dsn);
  chain.SetBranchAddress("c2_dsst", &c2_dsst);
  chain.SetBranchAddress("dgf_dsst", &dgf_dsst);
  chain.SetBranchAddress("pst_dsn", &pst_dsn);
  chain.SetBranchAddress("hel_phi", &hel_phi);
  chain.SetBranchAddress("c2_pi0", &c2_pi0);
  chain.SetBranchAddress("dgf_pi0", &dgf_pi0);

  TNtuple *h1 = new TNtuple("h1","Ds2460 mass delta","dm_2460");

  for (int i=0; i<chain.GetEntries(); i++){
    chain.GetEntry(i);
    Double_t m_var = m_ds60-m_dsst;
    if(m_var>r0_2460 && m_var<r1_2460 && TMath::Abs(hel_phi)>0.35 && cn_dsn!=2 && pst_ds60>3.5 && TMath::Prob(c2_ds60, dgf_ds60)>0.001 && TMath::Prob(c2_dsn, dgf_dsn)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && TMath::Prob(c2_dsst, dgf_dsst)>0.001){
      h1->Fill(m_var);
      }}

  TFile *MyFile = new TFile("../Ds2460_selected.root","RECREATE");
  h1->Write();
  MyFile->Close();
}
