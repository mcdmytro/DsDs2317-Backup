#include "RooRealVar.h"
#include "RooDataSet.h"

using namespace RooFit;

void Ds2317_dataExportAndSelection(){
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
  chain.Add("../Ds2317.root");

  Float_t
    m_ds17, cn_dsii, cn_dsi, mctr, hel_phi, hel_ds17, pst_ds17, c2_ds17, dgf_ds17, c2_dsii, dgf_dsii, c2_pi0, dgf_pi0, pst_dsii, m_dsii;

  chain.SetBranchAddress("m_ds17", &m_ds17);
  chain.SetBranchAddress("m_dsii", &m_dsii);
  chain.SetBranchAddress("cn_dsi", &cn_dsi);
  chain.SetBranchAddress("cn_dsii", &cn_dsii);
  chain.SetBranchAddress("mctr", &mctr);
  chain.SetBranchAddress("hel_phi", &hel_phi);
  chain.SetBranchAddress("hel_ds17", &hel_ds17);
  chain.SetBranchAddress("pst_ds17", &pst_ds17);
  chain.SetBranchAddress("c2_ds17", &c2_ds17);
  chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
  chain.SetBranchAddress("c2_dsii", &c2_dsii);
  chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
  chain.SetBranchAddress("c2_pi0", &c2_pi0);
  chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
  chain.SetBranchAddress("pst_dsii", &pst_dsii);

  TNtuple *h1 = new TNtuple("h1","Ds2317 mass delta","dm_2317");

  for (int i=0; i<chain.GetEntries(); i++){
    chain.GetEntry(i);
    Double_t m_var = m_ds17-m_dsii;
    if(m_var>r0_2317 && m_var<r1_2317 && TMath::Abs(hel_phi)>0.44 && cn_dsii!=2 && pst_ds17>2.79 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01){
      h1->Fill(m_var);
    }}

  TFile *MyFile = new TFile("../Ds2317_selected.root","RECREATE");
  h1->Write();
  MyFile->Close();
}
