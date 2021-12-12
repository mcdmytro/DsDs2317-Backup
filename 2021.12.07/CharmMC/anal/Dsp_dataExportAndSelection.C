#include "RooRealVar.h"
#include "RooDataSet.h"

using namespace RooFit;

void Dsp_dataExportAndSelection(){
  gROOT->Reset();

  Float_t range0=3.5, range1=10;
  Float_t x0_ds17=2.305, x1_ds17=2.329;

  RooRealVar *mass = new RooRealVar("mass", "", range0, range1);

  TChain chain("h1");
  chain.Add("../Ds2317.root");

  Float_t m_dsp,   hel_phi,   pst_ds17,   m_ds17,   m_dsii,   dgf_dsii,   c2_dsii,   cn_dsii,   dgf_pi0,   c2_pi0;

  chain.SetBranchAddress("m_dsp", &m_dsp);
  chain.SetBranchAddress("m_dsii", &m_dsii);
  chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
  chain.SetBranchAddress("c2_dsii", &c2_dsii);
  chain.SetBranchAddress("cn_dsii", &cn_dsii);
  chain.SetBranchAddress("m_ds17", &m_ds17);
  chain.SetBranchAddress("pst_ds17", &pst_ds17);
  chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
  chain.SetBranchAddress("c2_pi0", &c2_pi0);
  chain.SetBranchAddress("hel_phi", &hel_phi);

  TNtuple *h1 = new TNtuple("h1","Dsp mass","dsp");

  for (int i=0; i<chain.GetEntries(); i++){
    chain.GetEntry(i);
    float m_var=m_dsp-m_dsii-m_ds17+1.9683+2.3178;
    // if(m_var>range0 && m_var<range1 && m_ds17>x0_ds17 && m_ds17<x1_ds17 && TMath::Abs(hel_phi)>0.35 && cn_dsii!=2 && pst_ds17>3.50 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01){
    if(m_var>range0 && m_var<range1 && m_ds17>x0_ds17 && m_ds17<x1_ds17 && TMath::Abs(hel_phi)>0.44 && cn_dsii!=2 && pst_ds17>2.79 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01){
      h1->Fill(m_var);
    }}

  // TFile *MyFile = new TFile("../Dsp_tree_tough.root","RECREATE");
  TFile *MyFile = new TFile("../Dsp_histo_optimized.root","RECREATE");

  h1->Write();
  MyFile->Close();
}
