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

  RooRealVar *mass_2317 = new RooRealVar("mass_2317", "", r0_2317, r1_2317);
  RooRealVar *mass_2460 = new RooRealVar("mass_2460", "", r0_2460, r1_2460);

  TChain chain1("h1");
  TChain chain2("h1");
  chain1.Add("../Ds2317.root");
  chain2.Add("../Ds2460.root");

  Float_t
    m_ds17_1, cn_dsii_1, cn_dsi_1, mctr_1, hel_phi_1, hel_ds17_1, pst_ds17_1, c2_ds17_1, dgf_ds17_1, c2_dsii_1, dgf_dsii_1, c2_pi0_1, dgf_pi0_1, pst_dsii_1, m_dsii_1,
    cn_dsn_2, cn_dsw_2, m_ds60_2, mc_ds60_2, hel_phi_2, pstg_dss_2, hel_ds60_2, pst_ds60_2, c2_ds60_2, dgf_ds60_2, c2_dsn_2, dgf_dsn_2, c2_pi0_2, dgf_pi0_2, c2_dsst_2, dgf_dsst_2, pst_dsn_2, m_dsst_2;

  chain1.SetBranchAddress("m_ds17", &m_ds17_1);
  chain1.SetBranchAddress("m_dsii", &m_dsii_1);
  chain1.SetBranchAddress("cn_dsi", &cn_dsi_1);
  chain1.SetBranchAddress("cn_dsii", &cn_dsii_1);
  chain1.SetBranchAddress("mctr", &mctr_1);
  chain1.SetBranchAddress("hel_phi", &hel_phi_1);
  chain1.SetBranchAddress("hel_ds17", &hel_ds17_1);
  chain1.SetBranchAddress("pst_ds17", &pst_ds17_1);
  chain1.SetBranchAddress("c2_ds17", &c2_ds17_1);
  chain1.SetBranchAddress("dgf_ds17", &dgf_ds17_1);
  chain1.SetBranchAddress("c2_dsii", &c2_dsii_1);
  chain1.SetBranchAddress("dgf_dsii", &dgf_dsii_1);
  chain1.SetBranchAddress("c2_pi0", &c2_pi0_1);
  chain1.SetBranchAddress("dgf_pi0", &dgf_pi0_1);
  chain1.SetBranchAddress("pst_dsii", &pst_dsii_1);

  chain2.SetBranchAddress("cn_dsn", &cn_dsn_2);
  chain2.SetBranchAddress("cn_dsw", &cn_dsw_2);
  chain2.SetBranchAddress("m_dsst", &m_dsst_2);
  chain2.SetBranchAddress("m_ds60", &m_ds60_2);
  chain2.SetBranchAddress("mc_ds60", &mc_ds60_2);
  chain2.SetBranchAddress("pstg_dss", &pstg_dss_2);
  chain2.SetBranchAddress("hel_ds60", &hel_ds60_2);
  chain2.SetBranchAddress("pst_ds60", &pst_ds60_2);
  chain2.SetBranchAddress("c2_ds60", &c2_ds60_2);
  chain2.SetBranchAddress("dgf_ds60", &dgf_ds60_2);
  chain2.SetBranchAddress("c2_dsn", &c2_dsn_2);
  chain2.SetBranchAddress("dgf_dsn", &dgf_dsn_2);
  chain2.SetBranchAddress("c2_dsst", &c2_dsst_2);
  chain2.SetBranchAddress("dgf_dsst", &dgf_dsst_2);
  chain2.SetBranchAddress("pst_dsn", &pst_dsn_2);
  chain2.SetBranchAddress("hel_phi", &hel_phi_2);
  chain2.SetBranchAddress("c2_pi0", &c2_pi0_2);
  chain2.SetBranchAddress("dgf_pi0", &dgf_pi0_2);


  // TFile *MyDs2460 = new TFile("Ds2460_selected.root","RECREATE");

  TNtuple *nt2317 = new TNtuple("nt2317","Ds2317 mass","m_2317");
  TNtuple *nt2460 = new TNtuple("nt2460","Ds2460 mass","m_2460");

  // for (int i=0; i<chain2.GetEntries(); i++){
  for (int i=0; i<1E4; i++){
    chain1.GetEntry(i);
    Double_t m_var = m_ds17_1-m_dsii_1;
    if(m_var>r0_2317 && m_var<r1_2317 && TMath::Abs(hel_phi_1)>0.44 && cn_dsii_1!=2 && pst_ds17_1>2.79 && TMath::Prob(c2_ds17_1, dgf_ds17_1)>0.001 && TMath::Prob(c2_dsii_1, dgf_dsii_1)>0.001 && TMath::Prob(c2_pi0_1, dgf_pi0_1)>0.01){
      nt2317->Fill(m_var);
    }}

  for (int i=0; i<1E4;i++){
    chain2.GetEntry(i);
    Double_t m_var = m_ds60_2-m_dsst_2;
    if(m_var>r0_2460 && m_var<r1_2460 && TMath::Abs(hel_phi_2)>0.35 && cn_dsn_2!=2 && pst_ds60_2>3.5 && TMath::Prob(c2_ds60_2, dgf_ds60_2)>0.001 && TMath::Prob(c2_dsn_2, dgf_dsn_2)>0.001 && TMath::Prob(c2_pi0_2, dgf_pi0_2)>0.01 && TMath::Prob(c2_dsst_2, dgf_dsst_2)>0.001){
      nt2460->Fill(m_var);
      }}
  TFile *MyDs2317 = new TFile("../Ds2317_selected.root","RECREATE");
  nt2317->Write();
  MyDs2317->Close();
  // nt2460->Write();
  // MyDs2460->Close();

}
