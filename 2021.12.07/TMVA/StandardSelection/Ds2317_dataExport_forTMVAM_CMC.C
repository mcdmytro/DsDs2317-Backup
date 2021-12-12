#include "RooRealVar.h"
#include "RooDataSet.h"

using namespace RooFit;

void Ds2317_dataExport_forTMVAM_CMC(){
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
  chain.Add("../../../../DsDs2317-rootFiles/CharmMC/Ds2317_BCS.root");

  Float_t
    m_ds17, mc_ds17, cn_dsii, cn_dsi, mctr, hel_phi, hel_ds17, pst_ds17, c2_ds17, dgf_ds17, c2_dsii, dgf_dsii, c2_pi0, dgf_pi0, pst_dsii, m_dsii;

  chain.SetBranchAddress("m_ds17", &m_ds17);
  chain.SetBranchAddress("mc_ds17", &mc_ds17);
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

  TNtuple *TreeS = new TNtuple("TreeS", "Sgn tree", "dm_ds17:mc_ds17:pst_ds17:hel_phi:pc2_ds17:pc2_dsii:pc2_pi0");
  TNtuple *TreeB = new TNtuple("TreeB", "Bkg tree", "dm_ds17:mc_ds17:pst_ds17:hel_phi:pc2_ds17:pc2_dsii:pc2_pi0");
  TNtuple *Tree  = new TNtuple("Tree",  "All tree", "dm_ds17:mc_ds17:pst_ds17:hel_phi:pc2_ds17:pc2_dsii:pc2_pi0");

  for (int i=0; i<chain.GetEntries(); i++){
    chain.GetEntry(i);
    Double_t m_var = m_ds17-m_dsii;
    if(m_var>r0_2317 && m_var<r1_2317 && cn_dsii!=2){
      Tree->Fill(m_var, TMath::Abs(mc_ds17), pst_ds17, TMath::Abs(hel_phi), TMath::Prob(c2_ds17, dgf_ds17), TMath::Prob(c2_dsii, dgf_dsii), TMath::Prob(c2_pi0, dgf_pi0));
      if(TMath::Abs(mc_ds17)==10431)
        TreeS->Fill(m_var, TMath::Abs(mc_ds17), pst_ds17, TMath::Abs(hel_phi), TMath::Prob(c2_ds17, dgf_ds17), TMath::Prob(c2_dsii, dgf_dsii), TMath::Prob(c2_pi0, dgf_pi0));
      else
        TreeB->Fill(m_var, TMath::Abs(mc_ds17), pst_ds17, TMath::Abs(hel_phi), TMath::Prob(c2_ds17, dgf_ds17), TMath::Prob(c2_dsii, dgf_dsii), TMath::Prob(c2_pi0, dgf_pi0));
      }
    }

  TFile *MyFile = new TFile("CMC_Ds2317_forTMVA.root","RECREATE");
  TreeS->Write();
  TreeB->Write();
  Tree ->Write();

  MyFile->Close();
  return 0;
}
