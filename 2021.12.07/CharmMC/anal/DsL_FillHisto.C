#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"
#include "RooDataHist.h"

using namespace RooFit;

void DsL_FillHisto(){
  gROOT->Reset();
  gSystem->Load("libRooFit");

  Float_t
      r0_2317=0.10,
      r1_2317=0.60,
      r0_2460=0.10,
      r1_2460=0.60;

  Int_t
    Nbins=100;

  TChain chain("h1");
  chain.Add("../Ds2317.root");
  chain.Add("../Ds2460.root");

  TFile *f = new TFile("./root_files/DsJ.root","recreate");
  TH1D *ds17_del = new TH1D("ds17_del","Ds2317 mass delta", Nbins, r0_2317, r1_2317);
  TH1D *ds60_del = new TH1D("ds60_del","Ds2460 mass delta", Nbins, r0_2460, r1_2460);

  Float_t
    m_ds17, cn_dsii, cn_dsi, mctr, hel_phi, hel_ds17, pst_ds17, c2_ds17, dgf_ds17, c2_dsii, dgf_dsii, c2_pi0, dgf_pi0, pst_dsii, m_dsii,
    cn_dsn, cn_dsw, m_ds60, mc_ds60, pstg_dss, hel_ds60, pst_ds60, c2_ds60, dgf_ds60, c2_dsn, dgf_dsn, c2_dsst, dgf_dsst, pst_dsn, m_dsst;

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
  chain.SetBranchAddress("hel_ds17", &hel_ds17);

  chain.SetBranchAddress("cn_dsn", &cn_dsn);
  chain.SetBranchAddress("cn_dsw", &cn_dsw);
  chain.SetBranchAddress("m_dsst", &m_dsst);
  chain.SetBranchAddress("m_ds60", &m_ds60);
  chain.SetBranchAddress("mc_ds60", &mc_ds60);
  chain.SetBranchAddress("pstg_dss", &pstg_dss);
  chain.SetBranchAddress("hel_ds60", &hel_ds60);
  chain.SetBranchAddress("pst_ds60", &pst_ds60);
  chain.SetBranchAddress("pstg_dss", &pstg_dss);
  chain.SetBranchAddress("c2_ds60", &c2_ds60);
  chain.SetBranchAddress("dgf_ds60", &dgf_ds60);
  chain.SetBranchAddress("c2_dsn", &c2_dsn);
  chain.SetBranchAddress("dgf_dsn", &dgf_dsn);
  chain.SetBranchAddress("c2_dsst", &c2_dsst);
  chain.SetBranchAddress("dgf_dsst", &dgf_dsst);
  chain.SetBranchAddress("pst_dsn", &pst_dsn);
  chain.SetBranchAddress("hel_ds60", &hel_ds60);
  chain.SetBranchAddress("pstg_dss", &pstg_dss);

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if((m_ds17-m_dsii)>r0_2317 && (m_ds17-m_dsii)<r1_2317 && TMath::Abs(hel_phi)>0.35 && cn_dsii!=2 && pst_ds17>3.5 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01)
         ds17_del->Fill(m_ds17-m_dsii);
    }

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if((m_ds60-m_dsst)>r0_2460 && (m_ds60-m_dsst)<r1_2460 && TMath::Abs(hel_phi)>0.35 && cn_dsn!=2 && pst_ds60>3.5 && TMath::Prob(c2_ds60, dgf_ds60)>0.001 && TMath::Prob(c2_dsn, dgf_dsn)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && TMath::Prob(c2_dsst, dgf_dsst)>0.001)
       ds60_del->Fill(m_ds60-m_dsst);
      }

    ds17_del->Write();
    ds60_del->Write();
    f->Close();
}
