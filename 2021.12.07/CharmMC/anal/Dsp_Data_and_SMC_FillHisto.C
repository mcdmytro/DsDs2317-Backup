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

void Dsp_Data_and_SMC_FillHisto(){
  gROOT->Reset();
  gSystem->Load("libRooFit");

  Float_t range0=3.5, range1=10;
  Float_t x0_ds17=2.305, x1_ds17=2.329;

  Int_t
    Nbins=100;

  TChain chain("h1");
  TChain chain1("h1");
  chain.Add("../Ds2317.root");
  chain1.Add("../../SignalMC/Ds2317inDs2317_m4600w67s0.root");

  TFile *f = new TFile("../Dsp_Data_and_SMC.root","recreate");

  TH1D *dsp_data = new TH1D("dsp_data","Dsp on data", Nbins, range0, range1);
  TH1D *dsp_SMC = new TH1D("dsp_SMC","Dsp on SMC", Nbins, range0, range1);

  Float_t m_dsp,   hel_phi,   pst_ds17,   m_ds17,   m_dsii,   dgf_dsii,   c2_dsii,   cn_dsii,   dgf_pi0,   c2_pi0;
  Float_t m_dsp_1, hel_phi_1, pst_ds17_1, m_ds17_1, m_dsii_1, dgf_dsii_1, c2_dsii_1, cn_dsii_1, dgf_pi0_1, c2_pi0_1, hel_ds17_1;

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

  chain1.SetBranchAddress("m_dsp", &m_dsp_1);
  chain1.SetBranchAddress("m_dsii", &m_dsii_1);
  chain1.SetBranchAddress("dgf_dsii", &dgf_dsii_1);
  chain1.SetBranchAddress("c2_dsii", &c2_dsii_1);
  chain1.SetBranchAddress("cn_dsii", &cn_dsii_1);
  chain1.SetBranchAddress("m_ds17", &m_ds17_1);
  chain1.SetBranchAddress("pst_ds17", &pst_ds17_1);
  chain1.SetBranchAddress("dgf_pi0", &dgf_pi0_1);
  chain1.SetBranchAddress("c2_pi0", &c2_pi0_1);
  chain1.SetBranchAddress("hel_phi", &hel_phi_1);
  chain1.SetBranchAddress("hel_ds17", &hel_ds17_1);

  for (int i=0; i<chain.GetEntries();i++){
    float m_var=m_dsp-m_dsii-m_ds17+1.9683+2.3178;
    chain.GetEntry(i);
    if(m_var>range0 && m_var<range1 && m_ds17>x0_ds17 && m_ds17<x1_ds17 && pst_ds17>2.79 && TMath::Abs(hel_phi)>0.44 && cn_dsii!=2)
      dsp_data->Fill(m_var);
    }

  for (int i=0; i<chain1.GetEntries();i++){
    float m_var=m_dsp_1-m_dsii_1-m_ds17_1+1.9683+2.3178;
    chain1.GetEntry(i);
    if(m_var>range0 && m_var<range1 && TMath::Abs(hel_phi_1)>0.35 && hel_ds17_1<0.7 && cn_dsii_1!=2)
      dsp_SMC->Fill(m_var);
  }

  Double_t
    data_maxbin=dsp_data->GetBinContent(0),
    SMC_maxbin=dsp_SMC->GetBinContent(0);

  for(int i=0; i<Nbins; i++){
    if(dsp_data->GetBinContent(i)>data_maxbin)
      data_maxbin=dsp_data->GetBinContent(i);
    if(dsp_SMC->GetBinContent(i)>SMC_maxbin)
      SMC_maxbin=dsp_SMC->GetBinContent(i);
  }

    dsp_data->Scale(1.0);
    dsp_data->Write();
    dsp_SMC->Scale(data_maxbin/SMC_maxbin);
    dsp_SMC->Write();
    f->Close();
}
