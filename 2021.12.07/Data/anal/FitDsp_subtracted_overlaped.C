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

void FitDsp_subtracted_overlaped(){
gROOT->Reset();
gSystem->Load("libRooFit");
Float_t range0=3.5, range1=10;
Float_t x0=2.215, x1=2.35;
Float_t x0_ds17=2.305, x1_ds17=2.329;
Float_t width=0.008;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
TChain chain1("h1");
chain.Add("../Ds2317.root");
chain1.Add("../../SignalMC/Ds2317inDs2317_m4600w67s0.root");
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

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");
RooDataSet *data1 = new RooDataSet("data1","", RooArgSet(*mass),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  float m_var=m_dsp-m_dsii-m_ds17+1.9683+2.3178;
  chain.GetEntry(i);
  if(m_var>range0 && m_var<range1 && m_ds17>x0_ds17 && m_ds17<x1_ds17 && pst_ds17>2.79 && TMath::Abs(hel_phi)>0.44 && cn_dsii!=2){
    mass->setVal(m_var);
    data->add(RooArgSet(*mass));
  }}

for (int i=0; i<chain1.GetEntries();i++){
  float m_var=m_dsp_1-m_dsii_1-m_ds17_1+1.9683+2.3178;
  chain1.GetEntry(i);
  if(m_var>range0 && m_var<range1 && TMath::Abs(hel_phi_1)>0.35 && hel_ds17_1<0.7 && cn_dsii_1!=2){
  mass->setVal(m_var);
  data1->add(RooArgSet(*mass));
}}

 TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);
 RooPlot *frame = mass->frame();
 data->plotOn(frame);
 data1->plotOn(frame, Rescale(5E-3), MarkerColor(kRed));

 gPad->SetLeftMargin(0.2);
 gPad->SetBottomMargin(0.12) ;
 frame->SetTitle("");
 frame->GetXaxis()->SetTitle("M(D_{s}D_{s0}^{*}(2317)), GeV ");
 //frame->GetYaxis()->SetTitle("Entries");
 frame->GetXaxis()->SetTitleSize(0.05);
 frame->GetYaxis()->SetTitleSize(0.05);
 frame->GetXaxis()->SetLabelSize(0.05);
 frame->GetYaxis()->SetLabelSize(0.05);
 frame->GetXaxis()->SetNdivisions(505);
 frame->GetYaxis()->SetNdivisions(505);
 frame->GetYaxis()->SetTitleOffset(1.8);
 frame->SetMaximum(100);
 frame->Draw();

 cv->SaveAs("SMCandData_DspinDs2317_subtracted_overlapped.pdf","pdf");
}
