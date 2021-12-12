#include "RooGlobalFunc.h"
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

void Dsp_subtracted_m5w2s0(){
gROOT->Reset();
gSystem->Load("libRooFit");
Float_t range0=4., range1=8.;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
chain.Add("../Ds2317inDs2317_m5w2.root");

Float_t m_dsp, hel_phi, hel_ds17, pst_ds17, m_ds17, m_dsi, dgf_dsii, c2_dsii, cn_dsii, dgf_pi0, c2_pi0;

chain.SetBranchAddress("m_dsp", &m_dsp);
chain.SetBranchAddress("m_dsi", &m_dsi);
chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
chain.SetBranchAddress("c2_dsii", &c2_dsii);
chain.SetBranchAddress("cn_dsii", &cn_dsii);
chain.SetBranchAddress("m_ds17", &m_ds17);
chain.SetBranchAddress("pst_ds17", &pst_ds17);
chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
chain.SetBranchAddress("c2_pi0", &c2_pi0);
chain.SetBranchAddress("hel_phi", &hel_phi);
chain.SetBranchAddress("hel_ds17", &hel_ds17);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  float m_var=m_dsp-m_dsi-m_ds17+1.9683+2.3178;
  chain.GetEntry(i);
  if(m_var>range0 && m_var<range1 && TMath::Abs(hel_phi)>0.35 && hel_ds17<0.7 && cn_dsii!=2){
    mass->setVal(m_var);
    data->add(RooArgSet(*mass));
  }}

 TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);
 RooPlot *frame = mass->frame();

 data->plotOn(frame);

 gPad->SetLeftMargin(0.2);
 gPad->SetBottomMargin(0.12) ;
 frame->SetTitle("");
 frame->GetXaxis()->SetTitle("M(D_{s}D_{s0}*(2317)) [GeV/c^{2}] ");
 frame->GetYaxis()->SetTitle("Entries / 40 MeV/c^{2}");
 frame->GetXaxis()->SetTitleSize(0.05);
 frame->GetYaxis()->SetTitleSize(0.05);
 frame->GetXaxis()->SetLabelSize(0.05);
 frame->GetYaxis()->SetLabelSize(0.05);
 frame->GetXaxis()->SetNdivisions(505);
 frame->GetYaxis()->SetNdivisions(505);
 frame->GetYaxis()->SetTitleOffset(1.8);

 frame->Draw();

 cv->SaveAs("./pdf_files/SMC_DspinDs2317_subtracted_m5w2s0.pdf","pdf");
}
