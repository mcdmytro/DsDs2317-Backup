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

void FitDsp_direct(){
gROOT->Reset();
gSystem->Load("libRooFit");
Float_t
  // range0=10.4, range1=10.7,
  range0=4.2, range1=4.4;
Float_t x0_ds17=2.305, x1_ds17=2.329;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
// chain.Add("Ds2317inDs2317_direct.root");
// chain.Add("Ds2317inDs2317_m5w2.root");
// chain.Add("Ds2317inDs2317_m6w4.root");
// chain.Add("Ds2317inDs2317_m6w4s1.root");
chain.Add("../Ds2317inDs2317_X4274-like-identical.root");

 Float_t
  m_dsp, hel_phi,
  m_ds17, e_ds17, px_ds17, py_ds17, pz_ds17, pst_ds17,
  e_dsi, px_dsi, py_dsi, pz_dsi,
  m_dsii, dgf_dsii, c2_dsii, cn_dsii,
  dgf_pi0, c2_pi0;

chain.SetBranchAddress("m_dsp", &m_dsp);
chain.SetBranchAddress("m_dsii", &m_dsii);
chain.SetBranchAddress("e_dsi", &e_dsi);
chain.SetBranchAddress("px_dsi", &px_dsi);
chain.SetBranchAddress("py_dsi", &py_dsi);
chain.SetBranchAddress("pz_dsi", &pz_dsi);
chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
chain.SetBranchAddress("c2_dsii", &c2_dsii);
chain.SetBranchAddress("cn_dsii", &cn_dsii);
chain.SetBranchAddress("m_ds17", &m_ds17);
chain.SetBranchAddress("pst_ds17", &pst_ds17);
chain.SetBranchAddress("e_ds17", &e_ds17);
chain.SetBranchAddress("px_ds17", &px_ds17);
chain.SetBranchAddress("py_ds17", &py_ds17);
chain.SetBranchAddress("pz_ds17", &pz_ds17);
chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
chain.SetBranchAddress("c2_pi0", &c2_pi0);
chain.SetBranchAddress("hel_phi", &hel_phi);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  float m_var = sqrt(pow(e_dsi+e_ds17,2)-pow(px_dsi+px_ds17,2)-pow(py_dsi+py_ds17,2)-pow(pz_dsi+pz_ds17,2));
  chain.GetEntry(i);
  if(m_dsp>range0 && m_dsp<range1 && m_ds17>x0_ds17 && m_ds17<x1_ds17 && TMath::Abs(hel_phi)>0.35 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii!=2){
    mass->setVal(m_dsp);
    data->add(RooArgSet(*mass));
  }}

 TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);

 RooPlot *frame = mass->frame();
 data->plotOn(frame);

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
 frame->SetMinimum(350);
 //frame->GetXaxis()->SetTitleOffset(0.8);
 // frame->SetFillColor(kYellow-8);
 /*
 TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
 txt->SetBorderSize(0);
 txt->SetFillColor(0);
 txt->SetTextSize(0.05);
 txt->SetTextFont(80);
 txt->AddText("(d)") ;
 frame->addObject(txt);
 */
 //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
 //pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));
 //pdf->paramOn(frame, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
 //frame->getAttText()->SetTextSize(0.025);
 //frame->getAttLine()->SetLineWidth(0);

 //pdf->plotOn(frame);
 //pdf->plotOn(frame,Components(RooArgSet(*gaus1)),LineColor(kYellow-2),LineStyle(kDashed));
 //pdf->plotOn(frame,Components(RooArgSet(*Pol)),LineColor(kRed),LineStyle(kDashed));
 //pdf->plotOn(frame,Components(RooArgSet(*gaus2)),LineColor(kGreen+3),LineStyle(kDashed));
 frame->Draw();

 cv->SaveAs("SMC_DspinDs2317_X-4274-like-identical_dsp.pdf","pdf");
}
