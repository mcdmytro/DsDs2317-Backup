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

using namespace RooFit;

void FitData_PID(){
gROOT->Reset();
gSystem->Load("libRooFit");
//gStyle->SetCanvasPreferGL(1);
Float_t range0=1.9, range1=2.05;
//Float_t range0=2.12, range1=2.5;
Float_t x0=1.9, x1=2.05;
//Float_t x0=2.2, x1=2.5;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
 chain.Add("../Ds2317-continuum-PID.root");

 Float_t r2, m_dsp, pst_dsp, mctruth, pst_dsi, m_dsi, mc_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,pst_dsii, m_dsii, dgf_dsii, c2_dsii, cn_dsii, ppi, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, pipi, c2_iipi, dgf_iipi, piipi, mc_ds17, hel_phi, k1_id, k2_id, k3_id, k4_id, pi1_id, pi2_id;

chain.SetBranchAddress("m_dsi", &m_dsi);
chain.SetBranchAddress("mc_dsi", &mc_dsi);
chain.SetBranchAddress("pst_dsi", &pst_dsi);
chain.SetBranchAddress("dgf_dsi", &dgf_dsi);
chain.SetBranchAddress("c2_dsi", &c2_dsi);
chain.SetBranchAddress("cn_dsi", &cn_dsi);
chain.SetBranchAddress("m_dsii", &m_dsii);
chain.SetBranchAddress("pst_dsii", &pst_dsii);
chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
chain.SetBranchAddress("c2_dsii", &c2_dsii);
chain.SetBranchAddress("cn_dsii", &cn_dsii);
chain.SetBranchAddress("m_ds17", &m_ds17);
chain.SetBranchAddress("pst_ds17", &pst_ds17);
chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
chain.SetBranchAddress("c2_ds17", &c2_ds17);
chain.SetBranchAddress("p_pi0", &ppi);
chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
chain.SetBranchAddress("c2_pi0", &c2_pi0);
chain.SetBranchAddress("dgf_phi1", &dgf_phi1);
chain.SetBranchAddress("c2_phi1", &c2_phi1);
chain.SetBranchAddress("dgf_phi2", &dgf_phi2);
chain.SetBranchAddress("c2_phi2", &c2_phi2);
chain.SetBranchAddress("c2_ipi", &c2_ipi);
chain.SetBranchAddress("dgf_ipi", &dgf_ipi);
chain.SetBranchAddress("p_ipi", &pipi);
chain.SetBranchAddress("c2_iipi", &c2_iipi);
chain.SetBranchAddress("dgf_iipi", &dgf_iipi);
chain.SetBranchAddress("p_iipi", &piipi);
chain.SetBranchAddress("mctr", &mctruth);
chain.SetBranchAddress("mc_ds17", &mc_ds17);
chain.SetBranchAddress("hel_phi", &hel_phi);
chain.SetBranchAddress("r2", &r2);
chain.SetBranchAddress("m_dsp", &m_dsp);
chain.SetBranchAddress("pst_dsp", &pst_dsp);
chain.SetBranchAddress("k1_id", &k1_id);
chain.SetBranchAddress("k2_id", &k2_id);
chain.SetBranchAddress("pi1_id", &pi1_id);
chain.SetBranchAddress("k3_id", &k3_id);
chain.SetBranchAddress("k4_id", &k4_id);
chain.SetBranchAddress("pi2_id", &pi2_id);

RooDataSet *data1 = new RooDataSet("data","", RooArgSet(*mass),"GeV");
RooDataSet *data2 = new RooDataSet("data","", RooArgSet(*mass),"GeV");
RooDataSet *data3 = new RooDataSet("data","", RooArgSet(*mass),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(m_dsi>range0 && m_dsi<range1 && pst_ds17>2.79 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi==1 && cn_dsii==1){
    mass->setVal(m_dsi);
    data1->add(RooArgSet(*mass));
  }}

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(k3_id>0.6 && k4_id>0.6 && pi2_id<0.6 && m_dsi>range0 && m_dsi<range1 && pst_ds17>2.79 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi==1 && cn_dsii==1){
    mass->setVal(m_dsi);
    data2->add(RooArgSet(*mass));
  }}

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(k1_id>0.6 && k2_id>0.6 && k3_id>0.6 && k4_id>0.6 && pi1_id<0.6 && pi2_id<0.6 && m_dsi>range0 && m_dsi<range1 && pst_ds17>2.79 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi==1 && cn_dsii==1){
    mass->setVal(m_dsi);
    data3->add(RooArgSet(*mass));
  }}
//k1_id>0.5 && k2_id>0.2 && k3_id>0.5 && k4_id>0.2 && pi1_id<0.9 && pi2_id<0.9 &&

  RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",1.9685,1.9,2.0,"GeV");
  RooRealVar *sigma=new RooRealVar("sigma", "Sigma of gaussian",0.006,0.003,0.015,"GeV");
  RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *mass, *mean, *sigma);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1)) ;

  RooRealVar *sig = new RooRealVar("Nsig", "",0.,1E5);
  RooRealVar *bkg = new RooRealVar("Nbkg", "",0.,1E6);

  RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus, *Pol), RooArgList(*sig, *bkg));

  TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);
  // Define "signal" range in x as [x0,x1]
  mass->setRange("signal",x0,x1) ;
  // Fit pdf only to data in "signal" range
  RooFitResult *fitresult = pdf->fitTo(*data1,Save(kTRUE),Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = mass->frame(Bins(100));
  data1->plotOn(frame);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.12) ;
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("Mass, GeV ");
  //frame->GetYaxis()->SetTitle("Entries");
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetTitleOffset(1.8);
  //frame->GetXaxis()->SetTitleOffset(0.8);
  // frame->SetFillColor(kYellow-8);

  // pdf->paramOn(frame, Format("NELU", 1) ,Layout(0.65,0.99,0.99));
  // frame->getAttText()->SetTextSize(0.025);
  // frame->getAttLine()->SetLineWidth(0);
  //
  // pdf->plotOn(frame);
  // pdf->plotOn(frame,Components(RooArgSet(*gaus)),LineColor(kYellow-2),LineStyle(kDashed));
  // pdf->plotOn(frame,Components(RooArgSet(*Pol)),LineColor(kRed),LineStyle(kDashed));
  data2->plotOn(frame, Name("data2"), MarkerColor(kGreen+3));
  data3->plotOn(frame, Name("data3"), MarkerColor(kRed));

  frame->Draw();

  cv->SaveAs("CMC_Ds_continuum_PIDstudy.pdf","pdf");
}
