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

void FitData_17_PID(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  Float_t range0=2.1, range1=2.6;
  Float_t x0=2.2, x1=2.55;
  RooRealVar *mass = new RooRealVar("mass","",range0,range1);
  TChain chain("h1");
  chain.Add("../Ds2317-continuum-PID.root");

  Float_t m_dsp, pst_dsp, mctr, cos0, r2, hel, pst_dsi, m_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,  pst_dsii, m_dsii, mc_dsii, dgf_dsii, c2_dsii, cn_dsii, p_pi0, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, p_ipi, c2_iipi, dgf_iipi, p_iipi, c2_kst0, dgf_kst0, p_kst0, k1_id, k2_id, k3_id, k4_id, pi1_id, pi2_id, mc_ds17, cosgg, hel_phi;

  chain.SetBranchAddress("m_dsi", &m_dsi);
  chain.SetBranchAddress("pst_dsi", &pst_dsi);
  chain.SetBranchAddress("dgf_dsi", &dgf_dsi);
  chain.SetBranchAddress("c2_dsi", &c2_dsi);
  chain.SetBranchAddress("cn_dsi", &cn_dsi);
  chain.SetBranchAddress("m_dsii", &m_dsii);
  chain.SetBranchAddress("mc_dsii", &mc_dsii);
  chain.SetBranchAddress("pst_dsii", &pst_dsii);
  chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
  chain.SetBranchAddress("c2_dsii", &c2_dsii);
  chain.SetBranchAddress("cn_dsii", &cn_dsii);
  chain.SetBranchAddress("m_ds17", &m_ds17);
  chain.SetBranchAddress("pst_ds17", &pst_ds17);
  chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
  chain.SetBranchAddress("c2_ds17", &c2_ds17);
  chain.SetBranchAddress("p_pi0", &p_pi0);
  chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
  chain.SetBranchAddress("c2_pi0", &c2_pi0);
  chain.SetBranchAddress("dgf_phi1", &dgf_phi1);
  chain.SetBranchAddress("c2_phi1", &c2_phi1);
  chain.SetBranchAddress("dgf_phi2", &dgf_phi2);
  chain.SetBranchAddress("c2_phi2", &c2_phi2);
  chain.SetBranchAddress("c2_ipi", &c2_ipi);
  chain.SetBranchAddress("dgf_ipi", &dgf_ipi);
  chain.SetBranchAddress("p_ipi", &p_ipi);
  chain.SetBranchAddress("c2_iipi", &c2_iipi);
  chain.SetBranchAddress("dgf_iipi", &dgf_iipi);
  chain.SetBranchAddress("p_iipi", &p_iipi);
  chain.SetBranchAddress("c2_kst0", &c2_kst0);
  chain.SetBranchAddress("dgf_kst0", &dgf_kst0);
  chain.SetBranchAddress("p_kst0", &p_kst0);
  chain.SetBranchAddress("mctr", &mctr);
  chain.SetBranchAddress("cos0", &cos0);
  chain.SetBranchAddress("r2", &r2);
  chain.SetBranchAddress("mc_ds17", &mc_ds17);
  chain.SetBranchAddress("cosgg", &cosgg);
  chain.SetBranchAddress("hel_phi", &hel_phi);
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
    if(m_ds17>range0 && m_ds17<range1 && pst_ds17>2.79 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi==1 && cn_dsii==1){
      mass->setVal(m_ds17);
      data1->add(RooArgSet(*mass));
  }}

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if(k1_id>0.6 && k2_id>0.6 && pi1_id<0.6 && m_ds17>range0 && m_ds17<range1 && pst_ds17>2.79 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi==1 && cn_dsii==1){
      mass->setVal(m_ds17);
      data2->add(RooArgSet(*mass));
  }}

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if(k1_id>0.6 && k2_id>0.6 && k3_id>0.6 && k4_id>0.6 && pi1_id<0.6 && pi2_id<0.6 && m_ds17>range0 && m_ds17<range1 && pst_ds17>2.79 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi==1 && cn_dsii==1){
      mass->setVal(m_ds17);
      data3->add(RooArgSet(*mass));
  }}
  //k1_id>0.6 && k2_id>0.6 && k3_id>0.6 && k4_id>0.6 && pi1_id<0.6 && pi2_id<0.6 &&

  RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",2.3178, 2.31,2.32);
  RooRealVar *mean_g2 = new RooRealVar("bkg mean", "Mean of gaussian",2.3149, 2.31, 2.32);
  RooRealVar *sigma_gaus1=new RooRealVar("#sigma", "Sigma of gaussian",0.0068, 0.003, 0.02);
  RooRealVar *sigma_gaus2=new RooRealVar("bkg #sigma", "Sigma of gaussian",0.0152, 0.01, 0.02);

  RooGaussian *gaus1= new RooGaussian("gaus1", "Gaussian PDF", *mass, *mean, *sigma_gaus1);
  RooGaussian *gaus2= new RooGaussian("gaus2", "Gaussian PDF", *mass, *mean_g2, *sigma_gaus2);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1,*a2,*a3)) ;

  RooRealVar *bkg_pol = new RooRealVar("N pol bkg", "", 0.,1E7);
  RooRealVar *bkg_g2 = new RooRealVar("N peaking bkg", "", 1637);
  RooRealVar *sig = new RooRealVar("N sig", "",2194, 1, 1E5);
  RooRealVar *refl = new RooRealVar("N refl", "",0., 1E4);

  RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *Pol), RooArgList(*sig, *bkg_pol));

  TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);

  // mass->setRange("signal",x0,x1) ;

  // RooFitResult *fitresult = pdf->fitTo(*data, Extended(kTRUE), Minos(kFALSE), Save(kTRUE), Range("signal")) ;

  RooPlot *frame = mass->frame();
  data1->plotOn(frame, Name("data1"));
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

  // pdf->plotOn(frame);
  // pdf->plotOn(frame,Components(RooArgSet(*gaus1)),LineColor(kYellow-2),LineStyle(kDashed));
  // pdf->plotOn(frame,Components(RooArgSet(*Pol)),LineColor(kRed),LineStyle(kDashed));

  data2->plotOn(frame, Name("data2"), MarkerColor(kGreen+3));
  data3->plotOn(frame, Name("data3"), MarkerColor(kRed));

  frame->Draw();

  cv->SaveAs("CMC_Ds2317_continuum_PIDstudy.pdf","pdf");
  // cv->SaveAs("CMC_Ds2317_test.pdf","pdf");
}
