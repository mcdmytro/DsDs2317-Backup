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

void Delta_2317(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  Float_t range0=0.3, range1=0.4;
  //Float_t x0=0.22, x1=0.5; 
  Float_t x0=0.32, x1=0.38; 
  RooRealVar *mass = new RooRealVar("mass","",range0,range1);
  TChain chain("h1");
  chain.Add("../Ds2317.root");

  Float_t mctr, cos0, r2, hel, pst_dsi, m_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,  pst_dsii, m_dsii, mc_dsii, dgf_dsii, c2_dsii, cn_dsii, p_pi0, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, p_ipi, c2_iipi, dgf_iipi, p_iipi, c2_kst0, dgf_kst0, p_kst0, k1_id, k2_id, k3_id, k4_id, pi1_id, pi2_id, mc_ds17, mc_dsi, cosgg, hel_phi;

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
  chain.SetBranchAddress("mc_dsi", &mc_dsi);
  chain.SetBranchAddress("cosgg", &cosgg);
  chain.SetBranchAddress("hel_phi", &hel_phi);

  RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if(TMath::Abs(mc_ds17)==10431 && (m_ds17-m_dsii)>range0 && (m_ds17-m_dsii)<range1 && pst_ds17>3.5 && TMath::Abs(hel_phi)>0.35 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii!=2){
      mass->setVal(m_ds17-m_dsii);
      data->add(RooArgSet(*mass));
    }}
  
  RooRealVar *mean_gS  = new RooRealVar("sgn mean",  "sgn mean", 0.35, 0.2, 0.5);
  RooRealVar *mean_gBS = new RooRealVar("BS mean",   "Mean of gaussian", 0.35, 0.2, 0.5);
  RooRealVar *mean_gR  = new RooRealVar("Refl mean", "Mean of gaussian", 0.35,  0.2, 0.5);
  RooRealVar *sigma_gS  = new RooRealVar("sgn sigma",  "Sigma of gaussian", 0.006, 0.003, 0.01);
  RooRealVar *sigma_gBS = new RooRealVar("BS sigma",   "Sigma of gaussian", 0.0195, 0.017, 0.025);
  RooRealVar *sigma_gR  = new RooRealVar("Refl sigma", "Sigma of gaussian", 0.0123, 0.01, 0.02);

  RooGaussian *gausS  = new RooGaussian("gausS",    "Gaussian PDF", *mass, *mean_gS, *sigma_gS);
  RooGaussian *gausBS = new RooGaussian("gausBS",   "Gaussian PDF", *mass, *mean_gBS, *sigma_gBS);
  RooGaussian *gausR  = new RooGaussian("gausR", "Gaussian PDF", *mass, *mean_gR, *sigma_gR);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1,*a2,*a3));

  RooRealVar *sig  = new RooRealVar("N_sig",  "N true signal", 0, 5E3);
  RooRealVar *refl = new RooRealVar("N_refl", "N refl signal", 0, 5E3);
  RooRealVar *BS   = new RooRealVar("N_BS",   "N broken signal", 0, 5E3);
  RooRealVar *bkg  = new RooRealVar("N_bkg",  "N background", 0);//, 1E5);

  //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gausS, *gausBS, *gausRefl, *Pol), RooArgList(*sig, *BS, *refl, *bkg));
  //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gausBS, *gausRefl, *Pol), RooArgList(*BS, *refl, *bkg));
  RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol", RooArgList(*gausS, *Pol), RooArgList(*sig, *bkg));     

  TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);
  // Define "signal" range in x as [x0,x1]
  mass->setRange("signal",x0,x1) ;
  // Fit pdf only to data in "signal" range
  RooFitResult *fitresult = pdf->fitTo(*data, Extended(kTRUE), Minos(kFALSE), Save(kTRUE), Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = mass->frame();
  data->plotOn(frame);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.12) ;
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("M(D_{s}#pi^{0})-M(D_{s}), GeV ");
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
  
  TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  txt->SetBorderSize(0);
  txt->SetFillColor(0);
  txt->SetTextSize(0.05);
  txt->SetTextFont(80);
  txt->AddText("(a)") ;
  frame->addObject(txt);
  
  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  //pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));
  /*
  pdf->paramOn(frame, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
  frame->getAttText()->SetTextSize(0.025);
  frame->getAttLine()->SetLineWidth(0);
  */
  pdf->plotOn(frame);
  pdf->plotOn(frame,Components(RooArgSet(*gausS)),LineColor(kYellow-2),LineStyle(kDashed));
  //pdf->plotOn(frame,Components(RooArgSet(*gausBS)),LineColor(kGreen-3),LineStyle(kDashed));
  //pdf->plotOn(frame,Components(RooArgSet(*gausR)),LineColor(kMagenta+3),LineStyle(kDashed));
  pdf->plotOn(frame,Components(RooArgSet(*Pol)),LineColor(kRed),LineStyle(kDashed));
  
  frame->Draw();
  
  cv->SaveAs("CMC_Ds2317_delta_True_referrence.pdf","pdf");
  fitresult->Print("v");
}
