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

void Fit_Ds2460_mass(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  //Float_t range0=1.91, range1=2.05;
  Float_t range0=2.3, range1=2.6;
  //Float_t x0=2.32, x1=2.57;
  Float_t x0=2.44, x1=2.48;
  //Float_t x0=2.4, x1=2.5;
  RooRealVar *mass = new RooRealVar("mass","",range0,range1);
  TChain chain("h1");
  chain.Add("../Ds2460.root");

  Float_t m_dsw, pst_dsw, dgf_dsw, c2_dsw, cn_dsw, m_dsn, pst_dsn, dgf_dsn, c2_dsn, cn_dsn, cn_dsw, m_ds60, mc_ds60, pst_ds60, dgf_ds60, c2_ds60, cn_ds60, p_pi0, dgf_pi0, c2_pi0, m_phin, c2_phin, dgf_phin, c2_phiw, dgf_phiw, c2_pi0n, dgf_pi0n, p_pi0n, c2_pi0w, dgf_pi0w, p_pi0w, m_kst0, c2_kst0, dgf_kst0, p_kst0, m_dsst, pst_dsst, dgf_dsst, c2_dsst, mctr, pstg_dss, hel_phi, hel_ds60, cn_dsn, cos60ds; 

  chain.SetBranchAddress("m_dsw", &m_dsw);
  chain.SetBranchAddress("pst_dsw", &pst_dsw);
  chain.SetBranchAddress("dgf_dsw", &dgf_dsw);
  chain.SetBranchAddress("c2_dsw", &c2_dsw);
  chain.SetBranchAddress("cn_dsw", &cn_dsw);
  chain.SetBranchAddress("m_dsn", &m_dsn);
  chain.SetBranchAddress("pst_dsn", &pst_dsn);
  chain.SetBranchAddress("dgf_dsn", &dgf_dsn);
  chain.SetBranchAddress("c2_dsn", &c2_dsn);
  chain.SetBranchAddress("cn_dsn", &cn_dsn);
  chain.SetBranchAddress("m_ds60", &m_ds60);
  chain.SetBranchAddress("mc_ds60", &mc_ds60);
  chain.SetBranchAddress("pst_ds60", &pst_ds60);
  chain.SetBranchAddress("dgf_ds60", &dgf_ds60);
  chain.SetBranchAddress("c2_ds60", &c2_ds60);
  chain.SetBranchAddress("cn_ds60", &cn_ds60);
  chain.SetBranchAddress("p_pi0", &p_pi0);
  chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
  chain.SetBranchAddress("c2_pi0", &c2_pi0);
  chain.SetBranchAddress("m_phin", &m_phin);
  chain.SetBranchAddress("dgf_phin", &dgf_phin);
  chain.SetBranchAddress("c2_phin", &c2_phin);
  chain.SetBranchAddress("dgf_phiw", &dgf_phiw);
  chain.SetBranchAddress("c2_phiw", &c2_phiw);
  chain.SetBranchAddress("c2_pi0w", &c2_pi0w);
  chain.SetBranchAddress("dgf_pi0w", &dgf_pi0w);
  chain.SetBranchAddress("p_pi0w", &p_pi0w);
  chain.SetBranchAddress("c2_pi0n", &c2_pi0n);
  chain.SetBranchAddress("dgf_pi0n", &dgf_pi0n);
  chain.SetBranchAddress("p_pi0n", &p_pi0n);
  chain.SetBranchAddress("m_kst0", &m_kst0);
  chain.SetBranchAddress("c2_kst0", &c2_kst0);
  chain.SetBranchAddress("dgf_kst0", &dgf_kst0);
  chain.SetBranchAddress("p_kst0", &p_kst0);
  chain.SetBranchAddress("m_dsst", &m_dsst);
  chain.SetBranchAddress("pst_dsst", &pst_dsst);
  chain.SetBranchAddress("dgf_dsst", &dgf_dsst);
  chain.SetBranchAddress("c2_dsst", &c2_dsst);
  chain.SetBranchAddress("mctr", &mctr);
  chain.SetBranchAddress("pstg_dss", &pstg_dss);
  chain.SetBranchAddress("hel_phi", &hel_phi);
  chain.SetBranchAddress("hel_ds60", &hel_ds60);
  chain.SetBranchAddress("cn_dsn", &cn_dsn);
  chain.SetBranchAddress("cn_dsw", &cn_dsw);
  chain.SetBranchAddress("cos60ds", &cos60ds);

  RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if(TMath::Abs(mc_ds60)==20433 && m_ds60>range0 && m_ds60<range1 && pst_ds60>3.5 && pst_dsn>2.5 && TMath::Abs(hel_phi)>0.35 && cn_dsn==1 && cn_dsw==1 && TMath::Prob(c2_ds60, dgf_ds60)>0.001 && TMath::Prob(c2_dsn, dgf_dsn)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && TMath::Prob(c2_dsst, dgf_dsst)>0.001){
	mass->setVal(m_ds60);
	data->add(RooArgSet(*mass));
    }}
  //TMath::Abs(mc_ds60)==20433 
  
  RooRealVar *mean_gS  = new RooRealVar("sgn mean",  "sgn mean", 2.4595, 2.45, 2.47);
  RooRealVar *mean_gBS = new RooRealVar("BS mean",   "Mean of gaussian", 2.4595, 2.45, 2.47);
  RooRealVar *mean_gR  = new RooRealVar("Refl mean", "Mean of gaussian", 2.468,  2.45, 2.47);
  RooRealVar *sigma_gS  = new RooRealVar("sgn sigma",  "Sigma of gaussian", 0.0045, 0.003, 0.008);
  RooRealVar *sigma_gBS = new RooRealVar("BS sigma",   "Sigma of gaussian", 0.0195, 0.015, 0.025);
  RooRealVar *sigma_gR  = new RooRealVar("Refl sigma", "Sigma of gaussian", 0.0123, 0.01, 0.015);

  RooGaussian *gausS    = new RooGaussian("gausS",    "Gaussian PDF", *mass, *mean_gS, *sigma_gS);
  RooGaussian *gausBS   = new RooGaussian("gausBS",   "Gaussian PDF", *mass, *mean_gBS, *sigma_gBS);
  RooGaussian *gausRefl = new RooGaussian("gausRefl", "Gaussian PDF", *mass, *mean_gR, *sigma_gR);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0));//,*a1,*a2, *a3));

  RooRealVar *sig  = new RooRealVar("N_sig",  "N true signal", 0, 5E3);
  RooRealVar *refl = new RooRealVar("N_refl", "N refl signal", 0, 1E4);
  RooRealVar *BS   = new RooRealVar("N_BS",   "N broken signal", 0, 1E4);
  RooRealVar *bkg  = new RooRealVar("N_bkg",  "N ackground", 1E3, 1E5);

  //RooAddPdf *g2pol = new RooAddPdf ("g2pol", "Gaussian + Pol",RooArgList(*gaus2, *Pol), RooArgList(*bkg_g2, *bkg_pol));
  //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *g2pol), RooArgList(*sig, *bkg));

  //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gausS, *Pol), RooArgList(*sig, *bkg));
  RooAddPdf *pdf = new RooAddPdf ("pdf", "2*Gaussian + Pol", RooArgList(*gausS), RooArgList(*sig));     

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
  frame->GetXaxis()->SetTitle("M(D_{s}#pi^{0}), GeV ");
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
  /*
    TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
    txt->SetBorderSize(0);
    txt->SetFillColor(0);
    txt->SetTextSize(0.05);
    txt->SetTextFont(80);
    txt->AddText("(a)") ;
    frame->addObject(txt);
  */
  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  //pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));
  
  pdf->paramOn(frame, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
  frame->getAttText()->SetTextSize(0.025);
  frame->getAttLine()->SetLineWidth(0);
  
  pdf->plotOn(frame);
  pdf->plotOn(frame,Components(RooArgSet(*gausS)),LineColor(kYellow-2),LineStyle(kDashed));
  //pdf->plotOn(frame,Components(RooArgSet(*gausBS)),LineColor(kGreen-3),LineStyle(kDashed));
  //pdf->plotOn(frame,Components(RooArgSet(*gausRefl)),LineColor(kMagenta+3),LineStyle(kDashed));
  pdf->plotOn(frame,Components(RooArgSet(*Pol)),LineColor(kRed),LineStyle(kDashed));
  
  //pdf->plotOn(frame,Components(RooArgSet(*gaus2)),LineColor(kGreen+3),LineStyle(kDashed));
  frame->Draw();
  
  cv->SaveAs("Ds2460inDs2460_BCS_mass_True.pdf","pdf");
  fitresult->Print("v");
}
