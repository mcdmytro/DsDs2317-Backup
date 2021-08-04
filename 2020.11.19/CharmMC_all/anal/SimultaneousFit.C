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

void SimultaneousFit(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  gStyle->SetLineWidth(0.2);  
  TCanvas *cv = new TCanvas("cv","Just Canvas",5,5,1900,800);
  cv->Divide(2,1);
  
  Float_t
    /*
    r0_2317=2.17,
    r1_2317=2.50,
    x0_2317=2.2,
    x1_2317=2.45,
    r0_2460=2.30,
    r1_2460=2.60,
    x0_2460=2.36,
    x1_2460=2.55;
    */
    r0_2317=r0_2460=0.10,
    r1_2317=r1_2460=0.60,
    x0_2317=x0_2460=0.2,
    x1_2317=x1_2460=0.5;

  RooRealVar *mass_2317 = new RooRealVar("mass_2317", "", r0_2317, r1_2317);
  RooRealVar *mass_2460 = new RooRealVar("mass_2460", "", r0_2460, r1_2460);

  TChain chain("h1");
  chain.Add("../Ds2317.root");
  chain.Add("../Ds2460.root");

  Float_t 
    m_ds17, cn_dsii, cn_dsi, mctr, hel_phi, hel_ds17, pst_ds17, c2_ds17, dgf_ds17, c2_dsii, dgf_dsii, c2_pi0, dgf_pi0, pst_dsii, m_dsii, hel_ds17,
    cn_dsn, cn_dsw, m_ds60, mc_ds60, pstg_dss, hel_ds60, pst_ds60, pstg_dss, c2_ds60, dgf_ds60, c2_dsn, dgf_dsn, c2_dsst, dgf_dsst, pst_dsn, m_dsst, hel_ds60, pstg_dss;
  
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

  RooDataSet *data_2317 = new RooDataSet("data_2317", "", RooArgSet(*mass_2317),"GeV");
  RooDataSet *data_2460 = new RooDataSet("data_2460", "", RooArgSet(*mass_2460),"GeV");

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if((m_ds17-m_dsii)>r0_2317 && (m_ds17-m_dsii)<r1_2317 && TMath::Abs(hel_phi)>0.35 && cn_dsii!=2 && pst_ds17>3.5 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01){
      mass_2317->setVal(m_ds17-m_dsii);
      data_2317->add(RooArgSet(*mass_2317));
    }}
    
  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if((m_ds60-m_dsst)>r0_2460 && (m_ds60-m_dsst)<r1_2460 && TMath::Abs(hel_phi)>0.35 && cn_dsn!=2 && pst_ds60>3.5 && TMath::Prob(c2_ds60, dgf_ds60)>0.001 && TMath::Prob(c2_dsn, dgf_dsn)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && TMath::Prob(c2_dsst, dgf_dsst)>0.001){
      mass_2460->setVal(m_ds60-m_dsst);
	data_2460->add(RooArgSet(*mass_2460));
      }}
  //TMath::Abs(mc_ds60)==20433 
    
  //RooRealVar *mean_2317        = new RooRealVar("sgn mean", "sgn mean", 2.317, 2.31, 2.32);
  //RooRealVar *mean_g2_2317     = new RooRealVar("refl mean", "refl mean", 2.317, 2.3, 2.32);
  RooRealVar *mean_2317        = new RooRealVar("sgn mean", "sgn mean", 0.3487, 0.32, 0.36);
  RooRealVar *mean_g2_2317     = new RooRealVar("refl mean", "refl mean", 0.3444, 0.3, 0.4);
  RooRealVar *sigma_gaus1_2317 = new RooRealVar("sgn sigma", "Sigma of gaussian",0.0047. 0.003, 0.01);
  RooRealVar *sigma_gaus2_2317 = new RooRealVar("refl sigma", "Sigma of gaussian",0.01492, 0.01, 0.02);
  
  RooGaussian *gaus1_2317 = new RooGaussian("gaus1_2317", "Gaussian PDF", *mass_2317, *mean_2317, *sigma_gaus1_2317);
  RooGaussian *gaus2_2317 = new RooGaussian("gaus2_2317", "Gaussian PDF", *mass_2317, *mean_g2_2317, *sigma_gaus2_2317);

  RooRealVar *a0_2317 = new RooRealVar("a0","",0.,-100.,100.);
  RooRealVar *a1_2317 = new RooRealVar("a1","",0.,-100.,100.);
  RooRealVar *a2_2317 = new RooRealVar("a2","",0.,-100.,100.);
  RooRealVar *a3_2317 = new RooRealVar("a3","",0.,-100.,100.);

  RooPolynomial *pol_2317 = new RooPolynomial("pol_2317","Polynomial for background",*mass_2317, RooArgList(*a0_2317,*a1_2317,*a2_2317,*a3_2317));

  RooRealVar *bkg_pol_2317    = new RooRealVar("N pol bkg", "",0.,1E6);
  RooRealVar *sig_2317        = new RooRealVar("N_sig", "N true signal", 3893, 300, 5000);
  RooRealVar *sig_2460        = new RooRealVar("N_sig", "N true signal", 1300, 1000, 3700);
  //RooRealVar *refl_2317       = new RooRealVar("Nrefl", "Nrefl", 0, 5000);
  //My correlarion
  RooRealVar *refl_2317       = new RooRealVar("Nrefl", "Nrefl", sig_2460.getVal()*1.918);
  //Mikami correlation
  //RooRealVar *refl_2317       = new RooRealVar("Nrefl", "Nrefl", sig_2460.getVal()*1.32);

  //RooRealVar *mean_2460        = new RooRealVar("sgn mean", "sgn mean", 2.4595, 2.45, 2.47);
  //RooRealVar *mean_g2_2460     = new RooRealVar("BS mean", "Mean of gaussian",2.4595, 2.44, 2.47);
  //RooRealVar *mean_g3_2460     = new RooRealVar("refl mean", "Mean of gaussian",2.4595, 2.44, 2.47);
  RooRealVar *mean_2460        = new RooRealVar("sgn mean", "sgn mean", 0.3475, 0.3, 0.38);
  RooRealVar *mean_g2_2460     = new RooRealVar("refl mean", "Mean of gaussian", 0.3519, 0.3, 0.4);
  RooRealVar *mean_g3_2460     = new RooRealVar("BS mean", "Mean of gaussian", 0.3454, 0.3, 0.4);
  
  RooRealVar *sigma_gaus1_2460 = new RooRealVar("sgn sigma", "Sigma of gaussian",0.004, 0.003, 0.01);
  RooRealVar *sigma_gaus2_2460 = new RooRealVar("refl sigma", "Sigma of gaussian",0.0123, 0.01, 0.02);
  RooRealVar *sigma_gaus3_2460 = new RooRealVar("BS sigma", "Sigma of gaussian",0.01952, 0.015, 0.025);

  RooGaussian *gaus1_2460 = new RooGaussian("gaus1", "Gaussian PDF", *mass_2460, *mean_2460, *sigma_gaus1_2460);
  RooGaussian *gaus2_2460 = new RooGaussian("gaus2", "Gaussian PDF", *mass_2460, *mean_g2_2460, *sigma_gaus2_2460);
  RooGaussian *gaus3_2460 = new RooGaussian("gaus3", "Gaussian PDF", *mass_2460, *mean_g3_2460, *sigma_gaus3_2460);

  RooRealVar *a0_2460 = new RooRealVar("a0","",0.,-100.,100.);
  RooRealVar *a1_2460 = new RooRealVar("a1","",0.,-100.,100.);
  RooRealVar *a2_2460 = new RooRealVar("a2","",0.,-100.,100.);
  RooRealVar *a3_2460 = new RooRealVar("a3","",0.,-100.,100.);
  
  RooPolynomial *pol_2460 = new RooPolynomial("pol_2460","Polynomial for background",*mass_2460, RooArgList(*a0_2460,*a1_2460,*a2_2460, *a3_2460));

  RooRealVar *bkg_pol_2460    = new RooRealVar("N pol bkg", "",0.,1E6); 
  //My correlation
  RooRealVar *broken_sig_2460 = new RooRealVar("N_BS", "N BS", sig_2460.getVal()*0.450);
  RooRealVar *refl_2460       = new RooRealVar("Nrefl", "Nrefl", sig_2317.getVal()*0.259);
  //Mikami correlation
  //RooRealVar *broken_sig_2460 = new RooRealVar("N_BS", "N NS", sig_2460.getVal()*0.156);
  //RooRealVar *refl_2460       = new RooRealVar("Nrefl", "Nrefl", sig_2317.getVal()*0.092);
  //No correlation
  //RooRealVar *broken_sig_2460 = new RooRealVar("N_BS", "N NS", 200, 100, 800);
  //RooRealVar *refl_2460       = new RooRealVar("Nrefl", "Nrefl", 3000, 2000, 5000);
 
  RooAddPdf *pdf_2317 = new RooAddPdf ("pdf_2317", "Gaussian + Pol",RooArgList(*gaus1_2317, *gaus2_2317, *pol_2317), RooArgList(*sig_2317, *refl_2317, *bkg_pol_2317));
  RooAddPdf *pdf_2460 = new RooAddPdf ("pdf_2460", "Gaussian + Pol",RooArgList(*gaus1_2460, *gaus2_2460, *gaus3_2460, *pol_2460), RooArgList(*sig_2460, *refl_2460, *broken_sig_2460, *bkg_pol_2460));
 
  mass_2317->setRange("signal_2317",x0_2317,x1_2317);
  mass_2460->setRange("signal_2460",x0_2460,x1_2460);
  
  RooFitResult *fitresult_2317 = pdf_2317->fitTo(*data_2317, Extended(kTRUE), Minos(kFALSE), Save(kTRUE), Range("signal_2317")) ;
  RooFitResult *fitresult_2460 = pdf_2460->fitTo(*data_2460, Extended(kTRUE), Minos(kFALSE), Save(kTRUE), Range("signal_2460")) ;
  
  RooPlot *frame_2317 = mass_2317->frame();
  data_2317->plotOn(frame_2317);
  frame_2317->SetTitle("");
  frame_2317->GetXaxis()->SetTitle("M(D_{s}#pi^{0}), GeV ");
  frame_2317->GetYaxis()->SetTitle(" ");
  frame_2317->GetXaxis()->SetTitleSize(0.05);
  frame_2317->GetYaxis()->SetTitleSize(0.05);
  frame_2317->GetXaxis()->SetLabelSize(0.05);
  frame_2317->GetYaxis()->SetLabelSize(0.05);
  frame_2317->GetXaxis()->SetNdivisions(505);
  frame_2317->GetYaxis()->SetNdivisions(505);
  frame_2317->GetYaxis()->SetTitleOffset(1.2);
  frame_2317->GetXaxis()->SetTitleOffset(0.9);
  // frame_2317->SetFillColor(kYellow-8);
  
  RooPlot *frame_2460 = mass_2460->frame();
  data_2460->plotOn(frame_2460);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.12) ;
  frame_2460->SetTitle("");
  frame_2460->GetXaxis()->SetTitle("M(D_{s}*#pi^{0}), GeV ");
  frame_2460->GetYaxis()->SetTitle(" ");
  frame_2460->GetXaxis()->SetTitleSize(0.05);
  frame_2460->GetYaxis()->SetTitleSize(0.05);
  frame_2460->GetXaxis()->SetLabelSize(0.05);
  frame_2460->GetYaxis()->SetLabelSize(0.05);
  frame_2460->GetXaxis()->SetNdivisions(505);
  frame_2460->GetYaxis()->SetNdivisions(505);
  frame_2460->GetYaxis()->SetTitleOffset(1.2);
  frame_2460->GetXaxis()->SetTitleOffset(0.9);
  // frame_2460->SetFillColor(kYellow-8);
  /*
  pdf_2317->paramOn(frame_2317, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
  frame_2317->getAttText()->SetTextSize(0.025);
  frame_2317->getAttLine()->SetLineWidth(0);

  pdf_2460->paramOn(frame_2460, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
  frame_2460->getAttText()->SetTextSize(0.025);
  frame_2460->getAttLine()->SetLineWidth(0);
  */
  pdf_2317->plotOn(frame_2317, LineWidth(1.5));
  pdf_2317->plotOn(frame_2317,Components(RooArgSet(*gaus1_2317)),LineColor(kMagenta+3),LineStyle(kDashed), LineWidth(1.5));
  pdf_2317->plotOn(frame_2317,Components(RooArgSet(*gaus2_2317)),LineColor(kGreen-3),LineStyle(kDashed), LineWidth(1.5));
  pdf_2317->plotOn(frame_2317,Components(RooArgSet(*pol_2317)),LineColor(kRed),LineStyle(kDashed), LineWidth(1.5));  

  pdf_2460->plotOn(frame_2460, LineWidth(1.5));
  pdf_2460->plotOn(frame_2460,Components(RooArgSet(*gaus1_2460)),LineColor(kMagenta+3),LineStyle(kDashed), LineWidth(1.5));
  pdf_2460->plotOn(frame_2460,Components(RooArgSet(*gaus2_2460)),LineColor(kCyan+3),LineStyle(kDashed), LineWidth(1.5));
  pdf_2460->plotOn(frame_2460,Components(RooArgSet(*gaus3_2460)),LineColor(kGreen-3),LineStyle(kDashed), LineWidth(1.5));
  pdf_2460->plotOn(frame_2460,Components(RooArgSet(*pol_2460)),LineColor(kRed),LineStyle(kDashed), LineWidth(1.5));
  
  cv->cd(1); 
  
  TPaveText* txt1 = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");                                                                                                                                                     
  txt1->SetBorderSize(0); 
  txt1->SetFillColor(0);                                     
  txt1->SetTextSize(0.0);
  txt1->SetTextFont(80);
  txt1->AddText("(a)");                                                                                                                                                                                           
  frame_2317->addObject(txt1);
  
  frame_2317->Draw();
  
  cv->cd(2);
  
  TPaveText* txt2 = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  txt2->SetBorderSize(0);
  txt2->SetFillColor(0);
  txt2->SetTextSize(0.0);
  txt2->SetTextFont(80);
  txt2->AddText("(b)");
  frame_2460->addObject(txt2);
  
  frame_2460->Draw();

  cv->SaveAs("CMC_SimultaneousFit_ready.pdf","pdf");
  fitresult_2317->Print("v");
  fitresult_2460->Print("v");
}
