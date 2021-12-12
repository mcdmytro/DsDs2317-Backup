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

void Fit_Ds2460_massDelta(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  //Float_t range0=1.91, range1=2.05;
  Float_t range0=0.2, range1=0.6;
  Float_t x0=0.24, x1=0.5;
  //Float_t x0=2.37, x1=2.52;
  //Float_t x0=2.4, x1=2.5;
  RooRealVar *mass = new RooRealVar("mass","",range0,range1);
  TChain chain("h1");
  chain.Add("../Ds2317inDs2460.root");

Float_t 
  m_ds17, cn_dsii, mctr, hel_phi, hel_ds17, pst_ds17, c2_ds17, dgf_ds17, c2_dsii, dgf_dsii, c2_pi0, dgf_pi0, pst_dsii, m_dsii, hel_ds17,
  cn_dsn, m_ds60, mc_ds60, pstg_dss, hel_ds60, pst_ds60, pstg_dss, c2_ds60, dgf_ds60, c2_dsn, dgf_dsn, c2_dsst, dgf_dsst, pst_dsn, m_dsst, hel_ds60, pstg_dss;
  
 chain.SetBranchAddress("m_ds17", &m_ds17);
 chain.SetBranchAddress("m_dsii", &m_dsii);
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


  RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if((m_ds17-m_dsii)>range0 && (m_ds17-m_dsii)<range1 && TMath::Abs(hel_phi)>0.35 && hel_ds17<0.7 && cn_dsn!=2 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01){
    mass->setVal(m_ds17-m_dsii);  
    //if(TMath::Abs(mc_ds60)!=20433 && (m_ds60-m_dsst)>range0 && (m_ds60-m_dsst)<range1 && TMath::Abs(hel_phi)>0.35 && hel_ds60<0.7 && cn_dsn!=2 && TMath::Prob(c2_ds60, dgf_ds60)>0.001 && TMath::Prob(c2_dsn, dgf_dsn)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && TMath::Prob(c2_dsst, dgf_dsst)>0.001){
      //mass->setVal(m_ds60-m_dsst);
      data->add(RooArgSet(*mass));
      }}

  RooRealVar *mean = new RooRealVar("sgn mean", "sgn mean", 0.35, 0.2, 0.6);
  RooRealVar *mean_g2 = new RooRealVar("BS mean", "Mean of gaussian",0.35, 0.2, 0.6);
  RooRealVar *sigma_gaus1=new RooRealVar("sgn sigma", "Sigma of gaussian",0.005, 0.003, 0.03);
  RooRealVar *sigma_gaus2=new RooRealVar("BS sigma", "Sigma of gaussian",0.0195, 0.003, 0.03);

  RooGaussian *gaus1= new RooGaussian("gaus1", "Gaussian PDF", *mass, *mean, *sigma_gaus1);
  //RooGaussian *gaus2= new RooGaussian("gaus2", "Gaussian PDF", *mass, *mean_g2, *sigma_gaus2);
  RooGaussian *gaus2= new RooGaussian("gaus2", "Gaussian PDF", *mass, *mean, *sigma_gaus2);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1,*a2,*a3));

  RooRealVar *bkg_pol = new RooRealVar("N pol bkg", "",5E3, 1E2, 5E4); 
  RooRealVar *bkg_g2 = new RooRealVar("N peaking bkg", "",0.,1E4);
  RooRealVar *sig = new RooRealVar("N_sig", "N true signal", 5E3, 1E2, 1E4);
  //RooRealVar *sig = new RooRealVar("N_sig", "N true signal", 1E2, 1E4);
  RooRealVar *broken_sig = new RooRealVar("N_BS", "N broken signal", 1E2, 1E4);
 
  //RooAddPdf *g2pol = new RooAddPdf ("g2pol", "Gaussian + Pol",RooArgList(*gaus2, *Pol), RooArgList(*bkg_g2, *bkg_pol));
  //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *g2pol), RooArgList(*sig, *bkg));

  //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1), RooArgList(*sig));
  RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *gaus2, *Pol), RooArgList(*sig, *broken_sig, *bkg_pol));     

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
  txt->AddText("(b)") ;
  frame->addObject(txt);
  */
  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  //pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));
  
  pdf->paramOn(frame, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
  frame->getAttText()->SetTextSize(0.025);
  frame->getAttLine()->SetLineWidth(0);
  
  pdf->plotOn(frame);
  pdf->plotOn(frame,Components(RooArgSet(*gaus1)),LineColor(kYellow-2),LineStyle(kDashed));
  pdf->plotOn(frame,Components(RooArgSet(*gaus2)),LineColor(kGreen-3),LineStyle(kDashed));
  pdf->plotOn(frame,Components(RooArgSet(*Pol)),LineColor(kRed),LineStyle(kDashed));
  
  //pdf->plotOn(frame,Components(RooArgSet(*gaus2)),LineColor(kGreen+3),LineStyle(kDashed));
  frame->Draw();
  
  cv->SaveAs("Ds2317inDs2460_test2.pdf","pdf");
  fitresult->Print("v");
}
