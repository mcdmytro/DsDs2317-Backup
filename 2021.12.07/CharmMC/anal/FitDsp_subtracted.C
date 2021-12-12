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

void FitDsp_subtracted(){
gROOT->Reset();
gSystem->Load("libRooFit");
Float_t range0=3.8, range1=6.;
Float_t x0=2.215, x1=2.35;
Float_t x0_ds17=2.305, x1_ds17=2.329;
Float_t width=0.008;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
chain.Add("../Ds2317.root");

 Float_t m_dsp, hel_phi, pst_ds17, m_ds17, m_dsii, dgf_dsii, c2_dsii, cn_dsii, dgf_pi0, c2_pi0;

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

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  float m_var=m_dsp-m_dsii-m_ds17+1.9683+2.3178;
  chain.GetEntry(i);
  if(m_var>range0 && m_var<range1 && m_ds17>x0_ds17 && m_ds17<x1_ds17 && pst_ds17>3.5 && TMath::Abs(hel_phi)>0.35 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii!=2){
    mass->setVal(m_var);
    data->add(RooArgSet(*mass));
  }}

 RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",2.312, 2.31, 2.32);
 RooRealVar *mean_g2 = new RooRealVar("Peaking bkg mean", "Mean of gaussian",2.276, 2.2, 2.3);
 //RooRealVar *sigma_gaus1=new RooRealVar("sigma", "Sigma of gaussian",0.007,0.003,0.015,"GeV");
 RooRealVar *sigma_gaus1=new RooRealVar("sigma_gaus1", "Sigma of gaussian",0.00492, 0.003, 0.015);
 RooRealVar *sigma_gaus2=new RooRealVar("Peaking bkg sigma", "Sigma of gaussian",0.007,0.003,0.015,"GeV");
 RooRealVar *sigma_bw = new RooRealVar("sigma_bw","Sigma of Breit Wigner",width);
 //RooRealVar *sigma_bw = new RooRealVar("sigma_bw","Sigma of Breit Wigner",0..0214,0.05,"GeV");

 RooGaussian *gaus1= new RooGaussian("gaus1", "Gaussian PDF", *mass, *mean, *sigma_gaus1);
 RooGaussian *gaus2= new RooGaussian("gaus2", "Gaussian PDF", *mass, *mean_g2, *sigma_gaus2);
 RooBreitWigner *bw = new RooBreitWigner("bw","Breit-Wigner PDF",*mass,*mean, *sigma_bw);
 RooVoigtian *voigt = new RooVoigtian("voigt","Voigt PDF", *mass, *mean, *sigma_bw, *sigma_gaus1);

 //RooRealVar *a0 = new RooRealVar("a0","a0",83.034);
 //RooRealVar *a1 = new RooRealVar("a1","a1",-33.88);
 RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
 RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
 RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
 RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
 RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1,*a2)) ;

 //RooRealVar *bkg = new RooRealVar("Nbkg", "",9672);
 //RooRealVar *bkg = new RooRealVar("N bkg", "",0.,1E4);
 RooRealVar *bkg_pol = new RooRealVar("N pol bkg", "",0.,1E5);
 RooRealVar *bkg_g2 = new RooRealVar("N peaking bkg", "",0.,1E4);
 //RooRealVar *sig = new RooRealVar("Nsig", "",5964);
 RooRealVar *sig = new RooRealVar("N sig", "",0.,1E4);

 //RooAddPdf *g2pol = new RooAddPdf ("g2pol", "Gaussian + Pol",RooArgList(*gaus2, *Pol), RooArgList(*bkg_g2, *bkg_pol));
 //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *g2pol), RooArgList(*sig, *bkg));

 //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *gaus2, *Pol), RooArgList(*sig, *bkg_g2, *bkg_pol));
 RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *Pol), RooArgList(*sig, *bkg_pol));

 TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);
 // Define "signal" range in x as [x0,x1]
 mass->setRange("signal",x0,x1) ;
 // Fit pdf only to data in "signal" range
 //RooFitResult *fitresult = pdf->fitTo(*data, Extended(kTRUE), Minos(kFALSE), Save(kTRUE), Range("signal")) ;
 //RooFitResult *fitresult = pdf->fitTo(*data);
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
 //frame->SetMinimum(350);
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

 cv->SaveAs("CMC_DspinDs2317_subtracted_4-6range.pdf","pdf");
 //fitresult->Print("v");
}
