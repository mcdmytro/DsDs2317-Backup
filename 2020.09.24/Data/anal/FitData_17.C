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

void FitData_17(){
gROOT->Reset();
gSystem->Load("libRooFit");
//Float_t range0=1.91, range1=2.05;
Float_t range0=2.12, range1=2.5;
//Float_t x0=1.91, x1=2.05;
Float_t x0=2.2, x1=2.5;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
chain.Add("../3modesPID_all_new.root");

 Float_t mctruth, cos, R2, helicity,pst_dsi, m_dsi, dgf_dsi, chi2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, chi2_ds17,pst_dsii, m_dsii, dgf_dsii, chi2_dsii, cn_dsii, ppi, dgf_pi, chi2_pi, dgf_phi1, chi2_phi1, dgf_phi2, chi2_phi2, chi2_ipi, dgf_ipi, pipi, chi2_iipi, dgf_iipi, piipi, chi2_kst, dgf_kst, pkst;

chain.SetBranchAddress("m_dsi", &m_dsi);
chain.SetBranchAddress("pst_dsi", &pst_dsi);
chain.SetBranchAddress("dgf_dsi", &dgf_dsi);
chain.SetBranchAddress("c2_dsi", &chi2_dsi);
chain.SetBranchAddress("cn_dsi", &cn_dsi);
chain.SetBranchAddress("m_dsii", &m_dsii);
chain.SetBranchAddress("pst_dsii", &pst_dsii);
chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
chain.SetBranchAddress("c2_dsii", &chi2_dsii);
chain.SetBranchAddress("cn_dsii", &cn_dsii);
chain.SetBranchAddress("m_ds17", &m_ds17);
chain.SetBranchAddress("pst_ds17", &pst_ds17);
chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
chain.SetBranchAddress("c2_ds17", &chi2_ds17);
chain.SetBranchAddress("p_pi0", &ppi);
chain.SetBranchAddress("dgf_pi0", &dgf_pi);
chain.SetBranchAddress("c2_pi0", &chi2_pi);
chain.SetBranchAddress("dgf_phi1", &dgf_phi1);
chain.SetBranchAddress("c2_phi1", &chi2_phi1);
chain.SetBranchAddress("dgf_phi2", &dgf_phi2);
chain.SetBranchAddress("c2_phi2", &chi2_phi2);
chain.SetBranchAddress("c2_ipi", &chi2_ipi);
chain.SetBranchAddress("dgf_ipi", &dgf_ipi);
chain.SetBranchAddress("p_ipi", &pipi);
chain.SetBranchAddress("c2_iipi", &chi2_iipi);
chain.SetBranchAddress("dgf_iipi", &dgf_iipi);
chain.SetBranchAddress("p_iipi", &piipi);
 chain.SetBranchAddress("c2_kst0", &chi2_kst);
 chain.SetBranchAddress("dgf_kst0", &dgf_kst);
 chain.SetBranchAddress("p_kst0", &pkst);
chain.SetBranchAddress("mctr", &mctruth);
chain.SetBranchAddress("cos0", &cos);
chain.SetBranchAddress("r2", &R2);
chain.SetBranchAddress("hel", &helicity);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(m_ds17>range0 && m_ds17<range1 && TMath::Abs(helicity)>0. && pst_ds17>3.2 && pst_dsii>2.2 && TMath::Prob(chi2_ds17, dgf_ds17)>0.001 && TMath::Prob(chi2_dsii, dgf_dsii)>0.001 && TMath::Prob(chi2_pi, dgf_pi)>0.01 && ((cn_dsii!=2) || (cn_dsii==2 && TMath::Prob(chi2_iipi, dgf_iipi)>0.1 && TMath::Prob(chi2_kst,dgf_kst)>0.01 && piipi>0.3))){
  mass->setVal(m_ds17);
  data->add(RooArgSet(*mass));
  }}

 RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",2.3178);
 RooRealVar *mean_g2 = new RooRealVar("Peaking bkg mean", "Mean of gaussian",2.302, 2.29, 2.305);
 RooRealVar *sigma_gaus1=new RooRealVar("sigma", "Sigma of gaussian",0.007,0.003,0.010,"GeV");
 //RooRealVar *sigma_gaus1=new RooRealVar("sigma_gaus1", "Sigma of gaussian",0.00492);
 RooRealVar *sigma_gaus2=new RooRealVar("Peaking bkg sigma", "Sigma of gaussian",0.0056,0.005,0.009,"GeV");
 RooRealVar *sigma_bw = new RooRealVar("sigma_bw","Sigma of Breit Wigner",0.001, 0.,0.010,"GeV");
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
 RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1)) ;

 //RooRealVar *bkg = new RooRealVar("Nbkg", "",9672);
 //RooRealVar *bkg = new RooRealVar("N bkg", "",0.,1E4);
 RooRealVar *bkg_pol = new RooRealVar("N pol bkg", "",0.,1E5); 
 RooRealVar *bkg_g2 = new RooRealVar("N peaking bkg", "",0.,1E3);
 //RooRealVar *sig = new RooRealVar("Nsig", "",5964);
 RooRealVar *sig = new RooRealVar("N sig", "",0.,1E4);
 
 RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *gaus2, *Pol), RooArgList(*sig, *bkg_g2, *bkg_pol));
  
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

 //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
 //pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));
 pdf->paramOn(frame, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
 frame->getAttText()->SetTextSize(0.025);
 frame->getAttLine()->SetLineWidth(0);

 pdf->plotOn(frame);
 pdf->plotOn(frame,Components(RooArgSet(*gaus1)),LineColor(kYellow-2),LineStyle(kDashed));
 pdf->plotOn(frame,Components(RooArgSet(*Pol)),LineColor(kRed),LineStyle(kDashed));
 pdf->plotOn(frame,Components(RooArgSet(*gaus2)),LineColor(kGreen+3),LineStyle(kDashed));
 frame->Draw();

  std::cout<<"chi^2 = "<< frame->chiSquare() <<std::endl;
  cv->SaveAs("all_PID_17.pdf","pdf");
}
