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
Float_t range0=2.2, range1=2.4;
Float_t x0=2.2, x1=2.4;
//Float_t x0=2.3178-5*0.005, x1=2.3178+5*0.005;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
chain.Add("../Ds2317inDs2317_BCS.root");

Float_t hel_ds17, mctr, cos0, r2, hel, pst_dsi, m_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,  pst_dsii, m_dsii, mc_dsii, dgf_dsii, c2_dsii, cn_dsii, p_pi0, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, p_ipi, c2_iipi, dgf_iipi, p_iipi, c2_kst1, dgf_kst1, p_kst1, k1_id, k2_id, k3_id, k4_id, pi1_id, pi2_id, mc_ds17, mc_dsi, cosgg, hel_phi;

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
chain.SetBranchAddress("c2_kst1", &c2_kst1);
chain.SetBranchAddress("dgf_kst1", &dgf_kst1);
chain.SetBranchAddress("p_kst1", &p_kst1);
chain.SetBranchAddress("mctr", &mctr);
chain.SetBranchAddress("cos0", &cos0);
chain.SetBranchAddress("r2", &r2);
chain.SetBranchAddress("mc_ds17", &mc_ds17);
chain.SetBranchAddress("mc_dsi", &mc_dsi);
chain.SetBranchAddress("cosgg", &cosgg);
chain.SetBranchAddress("hel_phi", &hel_phi);
chain.SetBranchAddress("hel_ds17", &hel_ds17);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(m_ds17>range0 && m_ds17<range1  && hel_ds17<0.7 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi!=2 && cn_dsii!=2){
    mass->setVal(m_ds17);
    data->add(RooArgSet(*mass));
  }}

 RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",2.312, 2.31, 2.32);
 RooRealVar *mean_g2 = new RooRealVar("Peaking bkg mean", "Mean of gaussian",2.276, 2.2, 2.3);
 //RooRealVar *sigma_gaus1=new RooRealVar("sigma", "Sigma of gaussian",0.007,0.003,0.015,"GeV");
 RooRealVar *sigma_gaus1=new RooRealVar("sigma_gaus1", "Sigma of gaussian",0.00492, 0.003, 0.015);
 RooRealVar *sigma_gaus2=new RooRealVar("Peaking bkg sigma", "Sigma of gaussian",0.007,0.003,0.015,"GeV");

 RooGaussian *gaus1= new RooGaussian("gaus1", "Gaussian PDF", *mass, *mean, *sigma_gaus1);
 RooGaussian *gaus2= new RooGaussian("gaus2", "Gaussian PDF", *mass, *mean_g2, *sigma_gaus2);

 RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
 RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
 RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
 RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
 RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1,*a2));
 RooChebychev *cpol = new RooChebychev("cpol","Chebyshev polynomial for background",*mass, RooArgList(*a0,*a1,*a2)) ;

 //RooRealVar *bkg = new RooRealVar("Nbkg", "",9672);
 //RooRealVar *bkg = new RooRealVar("N bkg", "",0.,1E4);
 RooRealVar *bkg_pol = new RooRealVar("N pol bkg", "",0.,1E5);
 RooRealVar *bkg_g2 = new RooRealVar("N peaking bkg", "",0.,1E4);
 //RooRealVar *sig = new RooRealVar("Nsig", "",5964);
 RooRealVar *sig = new RooRealVar("N sig", "",0.,1E4);

 //RooAddPdf *g2pol = new RooAddPdf ("g2pol", "Gaussian + Pol",RooArgList(*gaus2, *Pol), RooArgList(*bkg_g2, *bkg_pol));
 //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *g2pol), RooArgList(*sig, *bkg));

 //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *gaus2, *Pol), RooArgList(*sig, *bkg_g2, *bkg_pol));
 RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus1, *cpol), RooArgList(*sig, *bkg_pol));

 // Define "signal" range in x as [x0,x1]
 mass->setRange("signal",x0,x1) ;
 // Fit pdf only to data in "signal" range
 RooFitResult *fitresult = pdf->fitTo(*data, Extended(kTRUE), Minos(kFALSE), Save(kTRUE), Range("signal")) ;
 //RooFitResult *fitresult = pdf->fitTo(*data);
 RooPlot *frame = mass->frame();
 data->plotOn(frame);
 frame->SetTitle("");
 frame->GetXaxis()->SetTitle("M(D_{s}#pi^{0}) [GeV/c^{2}] ");
 frame->GetXaxis()->SetNdivisions(507);
 frame->GetXaxis()->SetLabelFont(132);
 frame->GetXaxis()->SetLabelSize(0.05);
 frame->GetXaxis()->SetTitleSize(0.07);
 frame->GetXaxis()->SetTitleFont(132);
 frame->GetXaxis()->SetTitleOffset(0.9);
 frame->GetYaxis()->SetTitle("Events / 2 MeV/c^{2}");
 frame->GetYaxis()->SetNdivisions(507);
 frame->GetYaxis()->SetLabelFont(132);
 frame->GetYaxis()->SetLabelSize(0.05);
 frame->GetYaxis()->SetTitleSize(0.07);
 frame->GetYaxis()->SetTitleOffset(1.0);
 frame->GetYaxis()->SetTitleFont(132);
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
 // pdf->paramOn(frame, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
 // frame->getAttText()->SetTextSize(0.025);
 // frame->getAttLine()->SetLineWidth(0);

 pdf->plotOn(frame);
 pdf->plotOn(frame,Components(RooArgSet(*gaus1)),LineColor(kYellow-2),LineStyle(kDashed));
 pdf->plotOn(frame,Components(RooArgSet(*cpol)),LineColor(kRed),LineStyle(kDashed));
 //pdf->plotOn(frame,Components(RooArgSet(*gaus2)),LineColor(kGreen+3),LineStyle(kDashed));

 TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);
 cv->Range(-0.2653501,-33.51433,0.2344704,183.1811);
 cv->SetFillColor(0);
 cv->SetBorderMode(0);
 cv->SetBorderSize(2);
 cv->SetLeftMargin(0.1307471);
 cv->SetRightMargin(0.06896552);
 cv->SetTopMargin(0.04661017);
 cv->SetBottomMargin(0.154661);
 cv->SetFrameBorderMode(0);
 cv->SetFrameBorderMode(0);
 gPad->SetLeftMargin(0.15);
 gPad->SetBottomMargin(0.15);

 frame->Draw();

 cv->SaveAs("./pdf_files/SMC_BCS_Ds2317inDs2317.pdf","pdf");
  fitresult->Print("v");

  // Fitresult aсquired on 23.07.2021, added in the BNv4

  // RooFitResult: minimized FCN value: -31149.1, estimated distance to minimum: 7.19968e-05
  //               covariance matrix quality: Full, accurate covariance matrix
  //               Status : MINIMIZE=0 HESSE=0
  //
  //   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  // --------------------  ------------  --------------------------  --------
  //            N pol bkg    5.0000e+04    1.0115e+03 +/-  3.50e+01  <none>
  //                N sig    5.0000e+03    2.1673e+03 +/-  4.88e+01  <none>
  //                   a0    0.0000e+00   -6.6863e-01 +/-  5.06e-02  <none>
  //                   a1    0.0000e+00   -1.2206e-01 +/-  5.78e-02  <none>
  //                   a2    0.0000e+00    8.0820e-02 +/-  5.03e-02  <none>
  //                 mean    2.3120e+00    2.3177e+00 +/-  1.12e-04  <none>
  //          sigma_gaus1    4.9200e-03    4.7891e-03 +/-  9.28e-05  <none>


}
