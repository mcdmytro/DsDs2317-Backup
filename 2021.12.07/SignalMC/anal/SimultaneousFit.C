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

  Float_t
    r0_2317=0.25,
    r0_2460=0.25,
    r1_2317=0.45,
    r1_2460=0.45,
    x0_2317=0.25,
    x0_2460=0.25,
    x1_2317=0.45,
    x1_2460=0.45;

  RooRealVar *mass_2317 = new RooRealVar("mass_2317", "", r0_2317, r1_2317);
  RooRealVar *mass_2460 = new RooRealVar("mass_2460", "", r0_2460, r1_2460);

  TChain chain1("h1");
  TChain chain2("h1");
  chain1.Add("../Ds2317inDs2317_BCS.root");
  chain2.Add("../Ds2460inDs2460_BCS.root");

  Float_t
    m_ds17, c2_ds17, dgf_ds17, cn_dsii, c2_dsii, dgf_dsii, cn_dsi, c2_dsi, dgf_dsi, c2_pi0_1, dgf_pi0_1, hel_phi_1, hel_ds17, m_dsii,
    cn_dsn, c2_dsn, dgf_dsn, cn_dsw, c2_dsw, dgf_dsw, m_ds60, mc_ds60, c2_ds60, dgf_ds60, pstg_dss, hel_ds60, m_dsst, c2_dsst, dgf_dsst, c2_pi0_2, dgf_pi0_2, hel_phi_2;

  chain1.SetBranchAddress("m_ds17", &m_ds17);
  chain1.SetBranchAddress("c2_ds17", &c2_ds17);
  chain1.SetBranchAddress("dgf_ds17", &dgf_ds17);
  chain1.SetBranchAddress("c2_pi0", &c2_pi0_1);
  chain1.SetBranchAddress("dgf_pi0", &dgf_pi0_1);
  chain1.SetBranchAddress("m_dsii", &m_dsii);
  chain1.SetBranchAddress("cn_dsii", &cn_dsii);
  chain1.SetBranchAddress("c2_dsii", &c2_dsii);
  chain1.SetBranchAddress("dgf_dsii", &dgf_dsii);
  chain1.SetBranchAddress("cn_dsi", &cn_dsi);
  chain1.SetBranchAddress("c2_dsi", &c2_dsi);
  chain1.SetBranchAddress("dgf_dsi", &dgf_dsi);
  chain1.SetBranchAddress("hel_phi", &hel_phi_1);
  chain1.SetBranchAddress("hel_ds17", &hel_ds17);

  chain2.SetBranchAddress("cn_dsn", &cn_dsn);
  chain2.SetBranchAddress("c2_dsn", &c2_dsn);
  chain2.SetBranchAddress("dgf_dsn", &dgf_dsn);
  chain2.SetBranchAddress("cn_dsw", &cn_dsw);
  chain2.SetBranchAddress("c2_dsw", &c2_dsw);
  chain2.SetBranchAddress("dgf_dsw", &dgf_dsw);
  chain2.SetBranchAddress("m_dsst", &m_dsst);
  chain2.SetBranchAddress("c2_dsst", &c2_dsst);
  chain2.SetBranchAddress("dgf_dsst", &dgf_dsst);
  chain2.SetBranchAddress("m_ds60", &m_ds60);
  chain2.SetBranchAddress("c2_ds60", &c2_ds60);
  chain2.SetBranchAddress("dgf_ds60", &dgf_ds60);
  chain2.SetBranchAddress("mc_ds60", &mc_ds60);
  chain2.SetBranchAddress("pstg_dss", &pstg_dss);
  chain2.SetBranchAddress("hel_ds60", &hel_ds60);
  chain2.SetBranchAddress("c2_pi0", &c2_pi0_2);
  chain2.SetBranchAddress("dgf_pi0", &dgf_pi0_2);
  chain2.SetBranchAddress("hel_phi", &hel_phi_2);

  RooDataSet *data_2317 = new RooDataSet("data_2317", "", RooArgSet(*mass_2317));
  RooDataSet *data_2460 = new RooDataSet("data_2460", "", RooArgSet(*mass_2460));

  for (int i=0; i<chain1.GetEntries();i++){
    chain1.GetEntry(i);
    if((m_ds17-m_dsii)>r0_2317 && (m_ds17-m_dsii)<r1_2317 && TMath::Abs(hel_phi_1)>0.44 && hel_ds17<0.7  && cn_dsi!=2 && cn_dsii!=2 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0_1, dgf_pi0_1)>0.01){
      mass_2317->setVal(m_ds17-m_dsii);
      data_2317->add(RooArgSet(*mass_2317));
    }}

  for (int i=0; i<chain2.GetEntries();i++){
    chain2.GetEntry(i);
    if((m_ds60-m_dsst)>r0_2460 && (m_ds60-m_dsst)<r1_2460 && TMath::Abs(hel_phi_2)>0.44 && hel_ds60<0.7 && cn_dsw!=2 && cn_dsn!=2 && TMath::Prob(c2_ds60, dgf_ds60)>0.001 && TMath::Prob(c2_dsn, dgf_dsn)>0.001 && TMath::Prob(c2_pi0_2, dgf_pi0_2)>0.01 && TMath::Prob(c2_dsst, dgf_dsst)>0.001){
      mass_2460->setVal(m_ds60-m_dsst);
	    data_2460->add(RooArgSet(*mass_2460));
      }}
  //TMath::Abs(mc_ds60)==20433

  //RooRealVar *mean_2317        = new RooRealVar("sgn mean", "sgn mean", 2.317, 2.31, 2.32);
  RooRealVar *mean_2317        = new RooRealVar("sgn mean", "sgn mean", 0.35, 0.2, 0.5);
  RooRealVar *sigma_gaus1_2317 = new RooRealVar("sgn sigma", "Sigma of gaussian",0.00492, 0.003,0.015);

  RooGaussian *gaus1_2317 = new RooGaussian("gaus1_2317", "Gaussian PDF", *mass_2317, *mean_2317, *sigma_gaus1_2317);

  RooRealVar *a0_2317 = new RooRealVar("a0","",0.,-100.,100.);
  RooRealVar *a1_2317 = new RooRealVar("a1","",0.,-100.,100.);
  RooRealVar *a2_2317 = new RooRealVar("a2","",0.,-100.,100.);
  RooRealVar *a3_2317 = new RooRealVar("a3","",0.,-100.,100.);

  RooPolynomial *pol_2317 = new RooPolynomial("pol_2317","Polynomial for background",*mass_2317, RooArgList(*a0_2317,*a1_2317,*a2_2317,*a3_2317));
  RooChebychev *cpol_2317 = new RooChebychev("cpol_2317","Chebyshev polynomial for background",*mass_2317,RooArgList(*a0_2317,*a1_2317,*a2_2317));

  RooRealVar *bkg_pol_2317    = new RooRealVar("N pol bkg", "",0.,1E5);
  RooRealVar *sig_2317        = new RooRealVar("N_sig", "N true signal", 0.,1E4);


  //RooRealVar *mean_2460        = new RooRealVar("sgn mean", "sgn mean", 2.4595, 2.45, 2.47);
  //RooRealVar *mean_g2_2460     = new RooRealVar("BS mean", "Mean of gaussian",2.4595, 2.44, 2.47);
  RooRealVar *mean_2460        = new RooRealVar("sgn mean", "sgn mean", 0.35, 0.2, 0.5);
  RooRealVar *mean_g2_2460     = new RooRealVar("BS mean", "Mean of gaussian",0.35, 0.2, 0.5);
  RooRealVar *sigma_gaus1_2460 = new RooRealVar("sgn sigma", "Sigma of gaussian",0.00538, 0.003, 0.01);
  RooRealVar *sigma_gaus2_2460 = new RooRealVar("BS sigma", "Sigma of gaussian",0.0195, 0.01, 0.03);

  RooGaussian *gaus1_2460 = new RooGaussian("gaus1", "Gaussian PDF", *mass_2460, *mean_2460, *sigma_gaus1_2460);
  RooGaussian *gaus2_2460 = new RooGaussian("gaus2", "Gaussian PDF", *mass_2460, *mean_g2_2460, *sigma_gaus2_2460);

  RooRealVar *a0_2460 = new RooRealVar("a0","",0.,-100.,100.);
  RooRealVar *a1_2460 = new RooRealVar("a1","",0.,-100.,100.);
  RooRealVar *a2_2460 = new RooRealVar("a2","",0.,-100.,100.);
  RooRealVar *a3_2460 = new RooRealVar("a3","",0.,-100.,100.);

  RooPolynomial *pol_2460 = new RooPolynomial("pol_2460","Polynomial for background",*mass_2460, RooArgList(*a0_2460,*a1_2460,*a2_2460, *a3_2460));
  RooChebychev *cpol_2460 = new RooChebychev("cpol_2460","Chebyshev polynomial for background",*mass_2460,RooArgList(*a0_2460,*a1_2460,*a2_2460));

  RooRealVar *bkg_pol_2460    = new RooRealVar("N pol bkg", "",0.,1E5);
  RooRealVar *sig_2460        = new RooRealVar("N_sig", "N true signal", 2800, 1E3, 5E3);
  RooRealVar *broken_sig_2460 = new RooRealVar("N_BS", "N NS", sig_2460->getVal()*0.141);

  RooAddPdf *pdf_2317 = new RooAddPdf ("pdf_2317", "Gaussian + Pol",RooArgList(*gaus1_2317, *cpol_2317), RooArgList(*sig_2317, *bkg_pol_2317));
  RooAddPdf *pdf_2460 = new RooAddPdf ("pdf_2460", "Gaussian + Pol",RooArgList(*gaus1_2460, *gaus2_2460, *cpol_2460), RooArgList(*sig_2460, *broken_sig_2460, *bkg_pol_2460));

  mass_2317->setRange("signal_2317",x0_2317,x1_2317);
  mass_2460->setRange("signal_2460",x0_2460,x1_2460);

  RooFitResult *fitresult_2317 = pdf_2317->fitTo(*data_2317, Extended(kTRUE), Minos(kFALSE), Save(kTRUE), Range("signal_2317")) ;
  RooFitResult *fitresult_2460 = pdf_2460->fitTo(*data_2460, Extended(kTRUE), Minos(kFALSE), Save(kTRUE), Range("signal_2460")) ;

  RooPlot *frame_2317 = mass_2317->frame();
  data_2317->plotOn(frame_2317);
  frame_2317->SetTitle("");
  frame_2317->GetXaxis()->SetTitle("M(D_{s}#pi^{0}) [GeV/c^{2}]");
  frame_2317->GetXaxis()->SetNdivisions(507);
  frame_2317->GetXaxis()->SetLabelFont(132);
  frame_2317->GetXaxis()->SetLabelSize(0.05);
  frame_2317->GetXaxis()->SetTitleSize(0.07);
  frame_2317->GetXaxis()->SetTitleFont(132);
  frame_2317->GetXaxis()->SetTitleOffset(0.9);
  frame_2317->GetYaxis()->SetTitle("Events / 2 MeV/c^{2}");
  frame_2317->GetYaxis()->SetNdivisions(507);
  frame_2317->GetYaxis()->SetLabelFont(132);
  frame_2317->GetYaxis()->SetLabelSize(0.05);
  frame_2317->GetYaxis()->SetTitleSize(0.07);
  frame_2317->GetYaxis()->SetTitleOffset(1.0);
  frame_2317->GetYaxis()->SetTitleFont(132);
  // frame_2317->SetFillColor(kYellow-8);
  /*
    TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
    txt->SetBorderSize(0);
    txt->SetFillColor(0);
    txt->SetTextSize(0.05);
    txt->SetTextFont(80);
    txt->AddText("(a)") ;
    frame->addObject(txt);
  */

  RooPlot *frame_2460 = mass_2460->frame();
  data_2460->plotOn(frame_2460);
  frame_2460->SetTitle("");
  frame_2460->GetXaxis()->SetTitle("M(D_{s}*#pi^{0}) [GeV/c^{2}]");
  frame_2460->GetXaxis()->SetNdivisions(507);
  frame_2460->GetXaxis()->SetLabelFont(132);
  frame_2460->GetXaxis()->SetLabelSize(0.05);
  frame_2460->GetXaxis()->SetTitleSize(0.07);
  frame_2460->GetXaxis()->SetTitleFont(132);
  frame_2460->GetXaxis()->SetTitleOffset(0.9);
  frame_2460->GetYaxis()->SetTitle("Events / 2 MeV/c^{2}");
  frame_2460->GetYaxis()->SetNdivisions(507);
  frame_2460->GetYaxis()->SetLabelFont(132);
  frame_2460->GetYaxis()->SetLabelSize(0.05);
  frame_2460->GetYaxis()->SetTitleSize(0.07);
  frame_2460->GetYaxis()->SetTitleOffset(1.0);
  frame_2460->GetYaxis()->SetTitleFont(132);
  // frame_2460->SetFillColor(kYellow-8);


  // pdf_2317->paramOn(frame_2317, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
  // frame_2317->getAttText()->SetTextSize(0.025);
  // frame_2317->getAttLine()->SetLineWidth(0);
  //
  // pdf_2460->paramOn(frame_2460, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
  // frame_2460->getAttText()->SetTextSize(0.025);
  // frame_2460->getAttLine()->SetLineWidth(0);

  pdf_2317->plotOn(frame_2317, LineWidth(1.));
  pdf_2317->plotOn(frame_2317,Components(RooArgSet(*gaus1_2317)),LineColor(kYellow-2),LineStyle(kDashed), LineWidth(1.));
  pdf_2317->plotOn(frame_2317,Components(RooArgSet(*cpol_2317)),LineColor(kRed),LineStyle(kDashed), LineWidth(1.));

  pdf_2460->plotOn(frame_2460, LineWidth(1.));
  pdf_2460->plotOn(frame_2460,Components(RooArgSet(*gaus1_2460)),LineColor(kYellow-2),LineStyle(kDashed), LineWidth(1.));
  pdf_2460->plotOn(frame_2460,Components(RooArgSet(*gaus2_2460)),LineColor(kGreen-3),LineStyle(kDashed), LineWidth(1.));
  pdf_2460->plotOn(frame_2460,Components(RooArgSet(*cpol_2460)),LineColor(kRed),LineStyle(kDashed), LineWidth(1.));

  TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,1600,800);
  cv->Divide(2,1);
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
  // gStyle->SetLineWidth(0.2);

  cv->cd(1);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  /*
  TPaveText* txt1 = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  txt1->SetBorderSize(0);
  txt1->SetFillColor(0);
  txt1->SetTextSize(0.05);
  txt1->SetTextFont(80);
  txt1->AddText("(a)") ;
  frame_2317->addObject(txt1);
  */
  frame_2317->Draw();

  cv->cd(2);
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);
  /*
  TPaveText* txt2 = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  txt2->SetBorderSize(0);
  txt2->SetFillColor(0);
  txt2->SetTextSize(0.05);
  txt2->SetTextFont(80);
  txt2->AddText("(b)") ;
  frame_2460->addObject(txt2);
  */
  frame_2460->Draw();

  cv->SaveAs("./pdf_files/SMC_BCS_SimultaneousFit.pdf","pdf");
  fitresult_2317->Print("v");
  fitresult_2460->Print("v");

  // Fitresult aсquired on 23.07.2021, added in the BNv4

  // RooFitResult: minimized FCN value: -30134, estimated distance to minimum: 0.00010072
  //             covariance matrix quality: Full, accurate covariance matrix
  //             Status : MINIMIZE=0 HESSE=0
  //
  // Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  // --------------------  ------------  --------------------------  --------
  //          N pol bkg    5.0000e+04    8.9492e+02 +/-  3.34e+01  <none>
  //              N_sig    5.0000e+03    2.1607e+03 +/-  4.88e+01  <none>
  //                 a0    0.0000e+00   -8.0068e-01 +/-  5.21e-02  <none>
  //                 a1    0.0000e+00   -1.1860e-01 +/-  6.40e-02  <none>
  //                 a2    0.0000e+00    8.3329e-02 +/-  5.44e-02  <none>
  //           sgn mean    3.5000e-01    3.4913e-01 +/-  1.12e-04  <none>
  //          sgn sigma    4.9200e-03    4.7728e-03 +/-  9.27e-05  <none>
  //
  //
  // RooFitResult: minimized FCN value: -20894.8, estimated distance to minimum: 8.99259e-06
  //             covariance matrix quality: Full, accurate covariance matrix
  //             Status : MINIMIZE=0 HESSE=0
  //
  // Constant Parameter    Value
  // --------------------  ------------
  //               N_BS    3.9480e+02
  //
  // Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  // --------------------  ------------  --------------------------  --------
  //            BS mean    3.5000e-01    3.4338e-01 +/-  1.25e-03  <none>
  //           BS sigma    1.9500e-02    1.3693e-02 +/-  1.31e-03  <none>
  //          N pol bkg    5.0000e+04    8.2811e+02 +/-  4.14e+01  <none>
  //              N_sig    2.8000e+03    1.0346e+03 +/-  4.50e+01  <none>
  //                 a0    0.0000e+00   -8.2083e-01 +/-  5.49e-02  <none>
  //                 a1    0.0000e+00   -1.2016e-02 +/-  7.90e-02  <none>
  //                 a2    0.0000e+00   -5.5219e-03 +/-  6.12e-02  <none>
  //           sgn mean    3.5000e-01    3.4754e-01 +/-  1.89e-04  <none>
  //          sgn sigma    5.3800e-03    4.3520e-03 +/-  1.70e-04  <none>

}
