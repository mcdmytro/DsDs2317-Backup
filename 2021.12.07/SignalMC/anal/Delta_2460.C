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

void Delta_2460(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  Float_t range0=0.25, range1=0.45;
  Float_t x0=0.25, x1=0.45;
  // Float_t x0=0.25, x1=0.55;
  RooRealVar *mass = new RooRealVar("mass","",range0,range1);
  TChain chain("h1");
  chain.Add("../Ds2460inDs2460_BCS.root");
  // chain.Add("../Ds2460inDs2317_BCS.root");
  Float_t m_dsw, pst_dsw, dgf_dsw, c2_dsw, cn_dsw, m_dsn, pst_dsn, dgf_dsn, c2_dsn, cn_dsn, m_ds60, mc_ds60, pst_ds60, dgf_ds60, c2_ds60, cn_ds60, p_pi0, dgf_pi0, c2_pi0, c2_phin, dgf_phin, c2_phiw, dgf_phiw, c2_pi0n, dgf_pi0n, p_pi0n, c2_pi0w, dgf_pi0w, p_pi0w, c2_kstn, dgf_kstn, p_kstn, mctr, pstg_dss, hel_phi, hel_ds60,  cos60ds, m_dsst, mc_dsst, c2_dsst, dgf_dsst, egam_dss;

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
  chain.SetBranchAddress("c2_kstn", &c2_kstn);
  chain.SetBranchAddress("dgf_kstn", &dgf_kstn);
  chain.SetBranchAddress("p_kstn", &p_kstn);
  chain.SetBranchAddress("mctr", &mctr);
  chain.SetBranchAddress("pstg_dss", &pstg_dss);
  chain.SetBranchAddress("egam_dss", &egam_dss);
  chain.SetBranchAddress("hel_phi", &hel_phi);
  chain.SetBranchAddress("hel_ds60", &hel_ds60);
  chain.SetBranchAddress("cn_dsn", &cn_dsn);
  chain.SetBranchAddress("cos60ds", &cos60ds);
  chain.SetBranchAddress("m_dsst", &m_dsst);
  chain.SetBranchAddress("dgf_dsst", &dgf_dsst);
  chain.SetBranchAddress("c2_dsst", &c2_dsst);
  chain.SetBranchAddress("mc_dsst", &mc_dsst);

  RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if(TMath::Abs(mc_ds60)==20433 && (m_ds60-m_dsst)>range0 && (m_ds60-m_dsst)<range1 && hel_ds60<0.7 && TMath::Abs(hel_phi)>0.44 && cn_dsw!=2 && cn_dsn!=2 && TMath::Prob(c2_ds60, dgf_ds60)>0.001 && TMath::Prob(c2_dsn, dgf_dsn)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && TMath::Prob(c2_dsst, dgf_dsst)>0.001){
	mass->setVal(m_ds60-m_dsst);
	data->add(RooArgSet(*mass));
    }}
  //TMath::Abs(mc_ds60)==20433

  RooRealVar *mean_gS  = new RooRealVar("sgn mean",  "sgn mean", 0.35, 0.33, 0.37);
  RooRealVar *mean_gBS = new RooRealVar("BS mean",   "Mean of gaussian", 0.35, 0.33, 0.37);
  RooRealVar *mean_gR  = new RooRealVar("Refl mean", "Mean of gaussian", 0.35, 0.33, 0.37);
  RooRealVar *sigma_gS  = new RooRealVar("sgn sigma",  "Sigma of gaussian", 0.006, 0.003, 0.01);
  RooRealVar *sigma_gBS = new RooRealVar("BS sigma",   "Sigma of gaussian", 0.0195, 0.012, 0.025);
  RooRealVar *sigma_gR  = new RooRealVar("Refl sigma", "Sigma of gaussian", 0.0123, 0.01, 0.015);

  RooGaussian *gausS  = new RooGaussian("gausS",    "Gaussian PDF", *mass, *mean_gS, *sigma_gS);
  RooGaussian *gausBS = new RooGaussian("gausBS",   "Gaussian PDF", *mass, *mean_gBS, *sigma_gBS);
  RooGaussian *gausR  = new RooGaussian("gausR", "Gaussian PDF", *mass, *mean_gR, *sigma_gR);

  // RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
  // RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
  // RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  // RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  // RooRealVar *a4 = new RooRealVar("a4","a4",0.,-100.,100.);
  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-1,1);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-1,1);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-1,1);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-1,1);
  RooRealVar *a4 = new RooRealVar("a4","a4",0.,-1,1);
  RooPolynomial *pol = new RooPolynomial("pol","Polynomial for background",*mass, RooArgList(*a0,*a1,*a2));
  RooChebychev *cpol = new RooChebychev("cpol","Chebyshev polynomial for background",*mass,RooArgList(*a0,*a1,*a2));

  RooRealVar *sig  = new RooRealVar("N_sig",  "N true signal", 0, 2E3);
  RooRealVar *refl = new RooRealVar("N_refl", "N refl signal", 1E2, 100, 5E3);
  RooRealVar *BS   = new RooRealVar("N_BS",   "N broken signal", 0, 5E3);
  RooRealVar *Npkng_bkg  = new RooRealVar("Npkng_bkg",  "N peaking background", 0, 1E5);
  RooRealVar *Npol_bkg  = new RooRealVar("Npol_bkg",  "N background", 0, 1E3);

  RooRealVar *frac1 = new RooRealVar("frac1"," ", 0.1, 0.01, 1.0);
  RooRealVar *frac2 = new RooRealVar("frac2"," ", 0.1, 0.01, 1.0);

  //RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gausS, *gausBS, *gausRefl, *cpol), RooArgList(*sig, *BS, *refl, *bkg));
  // RooAddPdf *pkng_bkg = new RooAddPdf ("pkng_bkg", "2 bkg Gaussians",RooArgList(*gausBS, *gausR), *frac1);
  // RooAddPdf *pdf = new RooAddPdf ("pdf", "2 bkg Gaussian + Pol", RooArgList(*gausBS, *gausR, *cpol), RooArgList(*frac1, *frac2));
  // RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol", RooArgList(*gausR, *cpol), RooArgList(*refl, *Npol_bkg));
  // RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol", RooArgList(*gausS, *gausBS, *cpol), RooArgList(*sig, *BS, *Npol_bkg));
  RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol", RooArgList(*gausS, *cpol), RooArgList(*sig, *Npol_bkg));
  // RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol", RooArgList(*gausBS, *cpol), RooArgList(*BS, *Npol_bkg));

  // Define "signal" range in x as [x0,x1]
  mass->setRange("signal",x0,x1);
  // Fit pdf only to data in "signal" range
  RooFitResult *fitresult = pdf->fitTo(*data, Extended(kFALSE), Minos(kFALSE), Save(kTRUE), Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = mass->frame();
  data->plotOn(frame);
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("M(D_{s}*#pi^{0})-M(D_{s}*) [GeV/c^{2}] ");
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

  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  //pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));

  // pdf->paramOn(frame, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
  // frame->getAttText()->SetTextSize(0.025);
  // frame->getAttLine()->SetLineWidth(0);

  pdf->plotOn(frame);
  pdf->plotOn(frame,Components(RooArgSet(*gausS)),LineColor(kYellow-2),LineStyle(kDashed));
  // pdf->plotOn(frame,Components(RooArgSet(*gausBS)),LineColor(kGreen-3),LineStyle(kDashed));
  // pdf->plotOn(frame,Components(RooArgSet(*gausR)),LineColor(kMagenta+3),LineStyle(kDashed));
  pdf->plotOn(frame,Components(RooArgSet(*cpol)),LineColor(kRed),LineStyle(kDashed));

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

  // cv->SaveAs("./pdf_files/SMC_BCS_Ds2460inDs2317_delta_fitted.pdf","pdf");
  cv->SaveAs("./pdf_files/SMC_BCS_Ds2460inDs2460_delta_fitted_true.pdf","pdf");

  fitresult->Print("v");

  // Fitresult aсquired on 23.07.2021, added in the BNv4

      // Ds2460 in Ds2317

          // RooFitResult: minimized FCN value: -1613.07, estimated distance to minimum: 3.94358e-06
          //               covariance matrix quality: Full, accurate covariance matrix
          //               Status : MINIMIZE=0 HESSE=0
          //
          //   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
          // --------------------  ------------  --------------------------  --------
          //               N_refl    1.0000e+02    6.2298e+02 +/-  2.70e+03  <none>
          //             Npol_bkg    5.0000e+02    2.9307e+02 +/-  5.56e+02  <none>
          //            Refl mean    3.5000e-01    3.5133e-01 +/-  7.81e-04  <none>
          //           Refl sigma    1.2300e-02    1.4605e-02 +/-  6.74e-04  <none>
          //                   a0    0.0000e+00   -8.5073e-01 +/-  1.11e-01  <none>
          //                   a1    0.0000e+00    2.6574e-01 +/-  1.57e-01  <none>
          //                   a2    0.0000e+00    3.7684e-02 +/-  1.03e-01  <none>



      // Ds2460 in Ds2460

          // RooFitResult: minimized FCN value: -5588.22, estimated distance to minimum: 1.31966e-06
          //               covariance matrix quality: Full, accurate covariance matrix
          //               Status : MINIMIZE=0 HESSE=0
          //
          //   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
          // --------------------  ------------  --------------------------  --------
          //              BS mean    3.5000e-01    3.4417e-01 +/-  9.60e-04  <none>
          //             BS sigma    1.9500e-02    1.2954e-02 +/-  1.25e-03  <none>
          //                 N_BS    2.5000e+03    4.1167e+02 +/-  2.97e+03  <none>
          //                N_sig    1.0000e+03    6.8766e+02 +/-  1.63e+03  <none>
          //             Npol_bkg    5.0000e+02    5.8543e+02 +/-  5.92e+02  <none>
          //                   a0    0.0000e+00   -8.5284e-01 +/-  5.94e-02  <none>
          //                   a1    0.0000e+00    6.5575e-02 +/-  9.14e-02  <none>
          //                   a2    0.0000e+00   -4.2759e-02 +/-  6.55e-02  <none>
          //             sgn mean    3.5000e-01    3.4764e-01 +/-  1.98e-04  <none>
          //            sgn sigma    6.0000e-03    4.0262e-03 +/-  2.13e-04  <none>

      // Ds2460 in Ds2460 True

          // RooFitResult: minimized FCN value: -4374.3, estimated distance to minimum: 0.00018956
          //               covariance matrix quality: Full, accurate covariance matrix
          //               Status : MINIMIZE=0 HESSE=0
          //
          //   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
          // --------------------  ------------  --------------------------  --------
          //                N_sig    1.0000e+03    8.5563e+02 +/-  1.57e+03  <none>
          //             Npol_bkg    5.0000e+02    2.0467e+01 +/-  1.34e+02  <none>
          //                   a0    0.0000e+00   -9.9999e-01 +/-  2.38e-01  <none>
          //                   a1    0.0000e+00   -1.0000e+00 +/-  1.82e-01  <none>
          //                   a2    0.0000e+00    9.9999e-01 +/-  1.14e-01  <none>
          //             sgn mean    3.5000e-01    3.4720e-01 +/-  1.56e-04  <none>
          //            sgn sigma    6.0000e-03    5.0648e-03 +/-  1.27e-04  <none>

      // Ds2460 in Ds2460 NOT True

          // RooFitResult: minimized FCN value: -1995.09, estimated distance to minimum: 3.67023e-06
          //               covariance matrix quality: Full, accurate covariance matrix
          //               Status : MINIMIZE=0 HESSE=0
          //
          //   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
          // --------------------  ------------  --------------------------  --------
          //              BS mean    3.5000e-01    3.4392e-01 +/-  1.56e-03  <none>
          //             BS sigma    1.9500e-02    1.6939e-02 +/-  1.83e-03  <none>
          //                 N_BS    2.5000e+03    1.2066e+02 +/-  5.67e+02  <none>
          //             Npol_bkg    5.0000e+02    2.5036e+02 +/-  6.18e+02  <none>
          //                   a0    0.0000e+00   -8.7704e-01 +/-  6.62e-02  <none>
          //                   a1    0.0000e+00    1.3054e-01 +/-  1.12e-01  <none>
          //                   a2    0.0000e+00   -7.5543e-02 +/-  7.38e-02  <none>
}
