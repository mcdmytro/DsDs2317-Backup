#include "RooGlobalFunc.h"
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
#include "BWShape.h"

using namespace RooFit;

void FitDsp_subtracted(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  Float_t range0=4.286, range1=4.296;
  // Float_t range0=4.3, range1=5.1;
  Float_t x0=4.286, x1=4.296;
  Float_t x0_ds17=2.305, x1_ds17=2.329;
  RooRealVar *mass = new RooRealVar("mass","",range0,range1);
  TChain chain("h1");
  chain.Add("../Ds2317inDs2317_BCS.root");

  // TFile *f = new TFile("../Dsp.root","recreate");
  // TH1D *xmass = new TH1D("xmass","xmass",100,4.285,4.295);

  Float_t m_dsp, hel_ds17, mctr, cos0, r2, hel, pst_dsi, m_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,  pst_dsii, m_dsii, mc_dsii, dgf_dsii, c2_dsii, cn_dsii, p_pi0, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, p_ipi, c2_iipi, dgf_iipi, p_iipi, c2_kst1, dgf_kst1, p_kst1, k1_id, k2_id, k3_id, k4_id, pi1_id, pi2_id, mc_ds17, mc_dsi, cosgg, hel_phi;

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
  chain.SetBranchAddress("m_dsp", &m_dsp);

  RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

  for (int i=0; i<chain.GetEntries();i++){
  float m_var=m_dsp-m_dsi-m_ds17+1.9683+2.3178;
  chain.GetEntry(i);
  if(m_var>range0 && m_var<range1 && hel_ds17<0.7 && (m_ds17-m_dsii)<range1 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii!=2){
    mass->setVal(m_var);
    data->add(RooArgSet(*mass));
  }}

  RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",4.287, 4.285, 4.290);
  RooRealVar *width = new RooRealVar("width", "Sigma of gaussian",1E-4, 1E-6, 0.015);
  RooRealVar *gamma = new RooRealVar("gamma","width", 1E-4, 1E-6, 0.015);
  RooRealVar *mass1 = new RooRealVar("mass1", "mass1",1.96830);  //Ds mass
  RooRealVar *mass2 = new RooRealVar("mass2", "mass2",1.96830);  //Ds mass
  RooRealVar *mass3 = new RooRealVar("mass3", "mass3",0.13498);  //pi0 mass
  RooRealVar *lorb = new RooRealVar("lorb", "lorb",0);
  RooRealVar *decayType = new RooRealVar("decayType", "decayType",1);
  RooRealVar *useBPhaseSpace = new RooRealVar("useBPhaseSpace", "useBPhaseSpace",1);

  RooGaussian *gaus  = new RooGaussian("gaus","Gaussian PDF", *mass, *mean, *width);
  BWShape *brew = new BWShape("brew","Breit-Wigner",*mass,*mean,*width,*mass1,*mass2,*mass3,*lorb,*decayType,*useBPhaseSpace);
  RooVoigtian *voigt = new RooVoigtian("voigt","Reconstructed peak", *mass, *mean, *gamma, *width);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  RooRealVar *a4 = new RooRealVar("a4","a4",0.,-100.,100.);
  RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1,*a2,*a3,*a4));
  RooChebychev *cpol = new RooChebychev("cpol","Chebyshev polynomial for background",*mass,RooArgList(*a0,*a1,*a2,*a3,*a4));

  RooRealVar *sig = new RooRealVar("N sig", "",0.,1E4);
  RooRealVar *bkg = new RooRealVar("N pol bkg", "",0.,1E5);

  RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*voigt, *cpol), RooArgList(*sig, *bkg));

  // Define "signal" range in x as [x0,x1]
  mass->setRange("signal",x0,x1) ;
  // Fit pdf only to data in "signal" range
  RooFitResult *fitresult = pdf->fitTo(*data, Extended(kTRUE), Minos(kFALSE), Save(kTRUE), Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = mass->frame();
  data->plotOn(frame);

  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("M(D_{s}D_{s0}^{*}(2317)) [GeV/2^{2}]");
  frame->GetXaxis()->SetNdivisions(507);
  frame->GetXaxis()->SetLabelFont(132);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleSize(0.07);
  frame->GetXaxis()->SetTitleFont(132);
  frame->GetXaxis()->SetTitleOffset(0.9);
  frame->GetYaxis()->SetTitle("Events / 0.1 MeV/c^{2}");
  frame->GetYaxis()->SetNdivisions(507);
  frame->GetYaxis()->SetLabelFont(132);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.07);
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleFont(132);
  //frame->SetMinimum(350);
  //frame->GetXaxis()->SetTitleOffset(0.8);
  // frame->SetFillColor(kYellow-8);

  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  // pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));
  //pdf->paramOn(frame, Format("NELU", 1) ,Layout(0.6,0.99,0.99));
  // frame->getAttText()->SetTextSize(0.025);
  // frame->getAttLine()->SetLineWidth(0);

  pdf->plotOn(frame);
  pdf->plotOn(frame,Components(RooArgSet(*brew)),LineColor(kYellow-2),LineStyle(kDashed));
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

  cv->SaveAs("./pdf_files/SMC_BCS_DspinDs2317.pdf","pdf");
  fitresult->Print("v");

  // Fitresult aсquired on 23.07.2021, added in the BNv4

  // RooFitResult: minimized FCN value: -22322.7, estimated distance to minimum: 0.0490951
  //             covariance matrix quality: Full matrix, but forced positive-definite
  //             Status : MINIMIZE=-1 HESSE=4
  //
  // Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  // --------------------  ------------  --------------------------  --------
  //          N pol bkg    5.0000e+04    7.6150e+02 +/-  5.37e-01  <none>
  //              N sig    5.0000e+03    1.0938e+03 +/-  1.41e+01  <none>
  //                 a0    0.0000e+00   -9.4013e-02 +/-  1.35e-03  <none>
  //                 a1    0.0000e+00   -5.9691e-01 +/-  1.39e-03  <none>
  //                 a2    0.0000e+00    5.0097e-01 +/-  1.23e-03  <none>
  //                 a3    0.0000e+00   -4.2823e-01 +/-  1.13e-03  <none>
  //                 a4    0.0000e+00    2.0404e-01 +/-  8.99e-06  <none>
  //              gamma    1.0000e-04    9.5860e-05 +/-  1.81e-07  <none>
  //               mean    4.2870e+00    4.2868e+00 +/-  9.30e-08  <none>
  //              width    1.0000e-04    1.1015e-04 +/-  1.03e-07  <none>


}
