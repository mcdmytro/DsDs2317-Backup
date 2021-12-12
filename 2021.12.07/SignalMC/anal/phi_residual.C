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

void phi_residual(){
gROOT->Reset();
gSystem->Load("libRooFit");

Float_t range0=-0.01, range1=0.01;
//Float_t range0=2.12, range1=2.5;
Float_t x0=-0.01, x1=0.01;
//Float_t x0=2.2, x1=2.5;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
chain.Add("../Ds2317inDs2317_BCS_noDsKmvFit.root");

Float_t m_phi2, phi_gm;

chain.SetBranchAddress("m_phi2", &m_phi2);
chain.SetBranchAddress("phi_gm", &phi_gm);

Float_t mctr, cos0, r2, hel, pst_dsi, m_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,  pst_dsii, m_dsii, mc_dsii, dgf_dsii, c2_dsii, cn_dsii, p_pi0, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, p_ipi, c2_iipi, dgf_iipi, p_iipi, k1_id, k2_id, k3_id, k4_id, pi1_id, pi2_id, mc_ds17, mc_dsi, cosgg, hel_phi;

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
  if((m_phi2-phi_gm)>range0 && (m_phi2-phi_gm)<range1 && m_phi2>0 && phi_gm>0 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii!=2){
    mass->setVal(m_phi2-phi_gm);
    data->add(RooArgSet(*mass));
}}


  RooRealVar *mean = new RooRealVar("mean", "Mean value",0,-0.01,0.01,"GeV");
  RooRealVar *sigma_gaus = new RooRealVar("sigma_gaus", "Sigma of gaussian",0.00052,0.0003,0.015,"GeV");
  RooRealVar *sigma_bw = new RooRealVar("sigma_bw","Sigma of Breit Wigner",0.0006,0.0003,0.002,"GeV");

  RooRealVar *mean1 = new RooRealVar("mean1", "Mean value",0,-0.01,0.01,"GeV");
  RooRealVar *sigma_gaus1 = new RooRealVar("sigma_gaus1", "Sigma of gaussian",0.001,0.0001,0.015,"GeV");
  RooGaussian *gaus1= new RooGaussian("gaus1", "Gaussian for signal", *mass, *mean1, *sigma_gaus1);
  RooRealVar *frac = new RooRealVar("frac", "Fraction",0.25,0,1);

  RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *mass, *mean, *sigma_gaus);
  RooBreitWigner *bw = new RooBreitWigner("bw","Breit-Wigner function fit", *mass, *mean, *sigma_bw);
  RooVoigtian *voigt = new RooVoigtian("voigt","Reconstructed peak", *mass, *mean, *sigma_bw, *sigma_gaus);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-1.,1.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-1.,1.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-1.,1.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-1.,1.);

  RooPolynomial *pol = new RooPolynomial("pol","Polynomial for background",*mass, RooArgList(*a0,*a1));
  RooChebychev *cpol = new RooChebychev("cpol","Chebyshev polynomial for background",*mass,RooArgList(*a0,*a1,*a2));

  RooRealVar *sig = new RooRealVar("Nsig", "",0.,1E5);
  RooRealVar *bkg = new RooRealVar("Nbkg", "",0.,1E5);

  RooAddPdf *sig_pdf = new RooAddPdf ("sig_pdf", "2xGauss", RooArgList(*gaus, *gaus1), RooArgList(*frac));
  RooAddPdf *pdf = new RooAddPdf ("pdf", "Signal + background function",RooArgList(*sig_pdf, *cpol), RooArgList(*sig, *bkg));

  // RooAddPdf *pdf = new RooAddPdf ("pdf", "Signal + background function",RooArgList(*bw, *cpol), RooArgList(*sig, *bkg));
  // RooAddPdf *pdf = new RooAddPdf ("pdf", "Signal + background function",RooArgList(*gaus, *cpol), RooArgList(*sig, *bkg));

  // Define "signal" range in x as [x0,x1]
  mass->setRange("signal",x0,x1) ;
  // Fit pdf only to data in "signal" range
  RooFitResult *fitresult = pdf->fitTo(*data,Save(kTRUE),Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = mass->frame(Title(" "));
  data->plotOn(frame);
  frame->GetYaxis()->SetTitle("Entries / 0.2 MeV/c^{2} ");
  frame->GetXaxis()->SetTitle("M(#phi^{rec})-M(#phi^{gen}) [GeV/c^{2}] ");
  frame->GetXaxis()->SetNdivisions(507);
  frame->GetXaxis()->SetLabelFont(132);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleSize(0.07);
  frame->GetXaxis()->SetTitleFont(132);
  frame->GetXaxis()->SetTitleOffset(0.9);
  frame->GetYaxis()->SetNdivisions(507);
  frame->GetYaxis()->SetLabelFont(132);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.07);
  frame->GetYaxis()->SetTitleFont(132);
  frame->GetYaxis()->SetTitleOffset(1.0);
  //frame->GetXaxis()->SetTitleOffset(0.8);
  // frame->SetFillColor(kYellow-8);

  // TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  // txt->SetBorderSize(0);
  // txt->SetFillColor(0);
  // txt->SetTextSize(0.05);
  // txt->SetTextFont(80);
  // txt->AddText("(b)") ;
  // frame->addObject(txt);

  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  //pdf->paramOn(frame, Format("NELU", 1), Layout(0.6,0.9,0.9));
  //frame->getAttText()->SetTextSize(0.025);
  //frame->getAttLine()->SetLineWidth(0);

  pdf->plotOn(frame);
  /*
  pdf->plotOn(frame,
	      Components(RooArgSet(*bw)),
	      LineColor(kYellow-2),
	      LineStyle(kDashed));
  */
  pdf->plotOn(frame,
	      Components(RooArgSet(*cpol)),
	      LineColor(kRed),
	      LineStyle(kDashed));

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

  // //Find the Chi2 value
  // Double_t chi2 = frame->chiSquare();
  //
  // //Find the number of bins of the fit range
  // Int_t lower = frame->GetXaxis()->FindFixBin(x0);
  // Int_t upper = frame->GetXaxis()->FindFixBin(x1);
  // Int_t range = upper-lower;
  //
  // //Find the number of free parameters
  // Int_t dgf=(fitresult->floatParsFinal()).getSize();
  //
  // //Find the Chi2/nDOF value
  // Double_t new_chi2 = chi2/(range - dgf);
  //
  // cout << "\nChi2:" << chi2;
  // cout << "\nNumber of bins:" << range;
  // cout << "\nNumber of free parameters:" << dgf;
  // cout << "\nChi2/nDOF for (Breit Wigner + Chebyshev Polynom) fit:" << new_chi2 << endl;
  //
  // RooDataHist* datah = data->binnedClone();
  // pdf->chi2FitTo(*datah) ;
  // RooChi2Var chi2_lowstat("chi2_lowstat","chi2",*pdf,*datah) ;
  // cout << chi2_lowstat.getVal() << endl ;

  cv->SaveAs("./pdf_files/SMC_BCS_phi_residual_2xGauss_ready.pdf","pdf");
  fitresult->Print("v");

  // Fitresult aсquired on 23.07.2021, added in the BNv4

  // RooFitResult: minimized FCN value: -57513.2, estimated distance to minimum: 0.000100646
  //               covariance matrix quality: Full, accurate covariance matrix
  //               Status : MINIMIZE=0 HESSE=0
  //
  //   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  // --------------------  ------------  --------------------------  --------
  //                 Nbkg    5.0000e+04    1.0747e+03 +/-  7.34e+01  <none>
  //                 Nsig    5.0000e+04    3.5592e+03 +/-  8.87e+01  <none>
  //                   a0    0.0000e+00    4.9193e-02 +/-  4.76e-02  <none>
  //                   a1    0.0000e+00   -7.7397e-01 +/-  7.26e-02  <none>
  //                   a3    0.0000e+00    3.3374e-03 +/-  5.66e-02  <none>
  //                 frac    2.5000e-01    6.7265e-01 +/-  3.45e-02  <none>
  //                 mean    0.0000e+00   -2.3150e-05 +/-  1.48e-05  <none>
  //                mean1    0.0000e+00    1.2510e-04 +/-  8.03e-05  <none>
  //           sigma_gaus    5.2000e-04    5.0437e-04 +/-  2.37e-05  <none>
  //          sigma_gaus1    1.0000e-03    1.5918e-03 +/-  1.74e-04  <none>

}
