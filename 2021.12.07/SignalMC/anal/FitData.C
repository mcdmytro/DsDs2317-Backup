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

void FitData(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  gStyle->SetCanvasPreferGL(1);
  Float_t range0=1.94, range1=2.0;
  //Float_t range0=2.12, range1=2.5;
  Float_t x0=1.94, x1=2.0;
  //Float_t x0=2.2, x1=2.5;
  RooRealVar *mass = new RooRealVar("mass","",range0,range1);
  TChain chain("h1");
  chain.Add("../Ds2317inDs2317_BCS_noDsKmvFit.root");

  Float_t r2, m_dsp, pst_dsp, mctruth, pst_dsi, m_dsi, mc_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,pst_dsii, m_dsii, dgf_dsii, c2_dsii, cn_dsii, ppi, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, pipi, c2_iipi, dgf_iipi, piipi, mc_ds17, hel_phi;

  chain.SetBranchAddress("m_dsi", &m_dsi);
  chain.SetBranchAddress("mc_dsi", &mc_dsi);
  chain.SetBranchAddress("pst_dsi", &pst_dsi);
  chain.SetBranchAddress("dgf_dsi", &dgf_dsi);
  chain.SetBranchAddress("c2_dsi", &c2_dsi);
  chain.SetBranchAddress("cn_dsi", &cn_dsi);
  chain.SetBranchAddress("m_dsii", &m_dsii);
  chain.SetBranchAddress("pst_dsii", &pst_dsii);
  chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
  chain.SetBranchAddress("c2_dsii", &c2_dsii);
  chain.SetBranchAddress("cn_dsii", &cn_dsii);
  chain.SetBranchAddress("m_ds17", &m_ds17);
  chain.SetBranchAddress("pst_ds17", &pst_ds17);
  chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
  chain.SetBranchAddress("c2_ds17", &c2_ds17);
  chain.SetBranchAddress("p_pi0", &ppi);
  chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
  chain.SetBranchAddress("c2_pi0", &c2_pi0);
  chain.SetBranchAddress("dgf_phi1", &dgf_phi1);
  chain.SetBranchAddress("c2_phi1", &c2_phi1);
  chain.SetBranchAddress("dgf_phi2", &dgf_phi2);
  chain.SetBranchAddress("c2_phi2", &c2_phi2);
  chain.SetBranchAddress("c2_ipi", &c2_ipi);
  chain.SetBranchAddress("dgf_ipi", &dgf_ipi);
  chain.SetBranchAddress("p_ipi", &pipi);
  chain.SetBranchAddress("c2_iipi", &c2_iipi);
  chain.SetBranchAddress("dgf_iipi", &dgf_iipi);
  chain.SetBranchAddress("p_iipi", &piipi);
  chain.SetBranchAddress("mctr", &mctruth);
  chain.SetBranchAddress("mc_ds17", &mc_ds17);
  chain.SetBranchAddress("hel_phi", &hel_phi);
  chain.SetBranchAddress("r2", &r2);
  chain.SetBranchAddress("m_dsp", &m_dsp);
  chain.SetBranchAddress("pst_dsp", &pst_dsp);
  RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if(cn_dsii==3 && m_dsii>range0 && m_dsii<range1 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01){// && TMath::Abs(mc_dsi)==431){ -------> && cn_dsii!=2
      mass->setVal(m_dsii);
      data->add(RooArgSet(*mass));
    }}

  RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",1.9685,1.9,2.0,"GeV");
  RooRealVar *sigma=new RooRealVar("sigma", "Sigma of gaussian",0.006,0.001,0.015,"GeV");
  RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *mass, *mean, *sigma);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1)) ;
  RooChebychev *cpol = new RooChebychev("cpol","Chebyshev polynomial for background",*mass,RooArgList(*a0,*a1,*a2));

  RooRealVar *sig = new RooRealVar("Nsig", "",100,0.,5E4);
  RooRealVar *bkg = new RooRealVar("Nbkg", "",10,0.,1E3);

  RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus, *cpol), RooArgList(*sig, *bkg));

  // Define "signal" range in x as [x0,x1]
  mass->setRange("signal",x0,x1) ;
  // Fit pdf only to data in "signal" range
  RooFitResult *fitresult = pdf->fitTo(*data,Save(kTRUE),Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = mass->frame(Bins(100));
  data->plotOn(frame);
  frame->SetTitle("");
  // frame->GetXaxis()->SetTitle("M(#phi#pi, #phi#pi#pi^{0}, K*K) [GeV/c^{2}]");
  // frame->GetXaxis()->SetTitle("M(#phi#pi) [GeV/c^{2}]");
  // frame->GetXaxis()->SetTitle("M(#phi#pi#pi^{0}) [GeV/c^{2}]");
  frame->GetXaxis()->SetTitle("M(K*K) [GeV/c^{2}]");
  frame->GetXaxis()->SetNdivisions(507);
  frame->GetXaxis()->SetLabelFont(132);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleSize(0.07);
  frame->GetXaxis()->SetTitleFont(132);
  frame->GetXaxis()->SetTitleOffset(0.9);
  frame->GetYaxis()->SetTitle("Events / 0.6 MeV/c^{2}");
  frame->GetYaxis()->SetNdivisions(507);
  frame->GetYaxis()->SetLabelFont(132);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.07);
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleFont(132);
  //frame->GetXaxis()->SetTitleOffset(0.8);
  // frame->SetFillColor(kYellow-8);

  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  //pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));


  //pdf->paramOn(frame, Format("NEU", 1) ,Layout(0.64,0.99,0.99));
  //frame->getAttText()->SetTextSize(0.025);
  //frame->getAttLine()->SetLineWidth(0);



  Double_t
    sgnl = 1.948,
    sgnr = 1.982,
    lsbl = 1.92,
    lsbr = 1.937,
    rsbl = 1.993,
    rsbr = 2.01,
    ampl = 3500;
  /*
  TLine *vl_sgn_left = new TLine(sgnl, 0, sgnl, ampl);
  vl_sgn_left->SetLineColor(kRed);
  frame->addObject(vl_sgn_left);
  TLine *vl_sgn_right = new TLine(sgnr, 0, sgnr, ampl);
  vl_sgn_right->SetLineColor(kRed);
  frame->addObject(vl_sgn_right);
  TBox *box_sgn = new TBox(sgnl, 0, sgnr, ampl);
  box_sgn->SetFillColorAlpha(kRed-4, 0.15);
  frame->addObject(box_sgn);

  TLine *vl_lsb_left = new TLine(lsbl, 0, lsbl, ampl);
  vl_lsb_left->SetLineColor(kBlue);
  frame->addObject(vl_lsb_left);
  TLine *vl_lsb_right = new TLine(lsbr,0, lsbr, ampl);
  vl_lsb_right->SetLineColor(kBlue);
  frame->addObject(vl_lsb_right);
  TBox *box_lsb = new TBox(lsbl, 0, lsbr, ampl);
  box_lsb->SetFillColorAlpha(kBlue-4, 0.15);
  frame->addObject(box_lsb);

  TLine *vl_rsb_left = new TLine(rsbl, 0, rsbl, ampl);
  vl_rsb_left->SetLineColor(kBlue);
  frame->addObject(vl_rsb_left);
  TLine *vl_rsb_right = new TLine(rsbr, 0, rsbr, ampl);
  vl_rsb_right->SetLineColor(kBlue);
  frame->addObject(vl_rsb_right);
  TBox *box_rsb = new TBox(rsbl, 0, rsbr, ampl);
  box_rsb->SetFillColorAlpha(kBlue-4, 0.15);
  frame->addObject(box_rsb);
  */

  // TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  // txt->SetBorderSize(0);
  // txt->SetFillColor(0);
  // txt->SetTextSize(0.05);
  // txt->SetTextFont(80);
  // txt->AddText("(e)") ;
  // frame->addObject(txt);

  pdf->plotOn(frame);
  pdf->plotOn(frame, Components(RooArgSet(*gaus)), LineColor(kYellow-2), LineStyle(kDashed));
  pdf->plotOn(frame, Components(RooArgSet(*cpol)), LineColor(kRed), LineStyle(kDashed));

  //Create and customize canvas
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

  frame->Draw();
  /*
  mass->setRange("pureSGN", sgnl, sgnr);
  mass->setRange("lsb", lsbl, lsbr);
  mass->setRange("rsb", rsbl, rsbr);

  RooAbsReal* int_sig = Pol->createIntegral(*mass,Range("pureSGN"));
  RooAbsReal* int_lsb = Pol->createIntegral(*mass,Range("lsb"));
  RooAbsReal* int_rsb = Pol->createIntegral(*mass,Range("rsb"));

  cout << "Signal region integral: " << '\t' << int_sig->getVal() << '\n';
  cout <<"Sideband region integral: " << '\t' << int_lsb->getVal()+int_rsb->getVal() << '\n';
  */
  //std::cout<<"chi^2 = "<< frame->chiSquare() <<std::endl;
  cv->SaveAs("./pdf_files/SMC_BCS_ds_3.pdf","pdf");
  fitresult->Print("v");

  // Fitresult aсquired on 23.07.2021, added in the BNv4

  // 0. All modes included
  // RooFitResult: minimized FCN value: -97651.5, estimated distance to minimum: 0.00111453
  //               covariance matrix quality: Full, accurate covariance matrix
  //               Status : MINIMIZE=0 HESSE=0
  //
  //   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  // --------------------  ------------  --------------------------  --------
  //                 Nbkg    5.0000e+02    9.9999e+02 +/-  1.27e+01  <none>
  //                 Nsig    2.5000e+04    7.1168e+03 +/-  8.78e+01  <none>
  //                   a0    0.0000e+00   -4.1487e-01 +/-  3.82e-02  <none>
  //                   a1    0.0000e+00   -8.4138e-01 +/-  4.39e-02  <none>
  //                   a2    0.0000e+00    4.2170e-01 +/-  4.41e-02  <none>
  //                 mean    1.9685e+00    1.9687e+00 +/-  4.25e-05  <none>
  //                sigma    6.0000e-03    3.2180e-03 +/-  3.56e-05  <none>


  // 1. cn_dsii==1

  // RooFitResult: minimized FCN value: -49215.2, estimated distance to minimum: 5.63554e-05
  //               covariance matrix quality: Full, accurate covariance matrix
  //               Status : MINIMIZE=0 HESSE=0
  //
  //   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  // --------------------  ------------  --------------------------  --------
  //                 Nbkg    5.0000e+05    4.8999e+02 +/-  3.31e+01  <none>
  //                 Nsig    5.0000e+04    3.8210e+03 +/-  6.65e+01  <none>
  //                   a0    0.0000e+00   -3.4733e-01 +/-  6.77e-02  <none>
  //                   a1    0.0000e+00   -9.2269e-01 +/-  5.08e-02  <none>
  //                   a2    0.0000e+00    4.2269e-01 +/-  8.18e-02  <none>
  //                 mean    1.9685e+00    1.9687e+00 +/-  5.55e-05  <none>
  //                sigma    6.0000e-03    3.0590e-03 +/-  5.00e-05  <none>



  // 2. cn_dsii==2

  // RooFitResult: minimized FCN value: -1419.39, estimated distance to minimum: 0.000141376
  //             covariance matrix quality: Full, accurate covariance matrix
  //             Status : MINIMIZE=0 HESSE=0
  //
  // Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  // --------------------  ------------  --------------------------  --------
  //               Nbkg    5.0000e+05    5.4579e+01 +/-  5.04e+01  <none>
  //               Nsig    5.0000e+04    1.3741e+02 +/-  5.11e+01  <none>
  //                 a0    0.0000e+00   -1.0496e+00 +/-  6.08e-01  <none>
  //                 a1    0.0000e+00    7.2745e-02 +/-  9.45e-01  <none>
  //                 a2    0.0000e+00    2.5971e-01 +/-  3.27e-01  <none>
  //               mean    1.9685e+00    1.9694e+00 +/-  1.34e-03  <none>
  //              sigma    6.0000e-03    8.1815e-03 +/-  1.83e-03  <none>


  // 3. cn_dsii==3
  // RooFitResult: minimized FCN value: -40738.4, estimated distance to minimum: 6.77349e-06
  //               covariance matrix quality: Full, accurate covariance matrix
  //               Status : MINIMIZE=0 HESSE=0
  //
  //   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  // --------------------  ------------  --------------------------  --------
  //                 Nbkg    5.0000e+05    4.6566e+02 +/-  3.32e+01  <none>
  //                 Nsig    5.0000e+04    3.1873e+03 +/-  6.19e+01  <none>
  //                   a0    0.0000e+00   -4.5356e-01 +/-  5.78e-02  <none>
  //                   a1    0.0000e+00   -8.1463e-01 +/-  7.35e-02  <none>
  //                   a2    0.0000e+00    4.3973e-01 +/-  6.42e-02  <none>
  //                 mean    1.9685e+00    1.9688e+00 +/-  6.56e-05  <none>
  //                sigma    6.0000e-03    3.3107e-03 +/-  5.96e-05  <none>



}
