#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "RooAddPdf.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"

using namespace RooFit;

void pi0_mass_reco(){
gROOT->Reset();
gSystem->Load("libRooFit");
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);

Float_t range0=0.122, range1=0.148;
Float_t x0=0.122, x1=0.148;
RooRealVar *x = new RooRealVar("x","",range0,range1);
//x->setBins((int) (range1-range0)/0.0006);
TChain chain("h1");
chain.Add("../Ds2317inDs2317_BCS_noDsPi0KmVFit.root");

Float_t m_pi0, mctr;
chain.SetBranchAddress("m_pi0", &m_pi0);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*x),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(m_pi0<range1 && m_pi0>range0){
    x->setVal(m_pi0);
    data->add(RooArgSet(*x));
    }}

RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",0.135,0.122,0.148,"GeV");
RooRealVar *sigma=new RooRealVar("sigma", "Sigma of gaussian",0.0035,0.001,0.008,"GeV");
RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
RooRealVar *a2 = new RooRealVar("a2","a2",0.,-1.,1.);
RooRealVar *a3 = new RooRealVar("a3","a3",0.,-1.,1.);

RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *x, *mean, *sigma);
RooPolynomial *pol = new RooPolynomial("pol","Polynomial for background",*x, RooArgList(*a0,*a1));
RooChebychev *cpol = new RooChebychev("cpol","Chebyshev polynomial for background",*x,RooArgList(*a0));

RooRealVar *sig = new RooRealVar("Nsig", "",0.,1E5);
RooRealVar *bkg = new RooRealVar("Nbkg", "",0.,1E5);

RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus, *cpol), RooArgList(*sig, *bkg));

x->setRange("signal",x0,x1);
RooFitResult *fitresult = pdf->fitTo(*data,Save(kTRUE),Range("signal")) ;
RooPlot *frame = x->frame();
data->plotOn(frame);

Int_t ci;   // for color index setting
ci = TColor::GetColor("#000099");
frame->SetLineColor(ci);
frame->SetMarkerSize(0.8);
frame->GetXaxis()->SetTitle("M(#gamma#gamma) [GeV/c^{2}]");
frame->GetXaxis()->SetNdivisions(507);
frame->GetXaxis()->SetLabelFont(132);
frame->GetXaxis()->SetLabelSize(0.05);
frame->GetXaxis()->SetTitleSize(0.07);
frame->GetXaxis()->SetTitleFont(132);
frame->GetXaxis()->SetTitleOffset(0.9);
frame->GetYaxis()->SetTitle("Events / 0.26 keV/c^{2}");
frame->GetYaxis()->SetNdivisions(507);
frame->GetYaxis()->SetLabelFont(132);
frame->GetYaxis()->SetLabelSize(0.05);
frame->GetYaxis()->SetTitleSize(0.07);
frame->GetYaxis()->SetTitleOffset(1.0);
frame->GetYaxis()->SetTitleFont(132);

pdf->plotOn(frame);
pdf->plotOn(frame,Components(RooArgSet(*gaus)),LineColor(kYellow-2),LineStyle(kDashed));
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

cv->SaveAs("./pdf_files/SMC_BCS_pi0_mass_reco.pdf","pdf");
fitresult->Print("v");

// Fitresult aсquired on 23.07.2021, added in the BNv4

// RooFitResult: minimized FCN value: -150153, estimated distance to minimum: 1.17826e-06
//               covariance matrix quality: Full, accurate covariance matrix
//               Status : MINIMIZE=0 HESSE=0
//
//   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
// --------------------  ------------  --------------------------  --------
//                 Nbkg    5.0000e+04    6.3460e+03 +/-  1.93e+02  <none>
//                 Nsig    5.0000e+04    5.9911e+03 +/-  1.92e+02  <none>
//                   a0    0.0000e+00   -1.0378e-01 +/-  2.77e-02  <none>
//                 mean    1.3500e-01    1.3531e-01 +/-  9.26e-05  <none>
//                sigma    3.5000e-03    3.9097e-03 +/-  1.11e-04  <none>

}
