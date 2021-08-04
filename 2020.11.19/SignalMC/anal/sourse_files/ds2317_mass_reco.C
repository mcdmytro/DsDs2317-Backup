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

void ds2317_mass_reco(){
gROOT->Reset();
gSystem->Load("libRooFit");

Float_t range0=2.2, range1=2.45;
Float_t x0=2.25, x1=2.4;
RooRealVar *x = new RooRealVar("x","",range0,range1);
//x->setBins((int) (range1-range0)/0.0006);
TChain chain("h1");
chain.Add("../../latest.root");

TCanvas *cv = new TCanvas("cv", "",13,123,700,500);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(0);
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

Float_t m_ds17, mctr;
chain.SetBranchAddress("m_ds17", &m_ds17);
chain.SetBranchAddress("mctr", &mctr);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*x),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(m_ds17<range1 && m_ds17>range0 && mctr==1){
    x->setVal(m_ds17);
    data->add(RooArgSet(*x));
      }}

RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",2.317,2.3,2.35,"GeV");
RooRealVar *sigma=new RooRealVar("sigma", "Sigma of gaussian",0.005,0.003,0.007,"GeV");
RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *x, *mean, *sigma);
RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);

RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *x, *mean, *sigma);
RooPolynomial *pol = new RooPolynomial("pol","Polynomial for background",*x, RooArgList(*a0,*a1)) ;

RooRealVar *sig = new RooRealVar("Nsig", "",0.,1E4);
RooRealVar *bkg = new RooRealVar("Nbkg", "",0.,1E4);

RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian*2 + Pol",RooArgList(*gaus, *pol), RooArgList(*sig, *bkg));

x->setRange("signal",x0,x1);
RooFitResult *fitresult = pdf->fitTo(*data,Save(kTRUE),Range("signal")) ;
RooPlot *frame = x->frame();
data->plotOn(frame);

Int_t ci;   // for color index setting
ci = TColor::GetColor("#000099");
frame->SetLineColor(ci);
frame->SetMarkerSize(0.8);
frame->GetXaxis()->SetTitle("M(D_{s}#pi^{0}), (GeV/c^{2})");
frame->GetXaxis()->SetNdivisions(507);
frame->GetXaxis()->SetLabelFont(132);
frame->GetXaxis()->SetLabelSize(0.05);
frame->GetXaxis()->SetTitleSize(0.07);
frame->GetXaxis()->SetTitleFont(132);
frame->GetYaxis()->SetTitle("Events/2.5 MeV/c^{2}");
frame->GetYaxis()->SetNdivisions(507);
frame->GetYaxis()->SetLabelFont(132);
frame->GetYaxis()->SetLabelSize(0.05);
frame->GetYaxis()->SetTitleSize(0.07);
frame->GetYaxis()->SetTitleOffset(0.9);
frame->GetYaxis()->SetTitleFont(132);
//frame->GetZaxis()->SetLabelFont(132);
//frame->GetZaxis()->SetLabelSize(0.05);
//frame->GetZaxis()->SetTitleSize(0.07);
//frame->GetZaxis()->SetTitleFont(132);

pdf->plotOn(frame);
pdf->plotOn(frame,Components(RooArgSet(*gaus)),LineColor(kYellow-2),LineStyle(kDashed));
pdf->plotOn(frame,Components(RooArgSet(*pol)),LineColor(kRed),LineStyle(kDashed));
frame->Draw();

cv->SaveAs("../pdf_files/ds2317_mass_reco.pdf","pdf");
cv->SaveAs("../eps_files/ds2317_mass_reco.eps","eps");
}
