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

void phi_mass_gen(){
gROOT->Reset();
gSystem->Load("libRooFit");

Float_t range0=1.00946, range1=1.02946;
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

Float_t phi_gm;
chain.SetBranchAddress("phi_gm", &phi_gm);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*x),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(phi_gm<range1 && phi_gm>range0){
    x->setVal(phi_gm);
    data->add(RooArgSet(*x));
      }}

 RooRealVar *mean = new RooRealVar("mean", "Mean value",1.02,1.01,1.03,"GeV");
 RooRealVar *sigma_bw = new RooRealVar("sigma_bw","Sigma of Breit Wigner",0.002,0.001,0.006,"GeV");  

 RooBreitWigner *bw = new RooBreitWigner("bw","Breit-Wigner function fit", *x, *mean, *sigma_bw);
 
 RooRealVar *a0 = new RooRealVar("a0","a0",0.,-1.,1.);
 RooRealVar *a1 = new RooRealVar("a1","a1",0.,-1.,1.);
 
 RooPolynomial *pol = new RooPolynomial("pol","Polynomial for background",*x, RooArgList(*a0,*a1));
 
 RooRealVar *sig = new RooRealVar("Nsig", "",0.,5E4);
 RooRealVar *bkg = new RooRealVar("Nbkg", "",0.,1E3);

 RooAddPdf *pdf = new RooAddPdf ("pdf", "Signal + background function",RooArgList(*bw, *pol), RooArgList(*sig, *bkg));

 // Define "signal" range in x as [x0,x1]
 //x->setRange("signal",x0,x1) ;
 // Fit pdf only to data in "signal" range
 RooFitResult *fitresult = pdf->fitTo(*data,Save(kTRUE));//,Range("signal")) ;
 //RooFitResult *fitresult = pdf->fitTo(*data);
 RooPlot *frame = x->frame(Title(" "));
data->plotOn(frame);

Int_t ci;   // for color index setting
ci = TColor::GetColor("#000099");
frame->SetLineColor(ci);
frame->SetMarkerSize(0.8);
frame->GetXaxis()->SetTitle("M(KK), (GeV/c^{2})");
frame->GetXaxis()->SetNdivisions(507);
frame->GetXaxis()->SetLabelFont(132);
frame->GetXaxis()->SetLabelSize(0.05);
frame->GetXaxis()->SetTitleSize(0.07);
frame->GetXaxis()->SetTitleFont(132);
frame->GetYaxis()->SetTitle("Events");
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
  pdf->plotOn(frame,
        Components(RooArgSet(*bw)),
	      LineColor(kYellow-2),
	            LineStyle(kDashed));
 pdf->plotOn(frame,
	     Components(RooArgSet(*pol)),
	     LineColor(kRed),
	     LineStyle(kDashed));
frame->Draw();

cv->SaveAs("../pdf_files/phi_mass_gen.pdf","pdf");
cv->SaveAs("../eps_files/phi_mass_gen.eps","eps");
}
