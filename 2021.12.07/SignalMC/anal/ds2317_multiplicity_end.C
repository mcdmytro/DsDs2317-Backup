#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
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

void ds2317_multiplicity_end(){
gROOT->Reset();
gSystem->Load("libRooFit");

Float_t range0=-1, range1=21;
Float_t x0=0.1, x1=20.1;
RooRealVar *x = new RooRealVar("x","",range0,range1);
TChain chain("h2");
chain.Add("../Ds2317inDs2317_multiplicity.root");

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

Float_t mult_e;
chain.SetBranchAddress("mult_e", &mult_e);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*x));

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(mult_e>x0 && mult_e<x1){
    x->setVal(mult_e);
    data->add(RooArgSet(*x));
  }}

RooPlot *frame = x->frame(Bins(20));
data->plotOn(frame, DrawOption("B"), DataError(RooAbsData::None), XErrorSize(0), FillColor(kOrange));

Int_t ci;   // for color index setting
ci = TColor::GetColor("#000099");
frame->SetLineColor(ci);
frame->SetMarkerSize(0.8);
frame->GetXaxis()->SetTitle("Evt. multiplicity");
frame->GetXaxis()->SetNdivisions(507);
frame->GetXaxis()->SetLabelFont(132);
frame->GetXaxis()->SetLabelSize(0.05);
frame->GetXaxis()->SetTitleSize(0.07);
frame->GetXaxis()->SetTitleFont(132);
frame->GetYaxis()->SetTitle("Entries");
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

frame->Draw("");

cv->SaveAs("./pdf_files/SMC_ds2317_multiplicity_end.pdf","pdf");
}
