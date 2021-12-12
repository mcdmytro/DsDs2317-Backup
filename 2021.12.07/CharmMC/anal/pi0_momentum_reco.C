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

void pi0_momentum_reco(){
gROOT->Reset();
gSystem->Load("libRooFit");

Float_t range0=0.1, range1=0.8;
//Float_t x0=1.935, x1=2.01;
RooRealVar *x = new RooRealVar("x","",range0,range1);
TChain chain("h1");
chain.Add("../Ds2317_check.root");

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

Float_t p_pi0;
chain.SetBranchAddress("p_pi0", &p_pi0);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*x),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(p_pi0<range1 && p_pi0>range0){
    x->setVal(p_pi0);
    data->add(RooArgSet(*x));
  }}

RooPlot *frame = x->frame(Bins(100));
data->plotOn(frame);

Int_t ci;   // for color index setting
ci = TColor::GetColor("#000099");
frame->SetLineColor(ci);
frame->SetMarkerSize(0.8);
frame->GetXaxis()->SetTitle("p(#gamma#gamma), (GeV/c)");
frame->GetXaxis()->SetNdivisions(507);
frame->GetXaxis()->SetLabelFont(132);
frame->GetXaxis()->SetLabelSize(0.05);
frame->GetXaxis()->SetTitleSize(0.07);
frame->GetXaxis()->SetTitleFont(132);
frame->GetYaxis()->SetTitle("Events/7 MeV/c");
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

TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
txt->SetBorderSize(0);
txt->SetFillColor(0);
txt->SetTextSize(0.05);
txt->SetTextFont(80);
txt->AddText("(a)") ;
frame->addObject(txt); 

TArrow *arr = new TArrow(0.15,2.5E5,0.15,0,0.03,"|>");
arr->SetAngle(40);
arr->SetLineWidth(2);
frame->addObject(arr);

frame->Draw();

cv->SaveAs("../CMC_pi0_momentum_reco.pdf","pdf");
}
