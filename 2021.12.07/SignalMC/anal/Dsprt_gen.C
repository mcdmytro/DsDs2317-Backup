#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooBreitWigner.h"
#include "RooGaussian.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TText.h"
#include "TArrow.h"
#include "TFile.h"

using namespace RooFit;

void Dsprt_gen(){
gROOT->Reset();
gSystem->Load("libRooFit");
Float_t range0=4.286, range1=4.29;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
chain.Add("../Ds2317inDs2317_X4272-like.root");

Float_t dsp_gm;
chain.SetBranchAddress("dsp_gm", &dsp_gm);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass));

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(dsp_gm>range0 && dsp_gm<range1){ 
    mass->setVal(dsp_gm);  
    data->add(RooArgSet(*mass));}
 }

  TCanvas *cv = new TCanvas("cv","Just Canvas",5,5,800,800);

  RooPlot *frame = mass->frame();
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.12) ;
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("M^{gen}(D_{s}D_{s0}^{*}(2317)), GeV/c^{2} ");
  //frame->GetYaxis()->SetTitle("Entries");
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetTitleOffset(1.8);
  //frame->GetXaxis()->SetTitleOffset(0.8);
  // frame->SetFillColor(kYellow-8);

  data->plotOn(frame);
  frame->Draw();

  cv->SaveAs("SMC_Dsprt_gen.pdf","pdf");

}
