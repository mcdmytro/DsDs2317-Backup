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

void efficiencyStudy(){
gROOT->Reset();
gSystem->Load("libRooFit");
Float_t range0=4.2863, range1=4.29;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
chain.Add("../Ds2317inDs2317_X4272-like.root");

Float_t m_dsp, dsp_gm;
chain.SetBranchAddress("m_dsp", &m_dsp);
chain.SetBranchAddress("dsp_gm", &dsp_gm);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass));
RooDataSet *gend = new RooDataSet("gend","", RooArgSet(*mass));

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(m_dsp>range0 && m_dsp<range1){ 
    mass->setVal(m_dsp);  
    data->add(RooArgSet(*mass));}
  if(dsp_gm>range0 && dsp_gm<range1){
    mass->setVal(dsp_gm);
    gend->add(RooArgSet(*mass));}
 }

  RooDataHist *rdh = new RooDataHist("rdh", "binned version of data", RooArgSet(*mass), *data);
  RooDataHist *rgh = new RooDataHist("rgh", "binned version of gend", RooArgSet(*mass), *gend);
  
  TH1F *dh = rdh->createHistogram("  ", *mass, Binning(100,range0,range1));
  TH1F *gh = rgh->createHistogram("  ", *mass, Binning(100,range0,range1));

  TH1F *eff = dh->Clone();
  eff->Divide(gh);

  TCanvas *cv = new TCanvas("cv","Just Canvas",5,5,800,800);

  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.12) ;
  eff->GetXaxis()->SetTitle("M(D_{s}D_{s0}*(2317)), GeV/c^{2}");
  eff->GetYaxis()->SetTitle("Efficiency");
  eff->GetXaxis()->SetTitleSize(0.05);
  eff->GetYaxis()->SetTitleSize(0.05);
  eff->GetXaxis()->SetLabelSize(0.05);
  eff->GetYaxis()->SetLabelSize(0.05);
  eff->GetXaxis()->SetNdivisions(505);
  eff->GetYaxis()->SetNdivisions(505);
  eff->GetYaxis()->SetTitleOffset(1.8);
  eff->SetTitle("  ");
  eff->SetMarkerStyle(8); 
  eff->SetMarkerSize(1);
  eff->SetMarkerColor(1);
  eff->SetLineColor(1);
  eff->SetLineWidth(1.5);
  eff->SetStats(kFALSE);
  eff->Draw("E1");
  
  cv->SaveAs("SMC_efficiency","pdf");

}
