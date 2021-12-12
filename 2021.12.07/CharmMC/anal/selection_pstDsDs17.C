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
#include <vector>

using namespace RooFit;

void selection_pstDsDs17(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  gROOT->SetStyle("Plain");
  gStyle->SetLabelSize(0.025, "XY");

  TCanvas *cv = new TCanvas("cv","Just Canvas",5,5,800,800);

  TChain chain("h1");
  chain.Add("../Ds2317.root");

  std::vector<Float_t> pst_v, a_v;

  Float_t
    r0=0.0,
    r1=5.0,
    pst_max=0,
    a_max=0;

  for(Float_t pst=r0; pst<=r1; pst+=0.01){
    int s=chain.GetEntries(Form("pst_dsp>%f && cn_dsi!=2 && cn_dsii!=2 && TMath::Abs(mc_ds17)==10431", pst));
    int b=chain.GetEntries(Form("pst_dsp>%f && cn_dsi!=2 && cn_dsii!=2", pst));
    float a=s/TMath::Sqrt(b);
    pst_v.push_back(pst);
    a_v.push_back(a);
    //cout << pst << '\t' << s << '\t' << b << '\t' << a << '\n';

    if(a>a_max){
      pst_max=pst;
      a_max=a;
    }
  }

  cout << "pst_max = " << pst_max << '\t' << "a_max = " << a_max << '\n';

  TGraph *g = new TGraph(pst_v.size(), &pst_v[0], &a_v[0]);
  g->SetTitle("");
  g->GetXaxis()->SetTitle("p*, GeV/c ");
  g->GetYaxis()->SetTitle("S/#sqrt{S+B}");
  g->GetXaxis()->SetTitleSize(0.05);
  g->GetYaxis()->SetTitleSize(0.05);
  g->GetXaxis()->SetLabelSize(0.05);
  g->GetYaxis()->SetLabelSize(0.05);
  g->GetXaxis()->SetNdivisions(505);
  g->GetYaxis()->SetNdivisions(505);
  g->GetXaxis()->SetTitleOffset(0.95);
  g->GetYaxis()->SetTitleOffset(0.95);
  g->SetMarkerStyle(34);
  g->SetMarkerSize(1);
  g->Draw("AP");

  TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  txt->SetBorderSize(0);
  txt->SetFillColor(0);
  txt->SetTextSize(0.05);
  txt->SetTextFont(80);
  txt->AddText("(d)") ;
  txt->Draw("SAME");

  cv->SaveAs("CMC_selection_pstDsDs17.pdf","pdf");
}
