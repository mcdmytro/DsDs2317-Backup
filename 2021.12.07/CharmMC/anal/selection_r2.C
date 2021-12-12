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

void selection_r2(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  gROOT->SetStyle("Plain");
  gStyle->SetLabelSize(0.025, "XY");

  TCanvas *cv = new TCanvas("cv","Just Canvas",5,5,800,800);

  TChain chain("h1");
  chain.Add("../Ds2317.root");

  std::vector<Float_t> R2_v, a_v;

  Float_t
    r0=0.0,
    r1=0.75,
    R2_max=0,
    a_max=0;

  for(Float_t R2=r0; R2<=r1; R2+=0.005){
    int s=chain.GetEntries(Form("r2>%f && cn_dsi!=2 && cn_dsii!=2 && TMath::Abs(mc_ds17)==10431", R2));
    int b=chain.GetEntries(Form("r2>%f && cn_dsi!=2 && cn_dsii!=2", R2));
    float a=s/TMath::Sqrt(b);
    R2_v.push_back(R2);
    a_v.push_back(a);
    //cout << pst << '\t' << s << '\t' << b << '\t' << a << '\n';

    if(a>a_max){
      R2_max=R2;
      a_max=a;
    }
  }

  cout << "R2_max = " << R2_max << '\t' << "a_max = " << a_max << '\n';

  TGraph *g = new TGraph(R2_v.size(), &R2_v[0], &a_v[0]);
  g->SetTitle("");
  g->GetXaxis()->SetTitle("R_{2}");
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
  txt->AddText("(b)") ;
  txt->Draw("SAME");

  cv->SaveAs("CMC_selection_r2.pdf","pdf");
}
