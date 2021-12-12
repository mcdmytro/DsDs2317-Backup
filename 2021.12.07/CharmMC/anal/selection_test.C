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

void selection_test(){
  gROOT->Reset();

  TChain chain("h1");
  chain.Add("../Ds2317.root");

  std::vector<Float_t>
    xp_v, xp_a_v;

  Float_t
    xp_0=0.0,
    xp_1=0.45,
    xp_max=0,
    xp_aMax=0;

  for(Float_t xp=xp_0; xp_1<=xp_1; xp+=0.005){
    int s=chain.GetEntries(Form("x_p>%f && cn_dsi!=2 && cn_dsii!=2 && TMath::Abs(mc_ds17)==10431", xp));
    int b=chain.GetEntries(Form("x_p>%f && cn_dsi!=2 && cn_dsii!=2 ", xp));
    xp_v.push_back(xp);
    xp_a_v.push_back(s/TMath::Sqrt(b));
    if(TMath::Sqrt(b)>xp_aMax){
      xp_max=xp;
      xp_aMax=s/TMath::Sqrt(b);
    }
  }

  cout << "xp_max = " << xp_max << '\t' << "xp_aMax = " << xp_aMax << '\n';
  cout << "xp_max = " << xp_v[std::max_element(xp_a_v.begin(),xp_a_v.end()) - xp_a_v.begin()] << '\t' << "a_max = " << xp_a_v[std::max_element(xp_a_v.begin(),xp_a_v.end()) - xp_a_v.begin()] << '\n';
}
