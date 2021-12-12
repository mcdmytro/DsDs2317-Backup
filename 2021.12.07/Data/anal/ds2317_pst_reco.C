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

void ds2317_pst_reco(){
gROOT->Reset();
gSystem->Load("libRooFit");

Float_t range0=0.0, range1=5.0;
//Float_t x0=1.935, x1=2.01;
RooRealVar *x = new RooRealVar("x","",range0,range1);
TChain chain("h1");
chain.Add("../Ds2317.root");

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

 Float_t m_ds17, cn_dsii, mctr, hel_phi, hel_ds17, pst_ds17, c2_ds17, dgf_ds17, c2_dsii, dgf_dsii, c2_pi0, dgf_pi0, pst_dsii, m_dsii, hel_ds17;
 chain.SetBranchAddress("m_ds17", &m_ds17);
 chain.SetBranchAddress("m_dsii", &m_dsii);
 chain.SetBranchAddress("cn_dsii", &cn_dsii);
 chain.SetBranchAddress("mctr", &mctr);
 chain.SetBranchAddress("hel_phi", &hel_phi);
 chain.SetBranchAddress("hel_ds17", &hel_ds17);
 chain.SetBranchAddress("pst_ds17", &pst_ds17);
 chain.SetBranchAddress("c2_ds17", &c2_ds17);
 chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
 chain.SetBranchAddress("c2_dsii", &c2_dsii);
 chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
 chain.SetBranchAddress("c2_pi0", &c2_pi0);
 chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
 chain.SetBranchAddress("pst_dsii", &pst_dsii);
 
RooDataSet *data = new RooDataSet("data","", RooArgSet(*x),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(pst_ds17<range1 && pst_ds17>range0 && TMath::Abs(hel_phi)>0.35 && cn_dsii!=2 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01){
    x->setVal(pst_ds17);
    data->add(RooArgSet(*x));
  }}

RooPlot *frame = x->frame(Bins(100));
data->plotOn(frame);

Int_t ci;   // for color index setting
ci = TColor::GetColor("#000099");
frame->SetLineColor(ci);
frame->SetMarkerSize(0.8);
frame->GetXaxis()->SetTitle("p*(D_{s}#pi^{0}), (GeV/c)");
frame->GetXaxis()->SetNdivisions(507);
frame->GetXaxis()->SetLabelFont(132);
frame->GetXaxis()->SetLabelSize(0.05);
frame->GetXaxis()->SetTitleSize(0.07);
frame->GetXaxis()->SetTitleFont(132);
frame->GetYaxis()->SetTitle("Events/5.0 MeV/c");
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

TArrow *arr = new TArrow(3.5,2.5E5,3.5,0,0.03,"|>");
arr->SetAngle(40);
arr->SetLineWidth(2);
frame->addObject(arr);

frame->Draw();

cv->SaveAs("Data_ds2317_pst.pdf","pdf");
//cv->SaveAs("../eps_files/ds2317_pst_reco.eps","eps");
}
