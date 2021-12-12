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

void ds2317_r2_reco(){
gROOT->Reset();
gSystem->Load("libRooFit");

Float_t range0=0.0, range1=0.9;
//Float_t x0=1.935, x1=2.01;
RooRealVar *x = new RooRealVar("x","",range0,range1);
TChain chain("h1");
chain.Add("../../Ds2317inDs2317.root");

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

 Float_t mctruth, cos, hel_ds17, hel_phi, pst_dsi, m_dsi, dgf_dsi, chi2_dsi, cn_dsi, pst_ds17, m_ds17, mc_ds17, dgf_ds17, chi2_ds17, cn_ds17, pst_dsii, m_dsii, dgf_dsii, chi2_dsii, cn_dsii, ppi, dgf_pi, chi2_pi, dgf_phi1, chi2_phi1, dgf_phi2, chi2_phi2, chi2_ipi, dgf_ipi, pipi, chi2_iipi, dgf_iipi, piipi, chi2_kst, dgf_kst, pkst, r2;

 chain.SetBranchAddress("m_dsi", &m_dsi);
 chain.SetBranchAddress("pst_dsi", &pst_dsi);
 chain.SetBranchAddress("dgf_dsi", &dgf_dsi);
 chain.SetBranchAddress("c2_dsi", &chi2_dsi);
 chain.SetBranchAddress("cn_dsi", &cn_dsi);
 chain.SetBranchAddress("m_dsii", &m_dsii);
 chain.SetBranchAddress("pst_dsii", &pst_dsii);
 chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
 chain.SetBranchAddress("c2_dsii", &chi2_dsii);
 chain.SetBranchAddress("cn_dsii", &cn_dsii);
 chain.SetBranchAddress("m_ds17", &m_ds17);
 chain.SetBranchAddress("mc_ds17", &mc_ds17);
 chain.SetBranchAddress("pst_ds17", &pst_ds17);
 chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
 chain.SetBranchAddress("c2_ds17", &chi2_ds17);
 chain.SetBranchAddress("cn_ds17", &cn_ds17);
 chain.SetBranchAddress("p_pi0", &ppi);
 chain.SetBranchAddress("dgf_pi0", &dgf_pi);
 chain.SetBranchAddress("c2_pi0", &chi2_pi);
 chain.SetBranchAddress("dgf_phi1", &dgf_phi1);
 chain.SetBranchAddress("c2_phi1", &chi2_phi1);
 chain.SetBranchAddress("dgf_phi2", &dgf_phi2);
 chain.SetBranchAddress("c2_phi2", &chi2_phi2);
 chain.SetBranchAddress("c2_ipi", &chi2_ipi);
 chain.SetBranchAddress("dgf_ipi", &dgf_ipi);
 chain.SetBranchAddress("p_ipi", &pipi);
 chain.SetBranchAddress("c2_iipi", &chi2_iipi);
 chain.SetBranchAddress("dgf_iipi", &dgf_iipi);
 chain.SetBranchAddress("p_iipi", &piipi);
 chain.SetBranchAddress("c2_kst0", &chi2_kst);
 chain.SetBranchAddress("dgf_kst0", &dgf_kst);
 chain.SetBranchAddress("p_kst0", &pkst);
 chain.SetBranchAddress("mctr", &mctruth);
 chain.SetBranchAddress("cos0", &cos);
 chain.SetBranchAddress("r2", &r2);
 chain.SetBranchAddress("hel_ds17", &hel_ds17);
 chain.SetBranchAddress("hel_phi", &hel_phi);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*x),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(r2<range1 && r2>range0 && TMath::Abs(hel_phi)>0.35 && hel_ds17<0.7 && cn_dsii!=2){
    x->setVal(r2);
    data->add(RooArgSet(*x));
  }}

RooPlot *frame = x->frame(Bins(100));
data->plotOn(frame);

Int_t ci;   // for color index setting
ci = TColor::GetColor("#000099");
frame->SetLineColor(ci);
frame->SetMarkerSize(0.8);
frame->GetXaxis()->SetTitle("R2");
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
/*
 TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
 txt->SetBorderSize(0);
 txt->SetFillColor(0);
 txt->SetTextSize(0.05);
 txt->SetTextFont(80);
 txt->AddText("(b)") ;
 frame->addObject(txt); 
*/
frame->Draw();

cv->SaveAs("../pdf_files/SMC_ds2317_r2_reco.pdf","pdf");
cv->SaveAs("../eps_files/SMC_ds2317_r2_reco.eps","eps");
}
