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

void dsds2317_xP_reco(){
gROOT->Reset();
gSystem->Load("libRooFit");

Float_t range0=0.0, range1=1.0;
//Float_t x0=1.935, x1=2.01;
RooRealVar *x = new RooRealVar("x","",range0,range1);
TChain chain("h1");
chain.Add("../../Ds2317inDs2317_DspNew.root");

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

 Float_t x_p, pst_dsp, mctr, cos0, r2, hel, pst_dsi, m_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,  pst_dsii, m_dsii, mc_dsii, dgf_dsii, c2_dsii, cn_dsii, p_pi0, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, p_ipi, c2_iipi, dgf_iipi, p_iipi, c2_kst0, dgf_kst0, p_kst0, k1_id, k2_id, k3_id, k4_id, pi1_id, pi2_id, mc_ds17, mc_dsi, cosgg, hel_phi;

 chain.SetBranchAddress("x_p", &x_p);
chain.SetBranchAddress("m_dsi", &m_dsi);
chain.SetBranchAddress("pst_dsi", &pst_dsi);
chain.SetBranchAddress("dgf_dsi", &dgf_dsi);
chain.SetBranchAddress("c2_dsi", &c2_dsi);
chain.SetBranchAddress("cn_dsi", &cn_dsi);
chain.SetBranchAddress("m_dsii", &m_dsii);
chain.SetBranchAddress("mc_dsii", &mc_dsii);
chain.SetBranchAddress("pst_dsii", &pst_dsii);
chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
chain.SetBranchAddress("c2_dsii", &c2_dsii);
chain.SetBranchAddress("cn_dsii", &cn_dsii);
chain.SetBranchAddress("m_ds17", &m_ds17);
chain.SetBranchAddress("pst_ds17", &pst_ds17);
chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
chain.SetBranchAddress("c2_ds17", &c2_ds17);
chain.SetBranchAddress("p_pi0", &p_pi0);
chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
chain.SetBranchAddress("c2_pi0", &c2_pi0);
chain.SetBranchAddress("dgf_phi1", &dgf_phi1);
chain.SetBranchAddress("c2_phi1", &c2_phi1);
chain.SetBranchAddress("dgf_phi2", &dgf_phi2);
chain.SetBranchAddress("c2_phi2", &c2_phi2);
chain.SetBranchAddress("c2_ipi", &c2_ipi);
chain.SetBranchAddress("dgf_ipi", &dgf_ipi);
chain.SetBranchAddress("p_ipi", &p_ipi);
chain.SetBranchAddress("c2_iipi", &c2_iipi);
chain.SetBranchAddress("dgf_iipi", &dgf_iipi);
chain.SetBranchAddress("p_iipi", &p_iipi);
chain.SetBranchAddress("c2_kst0", &c2_kst0);
chain.SetBranchAddress("dgf_kst0", &dgf_kst0);
chain.SetBranchAddress("p_kst0", &p_kst0);
chain.SetBranchAddress("mctr", &mctr);
chain.SetBranchAddress("cos0", &cos0);
chain.SetBranchAddress("r2", &r2);
chain.SetBranchAddress("mc_ds17", &mc_ds17);
chain.SetBranchAddress("mc_dsi", &mc_dsi);
chain.SetBranchAddress("cosgg", &cosgg);
chain.SetBranchAddress("hel_phi", &hel_phi);
chain.SetBranchAddress("pst_dsp", &pst_dsp);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*x),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(x_p<range1 && x_p>range0 && TMath::Abs(hel_phi)>0.35 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii!=2){
    x->setVal(x_p);
    data->add(RooArgSet(*x));
  }}

RooPlot *frame = x->frame(Bins(100));
data->plotOn(frame);

Int_t ci;   // for color index setting
ci = TColor::GetColor("#000099");
frame->SetLineColor(ci);
frame->SetMarkerSize(0.8);
frame->GetXaxis()->SetTitle("x_{P}");
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
TArrow *arr = new TArrow(3.5,2.5E5,3.5,0,0.03,"|>");
arr->SetAngle(40);
arr->SetLineWidth(2);
frame->addObject(arr);
*/
frame->Draw();

cv->SaveAs("../pdf_files/SMC_dsds2317_xP.pdf","pdf");
cv->SaveAs("../eps_files/SMC_dsds2317_xP.eps","eps");
}
