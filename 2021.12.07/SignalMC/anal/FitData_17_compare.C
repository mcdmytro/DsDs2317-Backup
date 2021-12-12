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
#include "TLegend.h"
#include "TLegendEntry.h"

using namespace RooFit;

void FitData_17_compare(){
gROOT->Reset();
gSystem->Load("libRooFit");
Float_t range0=2.17, range1=2.5;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain1("h1");
TChain chain2("h1");
TChain chain3("h1");
TChain chain4("h1");
TChain chain5("h1");

chain1.Add("../Ds2317inDs2317_noInclusiveness.root");
chain2.Add("../Ds2317inDs2317_withInclusiveness.root");
chain3.Add("../Ds2317inDs2317_direct.root");
chain4.Add("../../CharmMC_all/Ds2317.root");
chain5.Add("../../Data/Ds2317.root");

Float_t 
  m_ds17_1, cn_dsii_1, hel_ds17_1, hel_phi_1,
  m_ds17_2, cn_dsii_2, hel_ds17_2, hel_phi_2,
  m_ds17_3, cn_dsii_3, hel_ds17_3, hel_phi_3,
  m_ds17_4, cn_dsii_4, hel_ds17_4, hel_phi_4, pst_ds17_4,
  m_ds17_5, cn_dsii_5, hel_ds17_5, hel_phi_5, pst_ds17_5;

chain1.SetBranchAddress("m_ds17", &m_ds17_1);
chain1.SetBranchAddress("cn_dsii", &cn_dsii_1);
chain1.SetBranchAddress("hel_ds17", &hel_ds17_1);
chain1.SetBranchAddress("hel_phi", &hel_phi_1);

chain2.SetBranchAddress("m_ds17", &m_ds17_2);
chain2.SetBranchAddress("cn_dsii", &cn_dsii_2);
chain2.SetBranchAddress("hel_ds17", &hel_ds17_2);
chain2.SetBranchAddress("hel_phi", &hel_phi_2);

chain3.SetBranchAddress("m_ds17", &m_ds17_3);
chain3.SetBranchAddress("cn_dsii", &cn_dsii_3);
chain3.SetBranchAddress("hel_ds17", &hel_ds17_3);
chain3.SetBranchAddress("hel_phi", &hel_phi_3);

chain4.SetBranchAddress("m_ds17", &m_ds17_4);
chain4.SetBranchAddress("cn_dsii", &cn_dsii_4);
chain4.SetBranchAddress("hel_ds17", &hel_ds17_4);
chain4.SetBranchAddress("hel_phi", &hel_phi_4);
chain4.SetBranchAddress("pst_ds17", &pst_ds17_4);

chain5.SetBranchAddress("m_ds17", &m_ds17_5);
chain5.SetBranchAddress("cn_dsii", &cn_dsii_5);
chain5.SetBranchAddress("hel_ds17", &hel_ds17_5);
chain5.SetBranchAddress("hel_phi", &hel_phi_5);
chain5.SetBranchAddress("pst_ds17", &pst_ds17_5);

RooDataSet *data1 = new RooDataSet("data1","", RooArgSet(*mass),"GeV");
RooDataSet *data2 = new RooDataSet("data2","", RooArgSet(*mass),"GeV");
RooDataSet *data3 = new RooDataSet("data3","", RooArgSet(*mass),"GeV");
RooDataSet *data4 = new RooDataSet("data4","", RooArgSet(*mass),"GeV");
RooDataSet *data5 = new RooDataSet("data5","", RooArgSet(*mass),"GeV");

for(int i=0; i<chain1.GetEntries(); i++){
  chain1.GetEntry(i);
  if(m_ds17_1>range0 && m_ds17_1<range1 && TMath::Abs(hel_phi_1)>0.35 && hel_ds17_1<0.7 && cn_dsii_1!=2){
    mass->setVal(m_ds17_1);
    data1->add(RooArgSet(*mass));
  }}

for(int i=0; i<chain2.GetEntries(); i++){
  chain2.GetEntry(i);
  if(m_ds17_2>range0 && m_ds17_2<range1 && TMath::Abs(hel_phi_2)>0.35 && hel_ds17_2<0.7 && cn_dsii_2!=2){
    mass->setVal(m_ds17_2);
    data2->add(RooArgSet(*mass));
  }}

for(int i=0; i<chain3.GetEntries(); i++){
  chain3.GetEntry(i);
  if(m_ds17_3>range0 && m_ds17_3<range1 && TMath::Abs(hel_phi_3)>0.35 && hel_ds17_3<0.7 && cn_dsii_3!=2){
    mass->setVal(m_ds17_3);
    data3->add(RooArgSet(*mass));
  }}

for(int i=0; i<chain4.GetEntries(); i++){
  chain4.GetEntry(i);
  if(m_ds17_4>range0 && m_ds17_4<range1 && TMath::Abs(hel_phi_4)>0.35 && cn_dsii_4!=2 && pst_ds17_4>3.5){
    mass->setVal(m_ds17_4);
    data4->add(RooArgSet(*mass));
  }}

for(int i=0; i<chain5.GetEntries(); i++){
  chain5.GetEntry(i);
  if(m_ds17_5>range0 && m_ds17_5<range1 && TMath::Abs(hel_phi_5)>0.35 && cn_dsii_5!=2 && pst_ds17_5>3.5){
    mass->setVal(m_ds17_5);
    data5->add(RooArgSet(*mass));
  }}


 TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);

 RooPlot *frame = mass->frame();
 data1->plotOn(frame, MarkerSize(0.7));
 data2->plotOn(frame, MarkerColor(kGreen+3), MarkerSize(0.7), LineColor(kGreen+3));
 data3->plotOn(frame, MarkerColor(kRed+1), MarkerSize(0.7), LineColor(kRed+1));
 data4->plotOn(frame, MarkerColor(kMagenta+2), MarkerSize(0.7), LineColor(kMagenta+2));
 data5->plotOn(frame, MarkerColor(kBlue+2), MarkerSize(0.7), LineColor(kBlue+2));

 gPad->SetLeftMargin(0.2);
 gPad->SetBottomMargin(0.12) ;
 frame->SetTitle("");
 frame->GetXaxis()->SetTitle("M(D_{s}#pi^{0}), GeV ");
 //frame->GetYaxis()->SetTitle("Entries");
 frame->GetXaxis()->SetTitleSize(0.05);
 frame->GetYaxis()->SetTitleSize(0.05);
 frame->GetXaxis()->SetLabelSize(0.05);
 frame->GetYaxis()->SetLabelSize(0.05);
 frame->GetXaxis()->SetNdivisions(505);
 frame->GetYaxis()->SetNdivisions(505);
 frame->GetYaxis()->SetTitleOffset(1.8);
 //frame->GetXaxis()->SetTitleOffset(0.8);
 //frame->SetFillColor(kYellow-8);

 frame->Draw();

 TLegend *leg = new TLegend(0.55,0.65,0.89,0.89);
 leg->SetFillColor(kWhite);
 leg->SetLineColor(kWhite);
 leg->AddEntry(data1, "BLACK - Dsprt exclusive signal MC", "p");
 leg->AddEntry(data2, "GREEN - Dsprt inclusive signal MC", "p");
 leg->AddEntry(frame->FindObject("data3"), "RED - Direct signal MC", "p");
 leg->AddEntry(data4, "MAGENTA - Generic MC", "p");
 leg->AddEntry(data5, "BLUE - Data", "p");
 //frame->addObject(leg);
 leg->Draw("SAME");

  cv->SaveAs("Ds2317inDs2317_compare.pdf","pdf");
 }
