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

void dsds17_pst_overlap(){
gROOT->Reset();
gSystem->Load("libRooFit");
 Float_t range0=0.0, range1=5.0;
  //Float_t x0=0.25, x1=0.475;
RooRealVar *pst = new RooRealVar("pst","",range0,range1);
TChain chain("h1");
chain.Add("../Ds2317_DspNew.root");

 Float_t pst_dsp, mctr, cos0, r2, hel, pst_dsi, m_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,  pst_dsii, m_dsii, mc_dsii, dgf_dsii, c2_dsii, cn_dsii, p_pi0, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, p_ipi, c2_iipi, dgf_iipi, p_iipi, c2_kst0, dgf_kst0, p_kst0, k1_id, k2_id, k3_id, k4_id, pi1_id, pi2_id, mc_ds17, mc_dsi, cosgg, hel_phi;

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

RooDataSet *data = new RooDataSet("data","", RooArgSet(*pst),"GeV");
RooDataSet *seldata = new RooDataSet("seldata","", RooArgSet(*pst),"GeV");

 for (int i=0; i<chain.GetEntries();i++){
   chain.GetEntry(i);
   if(pst_dsp<range1 && pst_dsp>range0 && (TMath::Abs(mc_dsii)!=431 || TMath::Abs(mc_ds17)!=10431) && TMath::Abs(hel_phi)>0.35 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii!=2){
     pst->setVal(pst_dsp);
     data->add(RooArgSet(*pst));
   }}
 for (int i=0; i<chain.GetEntries();i++){
   chain.GetEntry(i);
   if(pst_dsp<range1 && pst_dsp>range0 && (TMath::Abs(mc_dsii)!=431 || TMath::Abs(mc_ds17)!=10431) && pst_ds17>3.5 && TMath::Abs(hel_phi)>0.35 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii!=2){
     pst->setVal(pst_dsp);
     seldata->add(RooArgSet(*pst));
   }}


  pst->setBins(100);
  RooDataHist *dh = new RooDataHist("dh", "binned version of data", RooArgSet(*pst), *data);
  RooDataHist *seldh = new RooDataHist("seldh", "binned version of seldata", RooArgSet(*pst), *seldata);
  /*
  RooRealVar *mean = new RooRealVar("mean", "Mean value",0.3461,0.34,0.36,"GeV");
  RooRealVar *mean2 = new RooRealVar("mean BS", "Mean value",0.3461,0.34,0.36,"GeV");
  RooRealVar *sigma_gaus=new RooRealVar("sigma_gaus", "Sigma of gaussian",0.006,0.003,0.25,"GeV");
  RooRealVar *sigma_gaus2=new RooRealVar("sigma_BS", "Sigma of gaussian",0.006,0.003,0.25,"GeV");

  RooGaussian *gaus= new RooGaussian("gaus", "Gaussian PDF", *mass, *mean, *sigma_gaus);
  RooGaussian *gaus2= new RooGaussian("gaus BS", "Gaussian PDF", *mass, *mean2, *sigma_gaus2);

  RooRealVar *a0 = new RooRealVar("a0","a0",1.,-10.,10.);
  RooRealVar *a1 = new RooRealVar("a1","a1",1.,-10.,10.);
  RooRealVar *a2 = new RooRealVar("a2","a2",1.,-10.,10.);
  RooRealVar *a3 = new RooRealVar("a3","a3",1.,-10.,10.);
  RooRealVar *a4 = new RooRealVar("a4","a4",1.,-10.,10.);
  RooPolynomial *pol = new RooPolynomial("pol","Polynomial PDF",*mass, RooArgList(*a0,*a1,*a3)) ;
  
  RooRealVar *sig = new RooRealVar("Nsig", "", 0, 1E5);
  RooRealVar *BS = new RooRealVar("BS", "", 500, 1E5);
  RooRealVar *bkg = new RooRealVar("Nbkg", "", 0, 1E6);

  RooAddPdf *pdf = new RooAddPdf ("pdf", "Sygnal + background PDF",RooArgList(*gaus ,gaus2, *pol), RooArgList(*sig, *BS, *bkg));
  //RooAddPdf *pdf = new RooAddPdf ("pdf", "Sygnal + background PDF",RooArgList(*gaus, *pol), RooArgList(*sig, *bkg));
  */
  TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);

  // Define "signal" range in x as [x0,x1]
  //mass->setRange("signal",x0,x1);  
  // Fit pdf only to data in "signal" range
  //RooFitResult *fitresult = pdf->fitTo(*datah,Save(kTRUE),Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = pst->frame(Bins(100),Title(" "));
  //data->plotOn(frame);
  dh->plotOn(frame);
  
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.12) ;
  frame->GetXaxis()->SetTitle("p^{*}(D_{s}D_{s0}(2317)), GeV/c ");
  //frame->GetYaxis()->SetTitle("Entries");
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetTitleOffset(1.8);
  //xframe->GetXaxis()->SetTitleOffset(0.8);
  // frame->SetFillColor(kYellow-8);
  
  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  //pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));
  //pdf->paramOn(frame, Format("NEU", 1) ,Layout(0.64,0.99,0.99));
  //frame->getAttText()->SetTextSize(0.02);
  //frame->getAttLine()->SetLineWidth(0);
  
  TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  txt->SetBorderSize(0);
  txt->SetFillColor(0);
  txt->SetTextSize(0.05);
  txt->SetTextFont(80);
  txt->AddText("(a)") ;
  frame->addObject(txt);
  
  //pdf->plotOn(frame);
  //pdf->plotOn(frame,Components(RooArgSet(*gaus)),LineColor(kYellow-2),LineStyle(kDashed));
  //pdf->plotOn(frame,Components(RooArgSet(*gaus2)),LineColor(kGreen+2),LineStyle(kDashed));
  //pdf->plotOn(frame,Components(RooArgSet(*pol)),LineColor(kRed),LineStyle(kDashed));

  seldh->plotOn(frame,LineColor(kGreen+3),MarkerColor(kGreen+3));
  frame->Draw();

  cv->SaveAs("CMC_DsDs2317_pst_overlaped_NOTtrue.pdf","pdf");

}
