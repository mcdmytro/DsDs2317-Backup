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

using namespace RooFit;

void SB_Ds(){
gROOT->Reset();
gSystem->Load("libRooFit");
gStyle->SetCanvasPreferGL(1);
Float_t range0=1.91, range1=2.025;
//Float_t range0=2.12, range1=2.5;
Float_t x0=1.91, x1=2.025;
//Float_t x0=2.2, x1=2.5;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
chain.Add("../Ds2317_noDsKmvFit.root");

Float_t mctruth, pst_dsi, m_dsi, mc_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,pst_dsii, m_dsii, dgf_dsii, c2_dsii, cn_dsii, ppi, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, pipi, c2_iipi, dgf_iipi, piipi, mc_ds17, mc_dsi, hel_phi;

chain.SetBranchAddress("m_dsi", &m_dsi);
chain.SetBranchAddress("mc_dsi", &mc_dsi);
chain.SetBranchAddress("pst_dsi", &pst_dsi);
chain.SetBranchAddress("dgf_dsi", &dgf_dsi);
chain.SetBranchAddress("c2_dsi", &c2_dsi);
chain.SetBranchAddress("cn_dsi", &cn_dsi);
chain.SetBranchAddress("m_dsii", &m_dsii);
chain.SetBranchAddress("pst_dsii", &pst_dsii);
chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
chain.SetBranchAddress("c2_dsii", &c2_dsii);
chain.SetBranchAddress("cn_dsii", &cn_dsii);
chain.SetBranchAddress("m_ds17", &m_ds17);
chain.SetBranchAddress("pst_ds17", &pst_ds17);
chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
chain.SetBranchAddress("c2_ds17", &c2_ds17);
chain.SetBranchAddress("p_pi0", &ppi);
chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
chain.SetBranchAddress("c2_pi0", &c2_pi0);
chain.SetBranchAddress("dgf_phi1", &dgf_phi1);
chain.SetBranchAddress("c2_phi1", &c2_phi1);
chain.SetBranchAddress("dgf_phi2", &dgf_phi2);
chain.SetBranchAddress("c2_phi2", &c2_phi2);
chain.SetBranchAddress("c2_ipi", &c2_ipi);
chain.SetBranchAddress("dgf_ipi", &dgf_ipi);
chain.SetBranchAddress("p_ipi", &pipi);
chain.SetBranchAddress("c2_iipi", &c2_iipi);
chain.SetBranchAddress("dgf_iipi", &dgf_iipi);
chain.SetBranchAddress("p_iipi", &piipi);
chain.SetBranchAddress("mctr", &mctruth);
chain.SetBranchAddress("mc_ds17", &mc_ds17);
chain.SetBranchAddress("mc_dsi", &mc_dsi);
chain.SetBranchAddress("hel_phi", &hel_phi);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(m_dsii>range0 && m_dsii<range1 && pst_ds17>3.5 && TMath::Abs(hel_phi)>0.35 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii!=2){// && TMath::Abs(mc_dsi)==431){
    mass->setVal(m_dsii);
    data->add(RooArgSet(*mass));
  }}

  RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",1.9685,1.9,2.0,"GeV");
  RooRealVar *sigma=new RooRealVar("sigma", "Sigma of gaussian",0.006,0.003,0.015,"GeV");
  RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *mass, *mean, *sigma);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1)) ;

  RooRealVar *sig = new RooRealVar("Nsig", "",0.,1E5.);
  RooRealVar *bkg = new RooRealVar("Nbkg", "",0.,1E6.);

  RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian + Pol",RooArgList(*gaus, *Pol), RooArgList(*sig, *bkg));

  TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);
  // Define "signal" range in x as [x0,x1]
  mass.setRange("signal",x0,x1) ;  
  // Fit pdf only to data in "signal" range
  RooFitResult *fitresult = pdf->fitTo(*data,Save(kTRUE),Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = mass->frame(Bins(100));
  data->plotOn(frame);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.12) ;
  frame->SetTitle("");
  frame->GetXaxis()->SetTitle("Mass, GeV ");
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

  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  //pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));
 
  
  pdf->paramOn(frame, Format("NEU", 1) ,Layout(0.64,0.99,0.99));
  frame->getAttText()->SetTextSize(0.025);
  frame->getAttLine()->SetLineWidth(0);
  

  
  Double_t 
    sgnl = 1.948,
    sgnr = 1.982,
    lsbl = 1.92,
    lsbr = 1.937,
    rsbl = 1.993,
    rsbr = 2.01,
    ampl = 3500;

  TLine *vl_sgn_left = new TLine(sgnl, 0, sgnl, ampl);
  vl_sgn_left->SetLineColor(kRed);
  frame->addObject(vl_sgn_left);
  TLine *vl_sgn_right = new TLine(sgnr, 0, sgnr, ampl);
  vl_sgn_right->SetLineColor(kRed);
  frame->addObject(vl_sgn_right);
  TBox *box_sgn = new TBox(sgnl, 0, sgnr, ampl);
  box_sgn->SetFillColorAlpha(kRed-4, 0.15);
  frame->addObject(box_sgn);

  TLine *vl_lsb_left = new TLine(lsbl, 0, lsbl, ampl);
  vl_lsb_left->SetLineColor(kBlue);
  frame->addObject(vl_lsb_left);
  TLine *vl_lsb_right = new TLine(lsbr,0, lsbr, ampl);
  vl_lsb_right->SetLineColor(kBlue);
  frame->addObject(vl_lsb_right);
  TBox *box_lsb = new TBox(lsbl, 0, lsbr, ampl);
  box_lsb->SetFillColorAlpha(kBlue-4, 0.15);
  frame->addObject(box_lsb);

  TLine *vl_rsb_left = new TLine(rsbl, 0, rsbl, ampl);
  vl_rsb_left->SetLineColor(kBlue);
  frame->addObject(vl_rsb_left);
  TLine *vl_rsb_right = new TLine(rsbr, 0, rsbr, ampl);
  vl_rsb_right->SetLineColor(kBlue);
  frame->addObject(vl_rsb_right);
  TBox *box_rsb = new TBox(rsbl, 0, rsbr, ampl);
  box_rsb->SetFillColorAlpha(kBlue-4, 0.15);
  frame->addObject(box_rsb);
  
  /*
  TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  txt->SetBorderSize(0);
  txt->SetFillColor(0);
  txt->SetTextSize(0.05);
  txt->SetTextFont(80);
  txt->AddText("(a)") ;
  frame->addObject(txt); 
  */
  pdf->plotOn(frame);
  pdf->plotOn(frame,Components(RooArgSet(*gaus)),LineColor(kYellow-2),LineStyle(kDashed));
  pdf->plotOn(frame,Components(RooArgSet(*Pol)),LineColor(kRed),LineStyle(kDashed));
  frame->Draw();
  
  mass->setRange("pureSGN", sgnl, sgnr);
  mass->setRange("lsb", lsbl, lsbr);
  mass->setRange("rsb", rsbl, rsbr);

  RooAbsReal* int_sig = Pol->createIntegral(*mass,Range("pureSGN"));
  RooAbsReal* int_lsb = Pol->createIntegral(*mass,Range("lsb"));
  RooAbsReal* int_rsb = Pol->createIntegral(*mass,Range("rsb"));

  cout << "Signal region integral: " << '\t' << int_sig->getVal() << '\n';
  cout <<"Sideband region integral: " << '\t' << int_lsb->getVal()+int_rsb->getVal() << '\n';
  
  //std::cout<<"chi^2 = "<< frame->chiSquare() <<std::endl;
  cv->SaveAs("CMC_sb_Ds.pdf","pdf");
}
