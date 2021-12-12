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

void residual(){
gROOT->Reset();
gSystem->Load("libRooFit");

Float_t range0=-0.01, range1=0.01;
//Float_t range0=2.12, range1=2.5;
Float_t x0=-0.01, x1=0.01;
//Float_t x0=2.2, x1=2.5;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
chain.Add("../latest.root");

Float_t m_phi2, phi_gen;

chain.SetBranchAddress("m_phi2", &m_phi2);
chain.SetBranchAddress("phi_gen", &phi_gen);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if((m_phi2-phi_gen)>range0 && (m_phi2-phi_gen)<range1 && m_phi2>0 && phi_gen>0){
    mass->setVal(m_phi2-phi_gen);
    data->add(RooArgSet(*mass));
}}


  RooRealVar *mean = new RooRealVar("mean", "Mean value",0,-0.01,0.01,"GeV");
  RooRealVar *sigma_gaus = new RooRealVar("sigma_gaus", "Sigma of gaussian",0.00052,0.0003,0.015,"GeV");
  RooRealVar *sigma_bw = new RooRealVar("sigma_bw","Sigma of Breit Wigner",0.0006,0.0003,0.002,"GeV");  

  RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *mass, *mean, *sigma_gaus);
  RooBreitWigner *bw = new RooBreitWigner("bw","Breit-Wigner function fit", *mass, *mean, *sigma_bw);
  RooVoigtian *voigt = new RooVoigtian("voigt","Reconstructed peak", *mass, *mean, *sigma_bw, *sigma_gaus);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-1.,1.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-1.,1.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-1.,1.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-1.,1.);

  RooPolynomial *pol = new RooPolynomial("pol","Polynomial for background",*mass, RooArgList(*a0,*a1));
  RooChebychev *cpol = new RooChebychev("cpol","Chebyshev polynomial for background",*mass,RooArgList(*a0,*a1));

  RooRealVar *sig = new RooRealVar("Nsig", "",0.,1E5);
  RooRealVar *bkg = new RooRealVar("Nbkg", "",0.,1E5);

  RooAddPdf *pdf = new RooAddPdf ("pdf", "Signal + background function",RooArgList(*bw, *cpol), RooArgList(*sig, *bkg));

  TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);
  
  // Define "signal" range in x as [x0,x1]
  mass.setRange("signal",x0,x1) ;
  // Fit pdf only to data in "signal" range
  RooFitResult *fitresult = pdf->fitTo(*data,Save(kTRUE),Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = mass->frame(Title(" "));
  data->plotOn(frame);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.12) ;
  frame->GetXaxis()->SetTitle("M(#varphi^{rec})-M(#varphi^{gen}), GeV ");
  frame->GetXaxis()->SetTitleSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.05);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetNdivisions(505);
  frame->GetYaxis()->SetTitleOffset(1.8);
  //frame->GetXaxis()->SetTitleOffset(0.8);
  // frame->SetFillColor(kYellow-8);

  TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  txt->SetBorderSize(0);
  txt->SetFillColor(0);
  txt->SetTextSize(0.05);
  txt->SetTextFont(80);
  txt->AddText("(b)") ;
  frame->addObject(txt);

  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  //pdf->paramOn(frame, Format("NELU", 1), Layout(0.6,0.9,0.9));
  //frame->getAttText()->SetTextSize(0.025);
  //frame->getAttLine()->SetLineWidth(0);
 
  pdf->plotOn(frame);
  /*
  pdf->plotOn(frame,
	      Components(RooArgSet(*bw)),
	      LineColor(kYellow-2),
	      LineStyle(kDashed));
  */
  pdf->plotOn(frame,
	      Components(RooArgSet(*cpol)),
	      LineColor(kRed),
	      LineStyle(kDashed));
  frame->Draw();

  //Find the Chi2 value
  Double_t chi2 = frame->chiSquare();

  //Find the number of bins of the fit range
  Int_t lower = frame->GetXaxis()->FindFixBin(x0);
  Int_t upper = frame->GetXaxis()->FindFixBin(x1);
  Int_t range = upper-lower;
 
  //Find the number of free parameters
  Int_t dgf=(fitresult->floatParsFinal())->getSize();

  //Find the Chi2/nDOF value
  Double_t new_chi2 = chi2/(range - dgf);

  cout << "\nChi2:" << chi2;
  cout << "\nNumber of bins:" << range;
  cout << "\nNumber of free parameters:" << dgf;
  cout << "\nChi2/nDOF for (Breit Wigner + Chebyshev Polynom) fit:" << new_chi2 << endl;

  RooDataHist* datah = data->binnedClone();
  pdf.chi2FitTo(*datah) ;
  RooChi2Var chi2_lowstat("chi2_lowstat","chi2",*pdf,*datah) ;
  cout << chi2_lowstat.getVal() << endl ;
  cv->SaveAs("residual_ready.pdf","pdf");
}
