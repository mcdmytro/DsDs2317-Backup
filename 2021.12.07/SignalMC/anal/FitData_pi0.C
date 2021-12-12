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

void FitData_pi0(){
gROOT->Reset();
gSystem->Load("libRooFit");
Float_t range0=0.122, range1=0.148;
//Float_t range0=2.12, range1=2.5;
Float_t x0=0.122, x1=0.148;
//Float_t x0=2.2, x1=2.5;
RooRealVar *mass = new RooRealVar("mass","",range0,range1);
TChain chain("h1");
chain.Add("../pi0.root");

Float_t m_pi0;
chain.SetBranchAddress("m_pi0", &m_pi0);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(m_pi0>range0 && m_pi0<range1){
    mass->setVal(m_pi0);
    data->add(RooArgSet(*mass));
}}

  RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",0.1349770, 0.13,0.14,"GeV");
  RooRealVar *sigma=new RooRealVar("sigma", "Sigma of gaussian",0.001,0.0005,0.005,"GeV");
  RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *mass, *mean, *sigma);

  RooRealVar *a0 = new RooRealVar("a0","a0",0.5, 0., 1.);
  RooRealVar *a1 = new RooRealVar("a1","a1",1., 0., 10.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1)) ;

  RooRealVar *GZ = new RooRealVar("GZ","Z width", 0.0042, 0.001, 0.01, "GeV");
  RooVoigtian *voigt = new RooVoigtian("voigt","Reconstructed peak", *mass, *mean, *GZ, *sigma);
 
  RooRealVar *sig = new RooRealVar("Nsig", "", 0, 1E5);
  RooRealVar *bkg = new RooRealVar("Nbkg", "", 0, 1E5);

  RooAddPdf *pdf = new RooAddPdf ("pdf", "Voigtian + Pol",RooArgList(*gaus, *Pol), RooArgList(*sig, *bkg));

  TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);
  // Define "signal" range in x as [x0,x1]
  mass.setRange("signal",x0,x1) ;  
  // Fit pdf only to data in "signal" range
  RooFitResult *fitresult = pdf->fitTo(*data,Save(kTRUE),Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = mass->frame();
  data->plotOn(frame);
  gPad->SetLeftMargin(0.2);
  gPad->SetBottomMargin(0.12) ;
  frame->SetTitle(" ");
  frame->GetXaxis()->SetTitle("M(#gamma#gamma), GeV ");
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
  //pdf->paramOn(frame, Format("NEU", 1) ,Layout(0.685,0.99,0.99));
  //frame->getAttText()->SetTextSize(0.025);
  //frame->getAttLine()->SetLineWidth(0);

  pdf->plotOn(frame);
  pdf->plotOn(frame,Components(RooArgSet(*gaus)),LineColor(kYellow-2),LineStyle(kDashed));
  pdf->plotOn(frame,Components(RooArgSet(*Pol)),LineColor(kRed),LineStyle(kDashed));
  frame->Draw();

  // std::cout<<"chi^2 = "<< frame->chiSquare() <<std::endl;
  //Find the Chi2 value
  Double_t chi2 = frame->chiSquare();

  //Find the number of bins of the fit range
  Int_t lower = frame->GetXaxis()->FindFixBin(x0);
  Int_t upper = frame->GetXaxis()->FindFixBin(x1);
  Int_t range = upper-lower;
 
  //Find the number of free parameters
  Int_t dgf=(fitresult.floatParsFinal())->getSize();

  //Find the Chi2/nDOF value
  Double_t new_chi2 = chi2/(range - dgf);

  cout << "\nChi2:" << chi2;
  cout << "\nNumber of bins:" << range;
  cout << "\nNumber of free parameters:" << dgf;
  cout << "\nChi2/nDOF for (Breit Wigner + Chebyshev Polynom) fit:" << new_chi2 << endl;

  cv->SaveAs("SMC_pi0_ready.pdf","pdf");
}
