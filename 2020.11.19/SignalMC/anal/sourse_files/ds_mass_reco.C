#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
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

void ds_mass_reco(){
gROOT->Reset();
gSystem->Load("libRooFit");

Float_t range0=1.94, range1=2.0;
Float_t x0=1.94, x1=2.00;
RooRealVar *x = new RooRealVar("x","",range0,range1);
//x->setBins((int) (range1-range0)/0.0006);
TChain chain("h1");
chain.Add("../../latest.root");

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

int entr=0;
Float_t m_dsi, mctr;
chain.SetBranchAddress("m_dsi", &m_dsi);
chain.SetBranchAddress("mctr", &mctr);

RooDataSet *data = new RooDataSet("data","", RooArgSet(*x),"GeV");

for (int i=0; i<chain.GetEntries();i++){
  chain.GetEntry(i);
  if(m_dsi<range1 && m_dsi>range0 && mctr==1){
    x->setVal(m_dsi);
    data->add(RooArgSet(*x));
    if(m_dsi<x1 && m_dsi>x0) entr++;
  }}

RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",1.9685,1.9,2.0,"GeV");
RooRealVar *sigma=new RooRealVar("sigma", "Sigma of gaussian",0.0035,0.002,0.004,"GeV");
RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *x, *mean, *sigma);
RooRealVar *sigma1=new RooRealVar("sigma1", "Sigma of gaussian",0.0068,0.006,0.008,"GeV");
RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);

RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *x, *mean, *sigma);
RooGaussian *gaus1= new RooGaussian("gaus1", "Gaussian for signal", *x, *mean, *sigma1);
RooPolynomial *pol = new RooPolynomial("pol","Polynomial for background",*x, RooArgList(*a0,*a1)) ;

RooRealVar *sig = new RooRealVar("Nsig", "",0.,1E5);
RooRealVar *bkg = new RooRealVar("Nbkg", "",0.,1E5);

RooRealVar *sigsum = new RooRealVar ("sigsum","signak",0.85,0.,1.);
RooAddPdf *sign = new RooAddPdf ("sign","Signal",RooArgList(*gaus,*gaus1), *sigsum);
RooAddPdf *pdf = new RooAddPdf ("pdf", "Gaussian*2 + Pol",RooArgList(*sign, *pol), RooArgList(*sig, *bkg));

x->setRange("signal",x0,x1);
RooFitResult *fitresult = pdf->fitTo(*data,Save(kTRUE),Range("signal")) ;
RooPlot *frame = x->frame();
data->plotOn(frame);

Int_t ci;   // for color index setting
ci = TColor::GetColor("#000099");
frame->SetLineColor(ci);
frame->SetMarkerSize(0.8);
frame->GetXaxis()->SetTitle("M(D_{s}), (GeV/c^{2})");
frame->GetXaxis()->SetNdivisions(507);
frame->GetXaxis()->SetLabelFont(132);
frame->GetXaxis()->SetLabelSize(0.05);
frame->GetXaxis()->SetTitleSize(0.07);
frame->GetXaxis()->SetTitleFont(132);
frame->GetYaxis()->SetTitle("Events/0.6 MeV/c^{2}");
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

pdf->plotOn(frame);
pdf->plotOn(frame,Components(RooArgSet(*sign)),LineColor(kYellow-2),LineStyle(kDashed));
pdf->plotOn(frame,Components(RooArgSet(*pol)),LineColor(kRed),LineStyle(kDashed));
frame->Draw();

int nFloatParams = fitresult->floatParsFinal().getSize();
int fcn = fitresult->minNll();
int bins = x->getBins();
int ndf = bins-nFloatParams;

std::cout<<"entr = "<< entr <<std::endl;
std::cout<<"Nbins = "<< bins <<std::endl;
std::cout<<"NFloatParams = "<< nFloatParams <<std::endl;
std::cout<<"ndf = "<< ndf <<std::endl;
std::cout<<"FCN = "<< fcn <<std::endl;
std::cout<<"chi^{2}/dgf = "<< TMath::Prob(fcn, -ndf) <<std::endl;

cv->SaveAs("../pdf_files/ds_mass_reco.pdf","pdf");
cv->SaveAs("../eps_files/ds_mass_reco.eps","eps");
}
