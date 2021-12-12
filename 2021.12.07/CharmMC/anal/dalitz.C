#include <stdlib.h>
#include <string.h>
#include <vector>
#include "TLorentzVector.h"
void dalitz(){
	gROOT->Reset();
	gROOT->SetStyle("Plain");
	gStyle->SetLabelSize(0.025, "XY");
	gStyle->SetOptStat(0);
	TChain *h1=new TChain("h1");
	h1->Add("../Ds2317_Dsp.root");

	TCanvas *cv = new TCanvas("cv","Just Canvas",5,5,1500,500);
	cv->Divide(3,1);

	Float_t range1min=0., range1max=190.; // Ds^{pr}Ds^{sec} range
	Float_t range2min=0., range2max=190.; // Ds^{pr}#pi0 range
	Float_t range3min=0., range3max=190.; // Ds^{sec}#pi0 range
	Int_t nch1=100; // Ds^{pr}Ds^{sec} nbins
	Int_t nch2=100; // Ds^{pr}#pi0 nbins
	Int_t nch3=100; // Ds^{sec}#pi0 nbins

	TH2F * dal1 = new TH2F("dal1"," ", nch1, range1min, range1max, nch2, range2min, range2max);
	TH2F * dal2 = new TH2F("dal2"," ", nch1, range1min, range1max, nch3, range3min, range3max);
	TH2F * dal3 = new TH2F("dal3"," ", nch2, range2min, range2max, nch3, range3min, range3max);

	dal1->GetXaxis()->SetTitle("M^{2}(D_{s}^{pr}D_{0}^{sec}), GeV^{2}/c^{4}");
	dal1->GetYaxis()->SetTitle("M^{2}(D_{s}^{pr}#pi^{0}), GeV^{2}/c^{4}");
	dal1->GetXaxis()->SetLabelSize(0.02);
	dal1->GetXaxis()->SetTitleSize(0.03);
	dal1->GetYaxis()->SetLabelSize(0.02);
	dal1->GetYaxis()->SetTitleSize(0.03);
	dal1->GetZaxis()->SetLabelSize(0.02);
	dal1->GetZaxis()->SetTitleSize(0.03);
	dal2->GetXaxis()->SetTitle("M^{2}(D_{s}^{pr}D_{0}^{sec}), GeV^{2}/c^{4}");
	dal2->GetYaxis()->SetTitle("M^{2}(D_{s}^{sec}#pi^{0}), GeV^{2}/c^{4}");
	dal2->GetXaxis()->SetLabelSize(0.02);
	dal2->GetXaxis()->SetTitleSize(0.03);
	dal2->GetYaxis()->SetLabelSize(0.02);
	dal2->GetYaxis()->SetTitleSize(0.03);
	dal2->GetZaxis()->SetLabelSize(0.02);
	dal2->GetZaxis()->SetTitleSize(0.03);
	dal3->GetXaxis()->SetTitle("M^{2}(D_{s}^{pr}#pi^{0}), GeV^{2}/c^{4}");
	dal3->GetYaxis()->SetTitle("M^{2}(D_{s}^{sec}#pi^{0}), GeV^{2}/c^{4}");
	dal3->GetXaxis()->SetLabelSize(0.02);
	dal3->GetXaxis()->SetTitleSize(0.03);
	dal3->GetYaxis()->SetLabelSize(0.02);
	dal3->GetYaxis()->SetTitleSize(0.03);
	dal3->GetZaxis()->SetLabelSize(0.02);
	dal3->GetZaxis()->SetTitleSize(0.03);

	h1->SetAlias("ds1ds2","(e_dsi+e_dsii)*(e_dsi+e_dsii)-(px_dsi+px_dsii)*(px_dsi+px_dsii)-(py_dsi+py_dsii)*(py_dsi+py_dsii)-(pz_dsi+pz_dsii)*(pz_dsi+pz_dsii)");
	h1->SetAlias("ds1pi0","(e_dsi+e_pi0)*(e_dsi+e_pi0)-(px_dsi+px_pi0)*(px_dsi+px_pi0)-(py_dsi+py_pi0)*(py_dsi+py_pi0)-(pz_dsi+pz_pi0)*(pz_dsi+pz_pi0)");
	h1->SetAlias("ds2pi0","(e_dsii+e_pi0)*(e_dsii+e_pi0)-(px_dsii+px_pi0)*(px_dsii+px_pi0)-(py_dsii+py_pi0)*(py_dsii+py_pi0)-(pz_dsii+pz_pi0)*(pz_dsii+pz_pi0)");

	cv->cd(1);
	h1->Draw("ds1ds2","pst_ds17>3.5 && abs(hel_phi)>0.35 && cn_dsi!=2 && cn_dsii!=2");
	//h1->Draw("ds1ds2:ds1pi0>>dal1","pst_ds17>3.5 && abs(hel_phi)>0.35 && cn_dsi!=2 && cn_dsii!=2","COLZ");
	cv->cd(2);
	h1->Draw("ds1pi0","pst_ds17>3.5 && abs(hel_phi)>0.35 && cn_dsi!=2 && cn_dsii!=2");
	//h1->Draw("ds1ds2:ds2pi0>>dal2","pst_ds17>3.5 && abs(hel_phi)>0.35 && cn_dsi!=2 && cn_dsii!=2","COLZ");
	cv->cd(3);
	h1->Draw("ds2pi0","pst_ds17>3.5 && abs(hel_phi)>0.35 && cn_dsi!=2 && cn_dsii!=2");
	//h1->Draw("ds1pi0:ds2pi0>>dal3","pst_ds17>3.5 && abs(hel_phi)>0.35 && cn_dsi!=2 && cn_dsii!=2","COLZ");

	cv->SaveAs("CMC_Dalitz.pdf","pdf");
}
