{
  gROOT->Reset();
  gROOT->SetStyle("Plain");
  gStyle->SetLabelSize(0.025, "XY");
  gStyle->SetOptStat(0);
  //gPad->SetBottomMargin(big);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.05);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadTopMargin(0.05);

  TCanvas *cv = new TCanvas("cv","Just Canvas",5,5,600,600);

  TChain * h1=new TChain("h1");
  h1->Add("../Ds2317_PID.root");

  Float_t range0=0, range1=1;
  Int_t nch=50;

  TH1F * h = new TH1F("h"," ", nch, range0, range1);
  h->SetLineColor(4);
  h->GetXaxis()->SetTitle("P(K/#pi)");
  h->GetYaxis()->SetTitle("Entries");
  h->SetMinimum(0.);
  h->SetTitle("");
  h->SetFillColor(kRed);

  h->SetMarkerSize(0.8);
  h->GetXaxis()->SetNdivisions(507);
  h->GetXaxis()->SetLabelFont(132);
  h->GetXaxis()->SetLabelSize(0.05);
  h->GetXaxis()->SetTitleSize(0.07);
  h->GetXaxis()->SetTitleFont(132);
  h->GetYaxis()->SetNdivisions(507);
  h->GetYaxis()->SetLabelFont(132);
  h->GetYaxis()->SetLabelSize(0.05);
  h->GetYaxis()->SetTitleSize(0.07);
  h->GetYaxis()->SetTitleOffset(0.9);
  h->GetYaxis()->SetTitleFont(132);

  h1->Draw("k1_id>>h","TMath::Abs(mc_k1)==321 && TMath::Abs(mc_dsii)==431 && TMath::Abs(mc_pi0)==111 && pst_ds17>3.5 && TMath::Abs(hel_phi)>0.35 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii!=2",""); // pi+:211, K+:321
  //h1->Draw("k2_id>>+h","TMath::Abs(mc_pi1)==211 && TMath::Abs(mc_dsii)==431 && TMath::Abs(mc_pi0)==111 && pst_ds17>3.5 && TMath::Abs(hel_phi)>0.35 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsii==3","");
  TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  txt->SetBorderSize(0);
  txt->SetFillColor(0);
  txt->SetTextSize(0.05);
  txt->SetTextFont(80);
  txt->AddText("(c)") ;
  txt->Draw("SAME");
  cv->SaveAs("CMC_TrueKaonID_inDs2317_trueDaughters.pdf","pdf");
}


