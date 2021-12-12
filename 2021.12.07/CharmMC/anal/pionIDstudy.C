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
  h1->Add("../Ds2317-continuum-PID.root");

  Float_t range0=0, range1=1;
  Int_t nch=50;

  TH1F* hist0 = new TH1F("hist0"," ", nch, range0, range1);
  hist0->SetLineColor(1);
  hist0->GetXaxis()->SetTitle("P(K/#pi)");
  hist0->GetYaxis()->SetTitle("Entries");
  hist0->SetMinimum(0.);
  hist0->SetTitle("");
  hist0->SetFillColor(29);

  hist0->SetMarkerSize(0.8);
  hist0->GetXaxis()->SetNdivisions(507);
  hist0->GetXaxis()->SetLabelFont(132);
  hist0->GetXaxis()->SetLabelSize(0.05);
  hist0->GetXaxis()->SetTitleSize(0.07);
  hist0->GetXaxis()->SetTitleFont(132);
  hist0->GetYaxis()->SetNdivisions(507);
  hist0->GetYaxis()->SetLabelFont(132);
  hist0->GetYaxis()->SetLabelSize(0.05);
  hist0->GetYaxis()->SetTitleSize(0.07);
  hist0->GetYaxis()->SetTitleOffset(0.9);
  hist0->GetYaxis()->SetTitleFont(132);

  TH1F * hist1 = new TH1F("hist1"," ", nch, range0, range1);
  TH1F * hist2 = new TH1F("hist2"," ", nch, range0, range1);
  hist1->SetFillColor(1);
  hist2->SetFillColor(2);

  //k1_id>0.5 && k2_id>0.2 && k3_id>0.5 && k4_id>0.2 && pi1_id<0.9 && pi2_id<0.9 &&
  h1->Draw("pi1_id>>hist0","abs(mc_k1)==321 && abs(mc_k2)==321 && abs(mc_k3)==321 && abs(mc_k4)==321 && abs(mc_pi1)==211 && abs(mc_pi2)==211 && pst_ds17>2.79 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi==1 && cn_dsii==1","");
  h1->Draw("pi1_id>>hist1","k1_id>0.5 && k2_id>0.2 && pi1_id<0.9 && abs(mc_k1)==321 && abs(mc_k2)==321 && abs(mc_k3)==321 && abs(mc_k4)==321 && abs(mc_pi1)==211 && abs(mc_pi2)==211 && pst_ds17>2.79 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi==1 && cn_dsii==1","SAME");
  h1->Draw("pi1_id>>hist2","k1_id>0.5 && k2_id>0.2 && k3_id>0.5 && k4_id>0.2 && pi1_id<0.9 && pi2_id<0.9 && abs(mc_k1)==321 && abs(mc_k2)==321 && abs(mc_k3)==321 && abs(mc_k4)==321 && abs(mc_pi1)==211 && abs(mc_pi2)==211 && pst_ds17>2.79 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi==1 && cn_dsii==1","SAME");
  cv->SaveAs("CMC_pionIDvar_comparison_0.5.0.2.0.9.pdf","pdf");
}
