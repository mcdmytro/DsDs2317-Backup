//***************************************************
//
//                     13.01.2020
//                   Dmytro Meleshko
//            Monte-Carlo simulation analisys
//                      Dalitz plot
//
//***************************************************

{
gROOT->Reset();
gROOT->SetStyle("Plain");
gStyle->SetLabelSize(0.025, "XY");
TChain *h1=new TChain("h1");
h1->Add("../latest.root");

TCanvas *cv1 = new TCanvas("cv1","Just Canvas",5,5,800,800);

Float_t rangeXmin=0.8, rangeXmax=2.0;
Float_t rangeYmin=0.5, rangeYmax=1.6;
Int_t nchX=200;
Int_t nchY=200;

TH2F * Dalitz = new TH2F("Dalitz","Dalitz plot for KK - K#pi", nchX, rangeXmin, rangeXmax, nchY, rangeYmin, rangeYmax);
Dalitz->GetXaxis()->SetTitle("KK mass, GeV");
Dalitz->GetYaxis()->SetTitle("K#pi mass, GeV");
Dalitz->GetXaxis()->SetLabelSize(0.02);
Dalitz->GetXaxis()->SetTitleSize(0.03);
Dalitz->GetYaxis()->SetLabelSize(0.02);
Dalitz->GetYaxis()->SetTitleSize(0.03);
Dalitz->GetZaxis()->SetLabelSize(0.02);
Dalitz->GetZaxis()->SetTitleSize(0.03);

h1->SetAlias("kk","sqrt((e_k3+e_k4)**2-(px_k3+px_k4)**2-(py_k3+py_k4)**2-(pz_k3+pz_k4)**2)");
h1->SetAlias("kpi","sqrt((e_k3+e_pi1)**2-(px_k3+px_pi1)**2-(py_k3+py_pi1)**2-(pz_k3+pz_pi1)**2)");
h1->Draw("kpi:kk>>Dalitz","pst_dsi>2.5 && TMath::Prob(c2_dsi, dgf_dsi)>0.001","COLZ");

TLine *line1 = new TLine(1.01946,0.5,1.01946,1.6);
TLine *line2 = new TLine(0.8,0.896,2.0,0.896);
line1->SetLineWidth(2);
line2->SetLineWidth(2);
line1->SetLineColor(kRed);
line2->SetLineColor(kRed);
line1->Draw();
line2->Draw();

TLatex l1, l2;
l1.SetTextSize(0.03);
l2.SetTextSize(0.03);
l1.DrawLatex(1.05,0.55, "#color[2]{#varphi(1020)}");
l1.DrawLatex(0.85,0.85, "#color[2]{#varphi(896)}");

cv1->SaveAs("CorrelationPlot.pdf","pdf");
}
