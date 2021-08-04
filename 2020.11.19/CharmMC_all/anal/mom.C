{
   gROOT->Reset();
   gROOT->SetStyle("Plain");
   gStyle->SetLabelSize(0.025, "XY");
   gStyle->SetOptStat(0);
   TCanvas *cv = new TCanvas("cv","Just Canvas",5,5,800,800);

   TChain *h1=new TChain("h1");
   h1->Add("./../gen.root");

   Float_t range0=0, range1=5;
   Int_t nch=100;

   TH1F *h = new TH1F("h","", nch, range0, range1);
   //h->SetLineColor(4);
   h->GetXaxis()->SetTitle("Mass, GeV ");
   h->GetYaxis()->SetTitle("Entries");
   h->GetXaxis()->SetTitleSize(0.05);
   h->GetYaxis()->SetTitleSize(0.05);
   h->GetXaxis()->SetLabelSize(0.05);
   h->GetYaxis()->SetLabelSize(0.05);
   h->GetXaxis()->SetNdivisions(505);
   h->GetYaxis()->SetNdivisions(505);
   h->GetYaxis()->SetTitleOffset(1.6);

   h1->Draw("p_phi2>>h","p_phi2>0");
   
   TLine *line1 = new TLine(0.7,0,0.7,50000);
   TLine *line2 = new TLine(1.5,0,1.5,50000);
   TLine *line3 = new TLine(2.5,0,2.5,50000);
   TLine *line4 = new TLine(3.5,0,3.5,50000);
   line1->SetLineWidth(2);
   line2->SetLineWidth(2);
   line3->SetLineWidth(2);
   line4->SetLineWidth(2);
   line1->SetLineColor(kRed);
   line2->SetLineColor(kRed);
   line3->SetLineColor(kRed);
   line4->SetLineColor(kRed);
   line1->Draw();
   line2->Draw();
   line3->Draw();
   line4->Draw();

   
 }
