using namespace RooFit;

void FitData_kst(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  gSystem->Load("libRooFit");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  Float_t range0=0.8, range1=1.0;
  //Float_t range0=2.12, range1=2.5;
  Float_t x0=0.8, x1=1;
  //Float_t x0=2.2, x1=2.5;
  RooRealVar *mass = new RooRealVar("mass","",range0,range1);
  TChain chain("h1");
  chain.Add("../Ds2317inDs2317_BCS_kstCheck.root");

  Float_t mctruth, cos, R2, helicity,pst_dsi, m_dsi, dgf_dsi, chi2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, chi2_ds17,pst_dsii, m_dsii, dgf_dsii, chi2_dsii, cn_dsii, ppi, dgf_pi, chi2_pi,  chi2_ipi, dgf_ipi, pipi, chi2_iipi, dgf_iipi, piipi, k1_id, k2_id, k3_id, k4_id, pi1_id, pi2_id, mc_ds17, mc_dsi, mc_phi1, m_phi1, dgf_phi1, chi2_phi1, mc_phi2, m_phi2, dgf_phi2, chi2_phi2, m_kst1;

  chain.SetBranchAddress("m_dsi", &m_dsi);
  chain.SetBranchAddress("pst_dsi", &pst_dsi);
  chain.SetBranchAddress("dgf_dsi", &dgf_dsi);
  chain.SetBranchAddress("c2_dsi", &chi2_dsi);
  chain.SetBranchAddress("cn_dsi", &cn_dsi);
  chain.SetBranchAddress("m_dsii", &m_dsii);
  chain.SetBranchAddress("pst_dsii", &pst_dsii);
  chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
  chain.SetBranchAddress("c2_dsii", &chi2_dsii);
  chain.SetBranchAddress("cn_dsii", &cn_dsii);
  chain.SetBranchAddress("m_ds17", &m_ds17);
  chain.SetBranchAddress("pst_ds17", &pst_ds17);
  chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
  chain.SetBranchAddress("c2_ds17", &chi2_ds17);
  chain.SetBranchAddress("m_phi1", &m_phi1);
  chain.SetBranchAddress("dgf_phi1", &dgf_phi1);
  chain.SetBranchAddress("c2_phi1", &chi2_phi1);
  chain.SetBranchAddress("m_phi2", &m_phi2);
  chain.SetBranchAddress("dgf_phi2", &dgf_phi2);
  chain.SetBranchAddress("c2_phi2", &chi2_phi2);
  chain.SetBranchAddress("p_pi0", &ppi);
  chain.SetBranchAddress("dgf_pi0", &dgf_pi);
  chain.SetBranchAddress("c2_pi0", &chi2_pi);
  chain.SetBranchAddress("c2_ipi", &chi2_ipi);
  chain.SetBranchAddress("dgf_ipi", &dgf_ipi);
  chain.SetBranchAddress("p_ipi", &pipi);
  chain.SetBranchAddress("c2_iipi", &chi2_iipi);
  chain.SetBranchAddress("dgf_iipi", &dgf_iipi);
  chain.SetBranchAddress("p_iipi", &piipi);
  chain.SetBranchAddress("mctr", &mctruth);
  chain.SetBranchAddress("cos0", &cos);
  chain.SetBranchAddress("r2", &R2);
  //chain.SetBranchAddress("hel", &helicity);
  //chain.SetBranchAddress("k1_id", &k1_id);
  //chain.SetBranchAddress("k2_id", &k2_id);
  //chain.SetBranchAddress("k3_id", &k3_id);
  //chain.SetBranchAddress("k4_id", &k4_id);
  //chain.SetBranchAddress("pi1_id", &pi1_id);
  //chain.SetBranchAddress("pi2_id", &pi2_id);
  chain.SetBranchAddress("mc_ds17", &mc_ds17);
  chain.SetBranchAddress("mc_dsi", &mc_dsi);
  chain.SetBranchAddress("mc_phi2", &mc_phi2);
  chain.SetBranchAddress("m_kst1", &m_kst1);
  RooDataSet *data = new RooDataSet("data","", RooArgSet(*mass),"GeV");

  for (int i=0; i<chain.GetEntries();i++){
    chain.GetEntry(i);
    if(m_kst1>range0 && m_kst1<range1){
      mass->setVal(m_kst1);
      data->add(RooArgSet(*mass));
  }}


  RooRealVar *mean = new RooRealVar("mean", "Mean of gaussian",0.895,0.845,0.945,"GeV");
  RooRealVar *sigma=new RooRealVar("sigma", "Sigma of gaussian",0.05,0.01,0.1,"GeV");
  RooGaussian *gaus= new RooGaussian("gaus", "Gaussian for signal", *mass, *mean, *sigma);
  RooRealVar *sigma_bw = new RooRealVar("sigma_bw","Sigma of Breit Wigner",0.05,0.01,0.1,"GeV");

  RooRealVar *a0 = new RooRealVar("a0","a0",0.,-100.,100.);
  RooRealVar *a1 = new RooRealVar("a1","a1",0.,-100.,100.);
  RooRealVar *a2 = new RooRealVar("a2","a2",0.,-100.,100.);
  RooRealVar *a3 = new RooRealVar("a3","a3",0.,-100.,100.);
  RooPolynomial *Pol = new RooPolynomial("Pol","Polynomial for background",*mass, RooArgList(*a0,*a1)) ;

  RooBreitWigner *bw = new RooBreitWigner("bw","Breit-Wigner function fit", *mass, *mean, *sigma_bw);
  RooVoigtian *voigt = new RooVoigtian("voigt","Reconstructed peak", *mass, *mean, *sigma_bw, *sigma);
  RooChebychev *cpol = new RooChebychev("cpol","Chebyshev polynomial for background", *mass, RooArgList(*a0));

  RooRealVar *sig = new RooRealVar("Nsig", "",0.,1E4);
  RooRealVar *bkg = new RooRealVar("Nbkg", "",0.,1E3);

  RooAddPdf *pdf = new RooAddPdf ("pdf", "Voigtian + Pol",RooArgList(*bw, *cpol), RooArgList(*sig, *bkg));

  // Define "signal" range in x as [x0,x1]
  mass->setRange("signal",x0,x1) ;
  // Fit pdf only to data in "signal" range
  RooFitResult *fitresult = pdf->fitTo(*data,Save(kTRUE),Range("signal")) ;
  //RooFitResult *fitresult = pdf->fitTo(*data);
  RooPlot *frame = mass->frame(Title(" "));
  data->plotOn(frame);
  frame->GetXaxis()->SetTitle("M(K#pi) [GeV/c^{2}]");
  frame->GetXaxis()->SetNdivisions(507);
  frame->GetXaxis()->SetLabelFont(132);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleSize(0.07);
  frame->GetXaxis()->SetTitleFont(132);
  frame->GetXaxis()->SetTitleOffset(0.9);
  frame->GetYaxis()->SetTitle("Events / 2 MeV/c^{2}");
  frame->GetYaxis()->SetNdivisions(507);
  frame->GetYaxis()->SetLabelFont(132);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.07);
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleFont(132);
  /*
  TPaveText* txt = new TPaveText(0.8,0.8,0.85,0.85,"blNDC");
  txt->SetBorderSize(0);
  txt->SetFillColor(0);
  txt->SetTextSize(0.05);
  txt->SetTextFont(80);
  txt->AddText("(b)") ;
  frame->addObject(txt);
  */
  //gaus->paramOn(frame,Layout(0.56,0.89,0.89));
  //pdf->paramOn(frame, Format("NELU", FixedPrecision(2)) ,Layout(0.69,0.99,0.99));
  //pdf->paramOn(frame, Format("NEU", 1) ,Layout(0.64,0.99,0.99));
  //frame->getAttText()->SetTextSize(0.025);
  //frame->getAttLine()->SetLineWidth(0);

  pdf->plotOn(frame);
  pdf->plotOn(frame,Components(RooArgSet(*bw)),LineColor(kYellow-2),LineStyle(kDashed));
  pdf->plotOn(frame,Components(RooArgSet(*cpol)),LineColor(kRed),LineStyle(kDashed));

  TArrow *arr1 = new TArrow(0.842,3.5E4,0.842,0,0.03,"|>");
  arr1->SetAngle(40);
  arr1->SetLineWidth(2);
  frame->addObject(arr1);

  TArrow *arr2 = new TArrow(0.942,3.5E4,0.942,0,0.03,"|>");
  arr2->SetAngle(40);
  arr2->SetLineWidth(2);
  frame->addObject(arr2);

  //Create and customize canvas
  TCanvas *cv=new TCanvas("cv","Just Canvas",5,5,800,800);
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
  gPad->SetLeftMargin(0.15);
  gPad->SetBottomMargin(0.15);

  frame->Draw();

  // //Find the Chi2 value
  // Double_t chi2 = frame->chiSquare();
  //
  // //Find the number of bins of the fit range
  // Int_t lower = frame->GetXaxis()->FindFixBin(x0);
  // Int_t upper = frame->GetXaxis()->FindFixBin(x1);
  // Int_t range = upper-lower;
  //
  // //Find the number of free parameters
  // Int_t dgf=(fitresult->floatParsFinal())->getSize();
  //
  // //Find the Chi2/nDOF value
  // Double_t new_chi2 = chi2/(range - dgf);
  //
  // cout << "\nChi2:" << chi2;
  // cout << "\nNumber of bins:" << range;
  // cout << "\nNumber of free parameters:" << dgf;
  // cout << "\nChi2/nDOF for (Breit Wigner + Chebyshev Polynom) fit:" << new_chi2 << endl;

  cv->SaveAs("./pdf_files/SMC_BCS_kst.pdf","pdf");
  fitresult->Print("v");

  // Fitresult aсquired on 23.07.2021, added in the BNv4

//   RooFitResult: minimized FCN value: -61220.1, estimated distance to minimum: 5.80716e-05
//               covariance matrix quality: Full, accurate covariance matrix
//               Status : MINIMIZE=0 HESSE=0
//
//   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
// --------------------  ------------  --------------------------  --------
//                 Nbkg    5.0000e+02    9.5741e+02 +/-  1.14e+02  <none>
//                 Nsig    5.0000e+03    5.3918e+03 +/-  1.38e+02  <none>
//                   a0    0.0000e+00    4.2925e-01 +/-  1.13e-01  <none>
//                 mean    8.9500e-01    8.9341e-01 +/-  4.83e-04  <none>
//             sigma_bw    5.0000e-02    4.6385e-02 +/-  1.62e-03  <none>

}
