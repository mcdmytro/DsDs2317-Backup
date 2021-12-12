using namespace RooFit;

void phi_mass_reco(){
  gROOT->Reset();
  gSystem->Load("libRooFit");
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  //Set range
  Float_t
    range0=1.009461,
    range1=1.029461,
    x0=1.009461,
    x1=1.029461;

  //Add a root file
  TChain chain("h1");
  chain.Add("../Ds2317inDs2317_BCS_noDsKmvFit.root");

  Float_t m_phi1, mctr, cos0, r2, hel, pst_dsi, m_dsi, dgf_dsi, c2_dsi, cn_dsi, pst_ds17, m_ds17, dgf_ds17, c2_ds17,  pst_dsii, m_dsii, mc_dsii, dgf_dsii, c2_dsii, cn_dsii, p_pi0, dgf_pi0, c2_pi0, dgf_phi1, c2_phi1, dgf_phi2, c2_phi2, c2_ipi, dgf_ipi, p_ipi, c2_iipi, dgf_iipi, p_iipi, k1_id, k2_id, k3_id, k4_id, pi1_id, pi2_id, mc_ds17, mc_dsi, cosgg, hel_phi;

  chain.SetBranchAddress("m_phi1", &m_phi1);
  chain.SetBranchAddress("m_dsi", &m_dsi);
  chain.SetBranchAddress("pst_dsi", &pst_dsi);
  chain.SetBranchAddress("dgf_dsi", &dgf_dsi);
  chain.SetBranchAddress("c2_dsi", &c2_dsi);
  chain.SetBranchAddress("cn_dsi", &cn_dsi);
  chain.SetBranchAddress("m_dsii", &m_dsii);
  chain.SetBranchAddress("mc_dsii", &mc_dsii);
  chain.SetBranchAddress("pst_dsii", &pst_dsii);
  chain.SetBranchAddress("dgf_dsii", &dgf_dsii);
  chain.SetBranchAddress("c2_dsii", &c2_dsii);
  chain.SetBranchAddress("cn_dsii", &cn_dsii);
  chain.SetBranchAddress("m_ds17", &m_ds17);
  chain.SetBranchAddress("pst_ds17", &pst_ds17);
  chain.SetBranchAddress("dgf_ds17", &dgf_ds17);
  chain.SetBranchAddress("c2_ds17", &c2_ds17);
  chain.SetBranchAddress("p_pi0", &p_pi0);
  chain.SetBranchAddress("dgf_pi0", &dgf_pi0);
  chain.SetBranchAddress("c2_pi0", &c2_pi0);
  chain.SetBranchAddress("dgf_phi1", &dgf_phi1);
  chain.SetBranchAddress("c2_phi1", &c2_phi1);
  chain.SetBranchAddress("dgf_phi2", &dgf_phi2);
  chain.SetBranchAddress("c2_phi2", &c2_phi2);
  chain.SetBranchAddress("c2_ipi", &c2_ipi);
  chain.SetBranchAddress("dgf_ipi", &dgf_ipi);
  chain.SetBranchAddress("p_ipi", &p_ipi);
  chain.SetBranchAddress("c2_iipi", &c2_iipi);
  chain.SetBranchAddress("dgf_iipi", &dgf_iipi);
  chain.SetBranchAddress("p_iipi", &p_iipi);
  chain.SetBranchAddress("mctr", &mctr);
  chain.SetBranchAddress("cos0", &cos0);
  chain.SetBranchAddress("r2", &r2);
  chain.SetBranchAddress("mc_ds17", &mc_ds17);
  chain.SetBranchAddress("mc_dsi", &mc_dsi);
  chain.SetBranchAddress("cosgg", &cosgg);
  chain.SetBranchAddress("hel_phi", &hel_phi);

  //Create a dataset and fill it with the data from a root file
  RooRealVar x("x","x variable", range0, range1);
  RooDataSet data("data","dataset", RooArgSet(x));

  for (int i=0; i<chain.GetEntries(); i++){
    chain.GetEntry(i);
    if(m_phi1<range1 && m_phi1>range0 && TMath::Abs(hel_phi)>0.44 && TMath::Prob(c2_ds17, dgf_ds17)>0.001 && TMath::Prob(c2_dsii, dgf_dsii)>0.001 && TMath::Prob(c2_pi0, dgf_pi0)>0.01 && cn_dsi!=2 && cn_dsii!=2){
      x.setVal(m_phi1);
      data.add(RooArgSet(x));
    }}

  //Perform fit
  RooRealVar
    mean("mean", "Mean of gaussian", 1.01, 1, 1.03),
    sigma_g("sigma", "Sigma of gaussian", 0.004, 0.0001, 0.04),
    sigma_bw("sigma_bw","Sigma of a Breit-Wigner ", 0.001, 0.0001, 0.04),
    a0("a0","a0 parameter of a polynomial", 0., -100., 100.),
    a1("a1","a1 parameter of a polynomial", 0., -100., 100.),
    a2("a2","a2 parameter of a polynomial", 0., -100., 100.),
    a3("a3","a3 parameter of a polynomial", 0., -100., 100.);

  RooVoigtian voigt("voigt","Reconstructed peak", x, mean, sigma_bw, sigma_g);
  RooGaussian gaus("gaus", "Gaussian for signal", x, mean, sigma_g);
  RooPolynomial pol("pol","Polynomial for background", x, RooArgList(a0, a1));
  RooChebychev cpol("cpol","Chebyshev polynomial for background",x,RooArgList(a0,a1,a2));

  RooRealVar
    sig("sig", "Number of signal events", 5E3, 5E2, 1E4),
    bkg("bkg", "Number of background events",1E2,1E2,1E3);

  RooAddPdf pdf("pdf", "Sum of signal and backgroung", RooArgList(voigt, cpol), RooArgList(sig, bkg));

  x.setRange("fit range", x0, x1);
  RooFitResult *fitresult = pdf.fitTo(data, Extended(kTRUE), Save(), Range("fit range"));

  //Create and customise frame
  RooPlot *frame = x.frame();
  data.plotOn(frame);
  frame->SetLineColor(TColor::GetColor("#000099"));
  frame->SetMarkerSize(0.8);
  frame->GetXaxis()->SetTitle("M(KK) [GeV/c^{2}]");
  frame->GetXaxis()->SetNdivisions(507);
  frame->GetXaxis()->SetLabelFont(132);
  frame->GetXaxis()->SetLabelSize(0.05);
  frame->GetXaxis()->SetTitleSize(0.07);
  frame->GetXaxis()->SetTitleFont(132);
  frame->GetXaxis()->SetTitleOffset(0.9);
  frame->GetYaxis()->SetTitle("Events / 0.2 MeV/c^{2}");
  frame->GetYaxis()->SetNdivisions(507);
  frame->GetYaxis()->SetLabelFont(132);
  frame->GetYaxis()->SetLabelSize(0.05);
  frame->GetYaxis()->SetTitleSize(0.07);
  frame->GetYaxis()->SetTitleOffset(1.0);
  frame->GetYaxis()->SetTitleFont(132);

  pdf.plotOn(frame);
  pdf.plotOn(frame, Components(RooArgSet(voigt)), LineColor(kYellow-2), LineStyle(kDashed));
  pdf.plotOn(frame, Components(RooArgSet(cpol)), LineColor(kRed), LineStyle(kDashed));

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

  cv->SaveAs("./pdf_files/SMC_BCS_phi_mass_reco.pdf","pdf");
  fitresult->Print("v");

  // Fitresult aсquired on 23.07.2021, added in the BNv4

  // RooFitResult: minimized FCN value: -54930.7, estimated distance to minimum: 1.31603
  //               covariance matrix quality: Full, accurate covariance matrix
  //               Status : MINIMIZE=-1 HESSE=4
  //
  //   Floating Parameter  InitialValue    FinalValue +/-  Error     GblCorr.
  // --------------------  ------------  --------------------------  --------
  //                   a0    0.0000e+00    1.2327e+00 +/-  5.75e-01  <none>
  //                   a1    0.0000e+00    6.8030e-01 +/-  6.47e-04  <none>
  //                   a2    0.0000e+00    4.4578e-01 +/-  7.04e-01  <none>
  //                  bkg    1.0000e+02    1.0001e+02 +/-  1.28e+02  <none>
  //                 mean    1.0100e+00    1.0194e+00 +/-  6.16e-05  <none>
  //                  sig    5.0000e+03    4.5812e+03 +/-  7.19e+01  <none>
  //                sigma    4.0000e-03    9.5259e-04 +/-  1.91e-04  <none>
  //             sigma_bw    1.0000e-03    4.0991e-03 +/-  2.44e-04  <none>

}
