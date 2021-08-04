{ gROOT->Reset();

  gStyle->SetOptFit(kTRUE);
  Color_t white=10;
  gStyle->SetCanvasColor(white);
  gSystem->Load("libRooFit");
  using namespace RooFit;




  //RooRealVar * delte = new RooRealVar("delte","M_{J/#psi #omega} (GeV/c^{2})",3.81,4.1);                                                                                                      
  RooRealVar *mbc = new RooRealVar("mbc","M_{#phi} [GeV/c^{2}]",1.0045,1.0390);
  
  
  TChain chain("ntuple1");
  chain.Add("data*root");

  Float_t x_delte, x_mbc, x_jmass, x_bflagc, x_bflagchi, x_mctruth,x_phimass;
  
  chain.SetBranchAddress("delte", &x_delte);
  chain.SetBranchAddress("mbc", &x_mbc);
  chain.SetBranchAddress("bflagchi",&x_bflagchi);
  chain.SetBranchAddress("phimass",&x_phimass);

  RooDataSet *data = new RooDataSet("data","", RooArgSet(*mbc));
  
  for (int i=0; i< chain.GetEntries();i++){
    chain.GetEntry(i);
    if(x_bflagchi>0&&x_delte>-0.03&&x_delte<0.03){
      mbc.setVal(x_phimass);
      data.add(RooArgSet(*mbc));
    }}
  RooRealVar *mean1 = new RooRealVar("mean1", "MEAN of gaussian",1.018,1.023);
  RooRealVar *sigma1=new RooRealVar("sigma1", "Sigma of gaussian",0.,0.05);//0.05
  RooRealVar *width1=new RooRealVar("width1", "width of gaussian",0.,0.1);//0.05                                                                       
  RooVoigtian *voig=new RooVoigtian("voigt","voigtian pdf",*mbc,*mean1,*width1,*sigma1);

  RooRealVar *slope1= new RooRealVar("slope1", "Slope of Polynomial", -5.,5.);
  RooRealVar *slope2= new RooRealVar("slope2", "Slope 2  of Polynomial", -1.,1.);
  RooChebychev *chebpol = new RooChebychev("chebpol","Chebychev Polynomial ", *mbc, RooArgList(*slope1));//,*slope2));
  
  
  RooRealVar *sig = new RooRealVar("Nsig", "",0.,400000.);
  RooRealVar *BKG = new RooRealVar("Nbkg", "",0.,500000.);
  
  RooAddPdf *depdf = new RooAddPdf ("depdf", "double Gauss+ Chebychev",RooArgList(*voig,*chebpol),RooArgList(*sig,*BKG));
  
  TCanvas *c1 =new TCanvas ("c1", "mbc Plot", 700, 600);
  RooFitResult *fitresult = depdf.fitTo(*data);//,Extended(true),Minos(false));
  RooPlot *mplot = mbc.frame(Bins(100));
  data.plotOn(mplot);
  depdf.paramOn(mplot,Layout(0.1,0.2,0.9));
  depdf.plotOn(mplot);
  depdf.plotOn(mplot,Components(RooArgSet(*voig)),LineColor(kGreen),LineStyle(kDashed));
  //depdf.plotOn(mplot,Components(RooArgSet(*gauss2)),LineColor(kCyan),LineStyle(kDashed));
  depdf.plotOn(mplot,Components(RooArgSet(*chebpol)),LineColor(kMagenta),LineStyle(kDashed));
  
  RooHist *hpull = mplot->pullHist();
  RooPlot *frame1 = mbc.frame(Title("Pull"));
  frame1->addPlotable(hpull,"P");
  TPad *pad1 = new TPad("pad1","pad1",0.,0.3,1.,1.);
  pad1->Draw();pad1->cd();
  mplot.SetTitle("");
  mplot->GetXaxis()->SetLabelSize(0.05);mplot->GetYaxis()->SetLabelSize(0.05);
  mplot->GetXaxis()->SetTitleSize(0.05);mplot->GetYaxis()->SetTitleSize(0.05);
  mplot->Draw();

  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2",0.0,0.05,1.,0.3);
  pad2->Draw();pad2->cd();
  frame1->GetXaxis()->SetLabelSize(0.125);frame1->GetYaxis()->SetLabelSize(0.125);
  frame1->GetXaxis()->SetTitleSize(0.125);frame1->GetYaxis()->SetTitleSize(0.125);
  frame1->Draw();



  Double_t chi2=mplot.chiSquare();
  // Int_t gdf=(fitresult->floatParsFinal())->getSize();                                                                                                                                        
  std::cout<<"chi^2 = "<<chi2<<std::endl;
  //std::cout<<"GDF = "<<gdf<<std::endl;                                                                                                                                                        
  //std::cout<<"result = "<<chi2/gdf<<std::endl;                                                                                                                                                
  c1->Print("phimass_datahadronbjchg_voigfit_7.eps","eps");
  return c1;

}
