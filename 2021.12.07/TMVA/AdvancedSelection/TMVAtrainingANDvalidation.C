#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

void TMVAtrainingANDvalidation( ){

  TMVA::Tools::Instance();

  TFile* outputFile = TFile::Open( "./TMVA.root", "RECREATE" );

  TMVA::Factory* factory = new TMVA::Factory( "MVAnalysis", outputFile,
                                            "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification");

  TFile* input = TFile::Open("./CMC_Ds2317_forTMVA.root");

  TMVA::DataLoader* dataloader = new TMVA::DataLoader("dataset");

  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;

  dataloader->    AddSignalTree((TTree*)input->Get("TreeS"),      signalWeight);
  dataloader->AddBackgroundTree((TTree*)input->Get("TreeB"), backgroundWeight );

  dataloader->AddVariable("pst_ds17",  'F');
  dataloader->AddVariable("c2_ds17",   'F');
  dataloader->AddVariable("pc2_ds17",  'F');
  dataloader->AddVariable("pst_dsii",  'F');
  dataloader->AddVariable("c2_dsii",   'F');
  dataloader->AddVariable("pc2_dsii",  'F');
  dataloader->AddVariable("pc2_pi0",   'F');
  dataloader->AddVariable("p_pi0",     'F');
  dataloader->AddVariable("c2_pi0",    'F');
  dataloader->AddVariable("egam1",     'F');
  dataloader->AddVariable("egam2",     'F');
  dataloader->AddVariable("hel_phi",   'F');

  dataloader->AddSpectator( "dm_ds17", 'F');
  dataloader->AddSpectator( "mc_ds17", 'F');

  TCut mycuts = "";
  TCut mycutb = "";

  dataloader->PrepareTrainingAndTestTree(mycuts, mycutb,
                                        "nTrain_Signal=2500:nTrain_Background=2500:SplitMode=Random:NormMode=NumEvents:!V");

  factory->BookMethod(dataloader, TMVA::Types::kMLP, "MLP",
                     "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=8:TestRate=5:!UseRegulator");

  factory->BookMethod(dataloader, TMVA::Types::kBDT, "BDT",
                     "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
  factory->BookMethod(dataloader, TMVA::Types::kFisher, "Fisher",
                     "H:!V:Fisher:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10");
  factory->BookMethod(dataloader, TMVA::Types::kLikelihood, "LikelihoodPCA",
                     "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA");
  /*
  // General layout.
  TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");

  // Training strategies.
  TString training0("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                   "ConvergenceSteps=30,BatchSize=256,TestRepetitions=10,"
                   "WeightDecay=1e-4,Regularization=None,"
                   "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
  TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                   "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                   "WeightDecay=1e-4,Regularization=L2,"
                   "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
  TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                   "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                   "WeightDecay=1e-4,Regularization=L2,"
                   "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
  TString trainingStrategyString ("TrainingStrategy=");
  trainingStrategyString += training0 + "|" + training1 + "|" + training2;

  // General Options.
  TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                     "WeightInitialization=XAVIERUNIFORM");
  dnnOptions.Append (":"); dnnOptions.Append (layoutString);
  dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);

  TString cpuOptions = dnnOptions + ":Architecture=CPU";

  factory->BookMethod(dataloader, TMVA::Types::kDL, "DNN_CPU", cpuOptions);
  */

  factory->TrainAllMethods();
  factory->TestAllMethods();
  factory->EvaluateAllMethods();

  outputFile->Close();
  delete factory;
  delete dataloader;

  cout << '\n' << ">>>>> Not bad, not bad ... What's next? <<<<<" << '\n';

  if (!gROOT->IsBatch()) TMVA::TMVAGui("TMVA.root");

  return 0;
}
