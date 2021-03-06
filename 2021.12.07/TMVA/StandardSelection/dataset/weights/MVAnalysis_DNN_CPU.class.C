// Class: ReadDNN_CPU
// Automatically generated by MethodBase::MakeClass
//

/* configuration options =====================================================

#GEN -*-*-*-*-*-*-*-*-*-*-*- general info -*-*-*-*-*-*-*-*-*-*-*-

Method         : DL::DNN_CPU
TMVA Release   : 4.2.1         [262657]
ROOT Release   : 6.22/02       [398850]
Creator        : mcdmytro
Date           : Thu Dec  2 18:07:24 2021
Host           : Linux MSI 5.10.16.3-microsoft-standard-WSL2 #1 SMP Fri Apr 2 22:23:49 UTC 2021 x86_64 x86_64 x86_64 GNU/Linux
Dir            : /mnt/c/shared/DsDs2317-GitSynch/DsDs2317-study/TMVA/PreSelectionApplied
Training events: 5000
Analysis type  : [Classification]


#OPT -*-*-*-*-*-*-*-*-*-*-*-*- options -*-*-*-*-*-*-*-*-*-*-*-*-

# Set by User:
V: "True" [Verbose output (short form of "VerbosityLevel" below - overrides the latter one)]
VarTransform: "N" [List of variable transformations performed before training, e.g., "D_Background,P_Signal,G,N_AllClasses" for: "Decorrelation, PCA-transformation, Gaussianisation, Normalisation, each for the given class of events ('AllClasses' denotes all events of all classes, if no class indication is given, 'All' is assumed)"]
H: "False" [Print method-specific help message]
Layout: "TANH|128,TANH|128,TANH|128,LINEAR" [Layout of the network.]
ErrorStrategy: "CROSSENTROPY" [Loss function: Mean squared error (regression) or cross entropy (binary classification).]
WeightInitialization: "XAVIERUNIFORM" [Weight initialization strategy]
Architecture: "CPU" [Which architecture to perform the training on.]
TrainingStrategy: "LearningRate=1e-2,Momentum=0.9,Repetitions=1,ConvergenceSteps=30,BatchSize=256,TestRepetitions=10,WeightDecay=1e-4,Regularization=None,DropConfig=0.0+0.5+0.5+0.5," [Defines the training strategies.]
# Default:
VerbosityLevel: "Verbose" [Verbosity level]
CreateMVAPdfs: "False" [Create PDFs for classifier outputs (signal and background)]
IgnoreNegWeightsInTraining: "False" [Events with negative weights are ignored in the training (but are included for testing and performance evaluation)]
InputLayout: "0|0|0" [The Layout of the input]
BatchLayout: "0|0|0" [The Layout of the batch]
RandomSeed: "0" [Random seed used for weight initialization and batch shuffling]
ValidationSize: "20%" [Part of the training data to use for validation. Specify as 0.2 or 20% to use a fifth of the data set as validation set. Specify as 100 to use exactly 100 events. (Default: 20%)]
##


#VAR -*-*-*-*-*-*-*-*-*-*-*-* variables *-*-*-*-*-*-*-*-*-*-*-*-

NVar 5
pst_ds17                      pst_ds17                      pst_ds17                      pst_ds17                                                        'F'    [0.0704832673073,4.73005151749]
hel_phi                       hel_phi                       hel_phi                       hel_phi                                                         'F'    [0.000175855340785,0.999992072582]
pc2_ds17                      pc2_ds17                      pc2_ds17                      pc2_ds17                                                        'F'    [0,0.999999940395]
pc2_dsii                      pc2_dsii                      pc2_dsii                      pc2_dsii                                                        'F'    [0,0.999359250069]
pc2_pi0                       pc2_pi0                       pc2_pi0                       pc2_pi0                                                         'F'    [5.58792194738e-08,1]
NSpec 2
dm_ds17                       dm_ds17                       dm_ds17                       F                                                               'F'    [0.250198841095,0.449937462807]
mc_ds17                       mc_ds17                       mc_ds17                       F                                                               'F'    [0,10431]


============================================================================ */

#include <array>
#include <vector>
#include <cmath>
#include <string>
#include <iostream>

#ifndef IClassifierReader__def
#define IClassifierReader__def

class IClassifierReader {

 public:

   // constructor
   IClassifierReader() : fStatusIsClean( true ) {}
   virtual ~IClassifierReader() {}

   // return classifier response
   virtual double GetMvaValue( const std::vector<double>& inputValues ) const = 0;

   // returns classifier status
   bool IsStatusClean() const { return fStatusIsClean; }

 protected:

   bool fStatusIsClean;
};

#endif

class ReadDNN_CPU : public IClassifierReader {

 public:

   // constructor
   ReadDNN_CPU( std::vector<std::string>& theInputVars )
      : IClassifierReader(),
        fClassName( "ReadDNN_CPU" ),
        fNvars( 5 )
   {
      // the training input variables
      const char* inputVars[] = { "pst_ds17", "hel_phi", "pc2_ds17", "pc2_dsii", "pc2_pi0" };

      // sanity checks
      if (theInputVars.size() <= 0) {
         std::cout << "Problem in class \"" << fClassName << "\": empty input vector" << std::endl;
         fStatusIsClean = false;
      }

      if (theInputVars.size() != fNvars) {
         std::cout << "Problem in class \"" << fClassName << "\": mismatch in number of input values: "
                   << theInputVars.size() << " != " << fNvars << std::endl;
         fStatusIsClean = false;
      }

      // validate input variables
      for (size_t ivar = 0; ivar < theInputVars.size(); ivar++) {
         if (theInputVars[ivar] != inputVars[ivar]) {
            std::cout << "Problem in class \"" << fClassName << "\": mismatch in input variable names" << std::endl
                      << " for variable [" << ivar << "]: " << theInputVars[ivar].c_str() << " != " << inputVars[ivar] << std::endl;
            fStatusIsClean = false;
         }
      }

      // initialize min and max vectors (for normalisation)
      fVmin[0] = -1;
      fVmax[0] = 1;
      fVmin[1] = -1;
      fVmax[1] = 1;
      fVmin[2] = -1;
      fVmax[2] = 1;
      fVmin[3] = -1;
      fVmax[3] = 0.99999988079071;
      fVmin[4] = -1;
      fVmax[4] = 1;

      // initialize input variable types
      fType[0] = 'F';
      fType[1] = 'F';
      fType[2] = 'F';
      fType[3] = 'F';
      fType[4] = 'F';

      // initialize constants
      Initialize();

      // initialize transformation
      InitTransform();
   }

   // destructor
   virtual ~ReadDNN_CPU() {
      Clear(); // method-specific
   }

   // the classifier response
   // "inputValues" is a vector of input values in the same order as the
   // variables given to the constructor
   double GetMvaValue( const std::vector<double>& inputValues ) const override;

 private:

   // method-specific destructor
   void Clear();

   // input variable transformation

   double fOff_1[3][5];
   double fScal_1[3][5];
   void InitTransform_1();
   void Transform_1( std::vector<double> & iv, int sigOrBgd ) const;
   void InitTransform();
   void Transform( std::vector<double> & iv, int sigOrBgd ) const;

   // common member variables
   const char* fClassName;

   const size_t fNvars;
   size_t GetNvar()           const { return fNvars; }
   char   GetType( int ivar ) const { return fType[ivar]; }

   // normalisation of input variables
   double fVmin[5];
   double fVmax[5];
   double NormVariable( double x, double xmin, double xmax ) const {
      // normalise to output range: [-1, 1]
      return 2*(x - xmin)/(xmax - xmin) - 1.0;
   }

   // type of input variable: 'F' or 'I'
   char   fType[5];

   // initialize internal variables
   void Initialize();
   double GetMvaValue__( const std::vector<double>& inputValues ) const;

   // private members (method specific)
inline double ReadDNN_CPU::GetMvaValue( const std::vector<double>& inputValues ) const
{
   // classifier response value
   double retval = 0;

   // classifier response, sanity check first
   if (!IsStatusClean()) {
      std::cout << "Problem in class \"" << fClassName << "\": cannot return classifier response"
                << " because status is dirty" << std::endl;
   }
   else {
         std::vector<double> iV(inputValues);
         Transform( iV, -1 );
         retval = GetMvaValue__( iV );
   }

   return retval;
}

//_______________________________________________________________________
inline void ReadDNN_CPU::InitTransform_1()
{
   double fMin_1[3][5];
   double fMax_1[3][5];
   // Normalization transformation, initialisation
   fMin_1[0][0] = 0.597169399261;
   fMax_1[0][0] = 4.73005151749;
   fScal_1[0][0] = 2.0/(fMax_1[0][0]-fMin_1[0][0]);
   fOff_1[0][0] = fMin_1[0][0]*fScal_1[0][0]+1.;
   fMin_1[1][0] = 0.0704832673073;
   fMax_1[1][0] = 4.62342357635;
   fScal_1[1][0] = 2.0/(fMax_1[1][0]-fMin_1[1][0]);
   fOff_1[1][0] = fMin_1[1][0]*fScal_1[1][0]+1.;
   fMin_1[2][0] = 0.0704832673073;
   fMax_1[2][0] = 4.73005151749;
   fScal_1[2][0] = 2.0/(fMax_1[2][0]-fMin_1[2][0]);
   fOff_1[2][0] = fMin_1[2][0]*fScal_1[2][0]+1.;
   fMin_1[0][1] = 0.012682701461;
   fMax_1[0][1] = 0.999992072582;
   fScal_1[0][1] = 2.0/(fMax_1[0][1]-fMin_1[0][1]);
   fOff_1[0][1] = fMin_1[0][1]*fScal_1[0][1]+1.;
   fMin_1[1][1] = 0.000175855340785;
   fMax_1[1][1] = 0.999483644962;
   fScal_1[1][1] = 2.0/(fMax_1[1][1]-fMin_1[1][1]);
   fOff_1[1][1] = fMin_1[1][1]*fScal_1[1][1]+1.;
   fMin_1[2][1] = 0.000175855340785;
   fMax_1[2][1] = 0.999992072582;
   fScal_1[2][1] = 2.0/(fMax_1[2][1]-fMin_1[2][1]);
   fOff_1[2][1] = fMin_1[2][1]*fScal_1[2][1]+1.;
   fMin_1[0][2] = 0;
   fMax_1[0][2] = 0.999999940395;
   fScal_1[0][2] = 2.0/(fMax_1[0][2]-fMin_1[0][2]);
   fOff_1[0][2] = fMin_1[0][2]*fScal_1[0][2]+1.;
   fMin_1[1][2] = 0;
   fMax_1[1][2] = 0.999998986721;
   fScal_1[1][2] = 2.0/(fMax_1[1][2]-fMin_1[1][2]);
   fOff_1[1][2] = fMin_1[1][2]*fScal_1[1][2]+1.;
   fMin_1[2][2] = 0;
   fMax_1[2][2] = 0.999999940395;
   fScal_1[2][2] = 2.0/(fMax_1[2][2]-fMin_1[2][2]);
   fOff_1[2][2] = fMin_1[2][2]*fScal_1[2][2]+1.;
   fMin_1[0][3] = 0;
   fMax_1[0][3] = 0.999359250069;
   fScal_1[0][3] = 2.0/(fMax_1[0][3]-fMin_1[0][3]);
   fOff_1[0][3] = fMin_1[0][3]*fScal_1[0][3]+1.;
   fMin_1[1][3] = 0;
   fMax_1[1][3] = 0.996803998947;
   fScal_1[1][3] = 2.0/(fMax_1[1][3]-fMin_1[1][3]);
   fOff_1[1][3] = fMin_1[1][3]*fScal_1[1][3]+1.;
   fMin_1[2][3] = 0;
   fMax_1[2][3] = 0.999359250069;
   fScal_1[2][3] = 2.0/(fMax_1[2][3]-fMin_1[2][3]);
   fOff_1[2][3] = fMin_1[2][3]*fScal_1[2][3]+1.;
   fMin_1[0][4] = 1.90725387483e-07;
   fMax_1[0][4] = 0.999999344349;
   fScal_1[0][4] = 2.0/(fMax_1[0][4]-fMin_1[0][4]);
   fOff_1[0][4] = fMin_1[0][4]*fScal_1[0][4]+1.;
   fMin_1[1][4] = 5.58792194738e-08;
   fMax_1[1][4] = 1;
   fScal_1[1][4] = 2.0/(fMax_1[1][4]-fMin_1[1][4]);
   fOff_1[1][4] = fMin_1[1][4]*fScal_1[1][4]+1.;
   fMin_1[2][4] = 5.58792194738e-08;
   fMax_1[2][4] = 1;
   fScal_1[2][4] = 2.0/(fMax_1[2][4]-fMin_1[2][4]);
   fOff_1[2][4] = fMin_1[2][4]*fScal_1[2][4]+1.;
}

//_______________________________________________________________________
inline void ReadDNN_CPU::Transform_1( std::vector<double>& iv, int cls) const
{
   // Normalization transformation
   if (cls < 0 || cls > 2) {
   if (2 > 1 ) cls = 2;
      else cls = 2;
   }
   const int nVar = 5;

   // get indices of used variables

   // define the indices of the variables which are transformed by this transformation
   static std::vector<int> indicesGet;
   static std::vector<int> indicesPut;

   if ( indicesGet.empty() ) {
      indicesGet.reserve(fNvars);
      indicesGet.push_back( 0);
      indicesGet.push_back( 1);
      indicesGet.push_back( 2);
      indicesGet.push_back( 3);
      indicesGet.push_back( 4);
   }
   if ( indicesPut.empty() ) {
      indicesPut.reserve(fNvars);
      indicesPut.push_back( 0);
      indicesPut.push_back( 1);
      indicesPut.push_back( 2);
      indicesPut.push_back( 3);
      indicesPut.push_back( 4);
   }

   static std::vector<double> dv;
   dv.resize(nVar);
   for (int ivar=0; ivar<nVar; ivar++) dv[ivar] = iv[indicesGet.at(ivar)];
   for (int ivar=0;ivar<5;ivar++) {
      double offset = fOff_1[cls][ivar];
      double scale  = fScal_1[cls][ivar];
      iv[indicesPut.at(ivar)] = scale*dv[ivar]-offset;
   }
}

//_______________________________________________________________________
inline void ReadDNN_CPU::InitTransform()
{
   InitTransform_1();
}

//_______________________________________________________________________
inline void ReadDNN_CPU::Transform( std::vector<double>& iv, int sigOrBgd ) const
{
   Transform_1( iv, sigOrBgd );
}
