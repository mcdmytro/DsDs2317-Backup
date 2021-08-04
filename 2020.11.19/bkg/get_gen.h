// Class genT defined to get the genHep information like momentum, energy, pT 
//of the generated particle and facilitate in understanding of the background
// Author  : Vishal Place: Nara W. University  
// Version : 0.08
// Limitation, upto 10 daughters are handleded using  this class if PHOTOS is used
// for PHOTOS uses upto additional 3 gammas. Otherwise it can check upto 13 
// daughters.

#ifndef GET_GEN_H
#define GET_GEN_H

#include "belle.h"
#include "particle/Particle.h"

#if defined(BELLE_NAMESPACE)
namespace Belle{
#endif
  
  
  class genT
  { 
  public:
    genT(){};
    ~genT(void){};
    int c2_body(int cda1, int cda2, int da1, int da2);
    int c3_body(int cda1, int cda2, int cda3, int da1, int da2, int da3);
    int c4_body(int cda1, int cda2, int cda3, int cda4, int da1, int da2, int da3, int da4);
    int c5_body(int cda1, int cda2, int cda3, int cda4, int cda5, int da1, int da2, int da3, int da4, int da5);
    int c6_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int da1, int da2, int da3, int da4, int da5, int da6);
    int c7_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7, int da1, int da2, int da3, int da4, int da5, int da6, int da7);
    int c8_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7, int cda8, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8);
    int c9_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7, int cda8, int cda9, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9);
    
    int c10_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7, int cda8,  int cda9, int cda10, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10);
    int c11_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7, int cda8,  int cda9, int cda10, int cda11, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10, int da11);
    int c12_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7, int cda8,  int cda9, int cda10, int cda11, int cda12, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10, int da11, int da12);
    int c13_body(int cda1, int cda2, int cda3, int cda4, int cda5, int cda6, int cda7, int cda8,  int cda9, int cda10, int cda11, int cda12, int cda13, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10, int da11, int da12, int da13);

    
    int dau(Gen_hepevt a);
    
        
    // Check that da1, da2, ... are daughters of particle a
    // and return 1 on success
    // has been increase till 13 to take into PHOTOS into account for 10
    int checkDecay(Gen_hepevt a, int da1, int da2);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10, int da11);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10, int da11, int da12);
    int checkDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10, int da11, int da12, int da13);
    
    
    // with PHOTOS, here P stands for PHOTOS
    // and return 1 on success
    int Pcheck(Gen_hepevt a, int num);
    int PcheckDecay(Gen_hepevt a, int da1, int da2);
    int PcheckDecay(Gen_hepevt a, int da1, int da2, int da3);
    int PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4);
    int PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5);
    int PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6);
    int PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7);
    int PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8);
    int PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9);
    int PcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10);
    
    // For Absolute, Insensitive to charge
    // Check that da1, da2, ... are daughters of particle a
    // and return 1 on success
    // has been increase till 13 to take into PHOTOS into account for 10
    int AbscheckDecay(Gen_hepevt a, int da1, int da2);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10, int da11);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10, int da11, int da12);
    int AbscheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10, int da11, int da12, int da13);
    
    // For Absolute, Insensitive to charge
    // with PHOTOS, here P stands for PHOTOS
    // and return 2 on success
    int AbsPcheckDecay(Gen_hepevt a, int da1, int da2);
    int AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3);
    int AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4);
    int AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5);
    int AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6);
    int AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7);
    int AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8);
    int AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9);
    int AbsPcheckDecay(Gen_hepevt a, int da1, int da2, int da3, int da4, int da5, int da6, int da7, int da8, int da9, int da10);
    
    //.............................................................................................................................
    
    //Get id of particle a
    int AmId(Gen_hepevt a);
    //Print all daughters of a
    void PrintDAU(Gen_hepevt a);
    //Get id of daughter of a with rank useful in cases with two daughters having same id
    // for example B+ --> chic1 K- pi+ pi- pi+ pi+ You have 3 pi+ so their rank is 1,2,3
    //return the daughter of a after matching with id and corresponding rank
    Gen_hepevt DaME(Gen_hepevt a, int id);
    Gen_hepevt DaME(Gen_hepevt a, int id, int rank);
    Gen_hepevt AbsDaME(Gen_hepevt a, int id);
    Gen_hepevt AbsDaME(Gen_hepevt a, int id, int rank);
 
   
    //return the fourmomentum of a
    HepLorentzVector Mom4DAU(Gen_hepevt a, int id);
    HepLorentzVector Mom4DAU(Gen_hepevt a, int id, int rank);
    HepLorentzVector AbsMom4DAU(Gen_hepevt a, int id);
    HepLorentzVector AbsMom4DAU(Gen_hepevt a, int id, int rank);
    HepLorentzVector Mom4(Gen_hepevt a);

    //return momentum of a's daughter
    double MomDAU(Gen_hepevt a, int id);
    double MomDAU(Gen_hepevt a, int id, int rank);
    double AbsMomDAU(Gen_hepevt a, int id);
    double AbsMomDAU(Gen_hepevt a, int id, int rank);
    double Mom(Gen_hepevt a);
    double Mom(HepLorentzVector a);
    

    //return pT of a's daughter
    double pTDAU(Gen_hepevt a, int id);
    double pTDAU(Gen_hepevt a, int id, int rank);
    double AbspTDAU(Gen_hepevt a, int id);
    double AbspTDAU(Gen_hepevt a, int id, int rank);
    double pT(Gen_hepevt a);			
    double pT(HepLorentzVector a);
   

    //return child of a according to which number child
    Gen_hepevt MyChild(Gen_hepevt a, int rank);

  };
  
#if defined(BELLE_NAMESPACE)
} // namespace Belle
#endif

#endif // GET_GEN_H
