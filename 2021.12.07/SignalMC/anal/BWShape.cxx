 /***************************************************************************** 
  * Project: RooFit                                                           * 
  *                                                                           * 
  * Copyright (c) 2000-2005, Regents of the University of California          * 
  *                          and Stanford University. All rights reserved.    * 
  *                                                                           * 
  * Redistribution and use in source and binary forms,                        * 
  * with or without modification, are permitted according to the terms        * 
  * listed in LICENSE (http://roofit.sourceforge.net/license.txt)             * 
  *****************************************************************************/ 

 // -- CLASS DESCRIPTION [PDF] -- 
 // Your description goes here... 

 #include <iostream> 

 #include "BWShape.h" 
 #include "RooAbsReal.h" 
 #include "TMath.h" 


 ClassImp(BWShape)

 // DecayType is not anymore used...
 BWShape::BWShape(const char *name, const char *title, 
                        RooAbsReal& _mass,
                        RooAbsReal& _mean,
		        RooAbsReal& _width,
		        RooAbsReal& _mass1,
		        RooAbsReal& _mass2,
		        RooAbsReal& _mass3,
		        RooAbsReal& _l,
	    	        RooAbsReal& _decayType,
		        RooAbsReal& _useBPhaseSpace
) :
   RooAbsPdf(name,title), 
   mass("mass","mass",this,_mass),
   mean("mean","mean",this,_mean),
   width("width","width",this,_width),
   mass1("mass1","mass1",this,_mass1),
   mass2("mass2","mass2",this,_mass2),
   mass3("mass3","mass3",this,_mass3),
   l("l","l",this,_l),
   decayType("decayType","decayType",this,_decayType),
   useBPhaseSpace("useBPhaseSpace","useBPhaseSpace",this,_useBPhaseSpace)
 { 
 } 


 BWShape::BWShape(const BWShape& other, const char* name) :  
   RooAbsPdf(other,name), 
   mass("mass",this,other.mass),
   mean("mean",this,other.mean),
   width("width",this,other.width),
   mass1("mass1",this,other.mass1),
   mass2("mass2",this,other.mass2),
   mass3("mass3",this,other.mass3),
   l("l",this,other.l),
   decayType("decayType",this,other.decayType),
   useBPhaseSpace("useBPhaseSpace",this,other.useBPhaseSpace)
{ 
 } 



 Double_t BWShape::evaluate() const 
 {   

   Double_t blattWeissKopfBarrier = 4.28629;

   Double_t thr = mass1 + mass2;
   if (mass <= thr) return 0;

   // The impulsion
   Double_t qMass = sqrt((mass*mass-(mass1+mass2)*(mass1+mass2)) * (mass*mass-(mass1-mass2)*(mass1-mass2)) )/(2*mass);
   Double_t qMean = sqrt((mean*mean-(mass1+mass2)*(mass1+mass2)) * (mean*mean-(mass1-mass2)*(mass1-mass2)) )/(2*mean);
   Double_t qMassForB = 0;
   if ((5.2*5.2-(mass+mass3)*(mass+mass3)) > 0) 
     qMassForB = sqrt((5.2*5.2-(mass+mass3)*(mass+mass3)) * (5.2*5.2-(mass-mass3)*(mass-mass3)) )/(2*5.2);
   
   // The width depending on the energy
   Double_t term1 = pow(qMass/qMean, 2*l+1);
   Double_t term2 = 0;
   if (l == 0)
     term2 = 1;
   else if (l == 1)
     term2 = ((1+pow(qMean*blattWeissKopfBarrier, 2))/(1+pow(qMass*blattWeissKopfBarrier, 2)));
       
   Double_t dependentWidth = term1 * term2 * width * mean / mass;

   //   Double_t inAndOut = 2 * pow(qMassForB * blattWeissKopfBarrier, 2) / (pow(qMassForB * blattWeissKopfBarrier, 2) + 1) * 
   //                       2 * pow(qMass * blattWeissKopfBarrier, 2) / (pow(qMass * blattWeissKopfBarrier, 2) + 1);

   // And the Breit-Wigner
   Double_t BW = ((mass * dependentWidth))/((mass*mass-mean*mean)*(mass*mass-mean*mean)+dependentWidth*dependentWidth*mean*mean );

   if (useBPhaseSpace) {
     return BW*qMassForB;
   }

   return BW;// * inAndOut;
 } 





