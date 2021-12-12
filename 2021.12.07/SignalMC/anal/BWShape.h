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

#ifndef BWSHAPE
#define BWSHAPE

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooAbsReal.h"
 
class BWShape : public RooAbsPdf {
public:
  BWShape(const char *name, const char *title,
	      RooAbsReal& _mass,
	      RooAbsReal& _mean,
	      RooAbsReal& _width,
	      RooAbsReal& _mass1,
	      RooAbsReal& _mass2,
	      RooAbsReal& _mass3,
	      RooAbsReal& _l,
	      RooAbsReal& _decayType,
	      RooAbsReal& _useBPhaseSpace);
  BWShape(const BWShape& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new BWShape(*this,newname); }
  inline virtual ~BWShape() { }
  

protected:

  RooRealProxy mass ;
  RooRealProxy mean ;
  RooRealProxy width ;
  RooRealProxy mass1 ;
  RooRealProxy mass2 ;
  RooRealProxy mass3 ;
  RooRealProxy l ;
  RooRealProxy decayType ;
  RooRealProxy useBPhaseSpace ;

  Double_t evaluate() const ;

private:

  ClassDef(BWShape,0) // Your description goes here...
};
 
#endif

