/*****************************************************************************
 * Project: RooFit                                                           *
 *                                                                           *
  * This code was autogenerated by RooClassFactory                            *
 *****************************************************************************/

#ifndef POL4EXPPDF
#define POL4EXPPDF

#include "RooAbsPdf.h"
#include "RooRealProxy.h"
#include "RooCategoryProxy.h"
#include "RooAbsReal.h"
#include "RooAbsCategory.h"

TString  nameParameters[] = {"p0","p1","p2","p3","p4","p5","N_bkg","N_sig","mean","width"};

class Pol4ExpPdf : public RooAbsPdf {
public:
  Pol4ExpPdf() {} ;
  Pol4ExpPdf(const char *name, const char *title,
	      RooAbsReal& _x,
	      RooAbsReal& _P0,
	      RooAbsReal& _P1,
        RooAbsReal& _P2,
        RooAbsReal& _P3,
        RooAbsReal& _P4,
        RooAbsReal& _P5);
  Pol4ExpPdf(const Pol4ExpPdf& other, const char* name=0) ;
  virtual TObject* clone(const char* newname) const { return new Pol4ExpPdf(*this,newname); }
  inline virtual ~Pol4ExpPdf() { }

protected:

  RooRealProxy x ;
  RooRealProxy P0 ;
  RooRealProxy P1 ;
  RooRealProxy P2 ;
  RooRealProxy P3 ;
  RooRealProxy P4 ;
  RooRealProxy P5 ;

  Double_t evaluate() const ;

private:

  ClassDef(Pol4ExpPdf,1) // Your description goes here...
};

#endif
