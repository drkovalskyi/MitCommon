//--------------------------------------------------------------------------------------------------
// $Id: Funcs.h,v 1.2 2009/07/17 10:37:19 loizides Exp $
//
// Funcs
//
// Function definitions.
//
// Authors: C.Loizides, G.Ceballos
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_MATHTOOLS_FUNCS_H
#define MITCOMMON_MATHTOOLS_FUNCS_H

#include "MitCommon/DataFormats/interface/Types.h"
#include <TMath.h>

class TF1;

namespace mithep
{
  class Funcs {
    public:
      static Double_t  BreitGaus(Double_t x, Double_t m, Double_t mwidth, Double_t msig, 
                                 Double_t fintf, Double_t xmin, Double_t xmax);
      static Double_t  BreitGaus(Double_t *x, Double_t *par);
      static Double_t  BreitWignerGamma(Double_t x, Double_t mean, Double_t gamma);
      static Double_t  BreitWignerGamma(Double_t *x, Double_t *par);
      static Double_t  BreitWignerZ(Double_t x, Double_t mean, Double_t gamma);
      static Double_t  BreitWignerZ(Double_t *x, Double_t *par);
      static Double_t  ExpRange(Double_t x, Double_t lambda, Double_t xmin, Double_t xmax);
      static Double_t  ExpRange(Double_t *x, Double_t *par);
      static Double_t  ZLineShapePlusBkg(Double_t *x, Double_t *par);
      static TF1      *CreateZLineShapePlusBkg(Double_t norm=1, Double_t xmin=60, Double_t xmax=110,
                                               const char *n="ZLineShapePlusBgk");

    ClassDef(Funcs, 0) // Function definitions
  };
}
#endif
