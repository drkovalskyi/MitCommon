//--------------------------------------------------------------------------------------------------
// $Id: MathUtils.h,v 1.2 2008/09/04 17:11:31 loizides Exp $
//
// MathUtils
//
// Math utility functions.
//
// Authors: S.Xie, C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_MATHTOOLS_MATHUTILS_H
#define MITCOMMON_MATHTOOLS_MATHUTILS_H

#include "MitCommon/DataFormats/interface/Types.h"
#include <TMath.h>

namespace mithep
{
  class MathUtils {
    public:
      static double AddInQuadrature(double a, double b);  
      static double DeltaPhi(double phi1, double phi2);
      static double DeltaR(const FourVector &v1, const FourVector &v2);
      static double DeltaR(double phi1, double eta1, double phi2, double eta2);
      static double Eta2Theta(double eta);
      static double Theta2Eta(double theta);
  };
}
#endif
