//--------------------------------------------------------------------------------------------------
// $Id: MathUtils.h,v 1.6 2009/03/03 18:12:25 bendavid Exp $
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
      static Double_t AddInQuadrature(Double_t a, Double_t b);  
      static Double_t DeltaPhi(Double_t phi1, Double_t phi2);
      template<class V1, class V2> 
      static Double_t DeltaPhi(const V1 &v1, const V2 &v2);
      static Double_t DeltaR(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2);
      template<class V1, class V2> 
      static Double_t DeltaR(const V1 &v1, const V2 &v2);
      static Double_t Eta2Theta(Double_t eta);
      static Double_t Theta2Eta(Double_t theta);
  };
}

//--------------------------------------------------------------------------------------------------
template<class V1, class V2>
Double_t mithep::MathUtils::DeltaPhi(const V1 &v1, const V2 &v2)
{
  // DeltaPhi between two given objects

  return mithep::MathUtils::DeltaPhi(v1.Phi(),v2.Phi());
}

//--------------------------------------------------------------------------------------------------
template<class V1, class V2>
Double_t mithep::MathUtils::DeltaR(const V1 &v1, const V2 &v2)
{
  // DeltaR between two given objects

  return mithep::MathUtils::DeltaR(v1.Phi(),v1.Eta(),v2.Phi(),v2.Eta());
}
#endif
