//--------------------------------------------------------------------------------------------------
// $Id: MathUtils.h,v 1.5 2009/02/18 15:38:26 loizides Exp $
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
      static Double_t DeltaPhi(const FourVector &v1, const FourVector &v2);
      static Double_t DeltaPhi(const FourVectorM &v1, const FourVectorM &v2);
      static Double_t DeltaR(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2);
      static Double_t DeltaR(const FourVector &v1, const FourVector &v2);
      static Double_t DeltaR(const FourVectorM &v1, const FourVectorM &v2);
      template<class V1, class V2> static Double_t DeltaR(const V1 &v1, const V2 &v2);
      static Double_t Eta2Theta(Double_t eta);
      static Double_t Theta2Eta(Double_t theta);
  };
}

//--------------------------------------------------------------------------------------------------
template<class V1, class V2>
Double_t mithep::MathUtils::DeltaR(const V1 &v1, const V2 &v2)
{

  return mithep::MathUtils::DeltaR(v1.Phi(),v1.Eta(),v2.Phi(),v2.Eta());

}

#endif
