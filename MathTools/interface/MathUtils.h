//--------------------------------------------------------------------------------------------------
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
    virtual ~MathUtils() {}

    public:
      static Double_t AddInQuadrature(Double_t a, Double_t b);  
      static void     CalcRatio(Double_t n1, Double_t n2, 
                                Double_t &r, Double_t &rlow, Double_t &rup, Int_t type);
      static Double_t DeltaPhi(Double_t phi1, Double_t phi2);
      template<class V1, class V2> 
      static Double_t DeltaPhi(const V1 &v1, const V2 &v2);
      template<class V1, class V2> 
      static Double_t DeltaPhi(const V1 *v1, const V2 *v2);
      static Double_t DeltaR(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2);
      template<class V1, class V2> 
      static Double_t DeltaR(const V1 *v1, const V2 *v2);
      template<class V1, class V2> 
      static Double_t DeltaR(const V1 &v1, const V2 &v2);
      static Double_t DeltaR2(Double_t phi1, Double_t eta1, Double_t phi2, Double_t eta2);
      template<class V1, class V2> 
      static Double_t DeltaR2(const V1 &v1, const V2 &v2);
      template<class V1, class V2> 
      static Double_t DeltaR2(const V1 *v1, const V2 *v2);
      static Double_t Eta2Theta(Double_t eta);
      static Double_t Theta2Eta(Double_t theta);
  };
}

//--------------------------------------------------------------------------------------------------
template<class V1, class V2>
Double_t mithep::MathUtils::DeltaPhi(const V1 &v1, const V2 &v2)
{
  // DeltaPhi between two given objects.

  return mithep::MathUtils::DeltaPhi(v1.Phi(),v2.Phi());
}

//--------------------------------------------------------------------------------------------------
template<class V1, class V2>
Double_t mithep::MathUtils::DeltaPhi(const V1 *v1, const V2 *v2)
{
  // DeltaPhi between two given objects.

  return mithep::MathUtils::DeltaPhi(v1->Phi(),v2->Phi());
}

//--------------------------------------------------------------------------------------------------
template<class V1, class V2>
Double_t mithep::MathUtils::DeltaR(const V1 &v1, const V2 &v2)
{
  // DeltaR between two given objects.

  return mithep::MathUtils::DeltaR(v1.Phi(),v1.Eta(),v2.Phi(),v2.Eta());
}

//--------------------------------------------------------------------------------------------------
template<class V1, class V2>
Double_t mithep::MathUtils::DeltaR(const V1 *v1, const V2 *v2)
{
  // DeltaR between two given objects

  return mithep::MathUtils::DeltaR(v1->Phi(),v1->Eta(),v2->Phi(),v2->Eta());
}

//--------------------------------------------------------------------------------------------------
template<class V1, class V2>
Double_t mithep::MathUtils::DeltaR2(const V1 &v1, const V2 &v2)
{
  // DeltaR between two given objects.

  return mithep::MathUtils::DeltaR2(v1.Phi(),v1.Eta(),v2.Phi(),v2.Eta());
}

//--------------------------------------------------------------------------------------------------
template<class V1, class V2>
Double_t mithep::MathUtils::DeltaR2(const V1 *v1, const V2 *v2)
{
  // DeltaR between two given objects.

  return mithep::MathUtils::DeltaR2(v1->Phi(),v1->Eta(),v2->Phi(),v2->Eta());
}
#endif
