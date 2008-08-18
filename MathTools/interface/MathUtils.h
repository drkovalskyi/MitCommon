//--------------------------------------------------------------------------------------------------
// $Id: $
//
// MathUtils
//
// Math Utility functions
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MATHTOOLS_MATHUTILS_H
#define MATHTOOLS_MATHUTILS_H

#include "MitAna/DataTree/interface/Types.h"
#include <TMath.h>

namespace mitMath
{

  /// Convert double into a string
  std::string ftoa(double x);

  //Math functions

  double deltaPhi(double phi1, double phi2);
  double phiEtaDeltaR(double phi1, double eta1, double phi2, double eta2);
  double deltaR(mithep::FourVector v1, mithep::FourVector v2);
  double addInQuadrature(double a, double b);  
  double Eta2Theta(double eta);
  double Theta2Eta(double theta);
}
#endif
