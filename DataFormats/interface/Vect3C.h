//--------------------------------------------------------------------------------------------------
// $Id: Vect3C.h,v 1.1 2009/03/08 12:00:54 loizides Exp $
//
// Vect3C
//
// Implementation of our own ThreeVectorC32.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_DATAFORMATS_VECT3C_H
#define MITCOMMON_DATAFORMATS_VECT3C_H

#include "MitCommon/DataFormats/interface/Types.h"

namespace mithep 
{
  class Vect3C
  {
    public:
      Vect3C() : 
        fRho(0), fEta(0), fPhi(0) {}
      Vect3C(Double_t rho, Double_t eta, Double_t phi) : 
        fRho(rho), fEta(eta), fPhi(phi) {}
      Vect3C(const ThreeVector pos) :
        fRho(pos.Rho()), fEta(pos.Eta()), fPhi(pos.Phi()) {}
      Vect3C(const ThreeVectorC pos) :
        fRho(pos.Rho()), fEta(pos.Eta()), fPhi(pos.Phi()) {}

      Double_t            Eta()        const { return fEta; }
      Double_t            Phi()        const { return fPhi; }
      Double_t            Rho()        const { return fRho; }
      const ThreeVectorC  V()          const { return ThreeVectorC(fRho,fEta,fPhi); }
      void                SetXYZ(Double_t x, Double_t y, Double_t z);

    protected:
      Double32_t          fRho; //[0,0,12]rho-component
      Double32_t          fEta; //[0,0,10]eta-component
      Double32_t          fPhi; //[0,0,10]phi-component

    ClassDef(Vect3C, 1)
  };
}

//--------------------------------------------------------------------------------------------------
inline void mithep::Vect3C::SetXYZ(Double_t x, Double_t y, Double_t z)
{ 
  // Set four vector.

  ThreeVector tmp(x, y, z);
  fRho=tmp.Rho();
  fEta=tmp.Eta();
  fPhi=tmp.Phi();
}
#endif
