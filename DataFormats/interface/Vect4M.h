//--------------------------------------------------------------------------------------------------
// $Id: Vect4M.h,v 1.5 2009/03/25 05:16:12 loizides Exp $
//
// Vect4M
//
// Implementation of our own FourVectorM32.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_DATAFORMATS_VECT4M_H
#define MITCOMMON_DATAFORMATS_VECT4M_H

#include "MitCommon/DataFormats/interface/Types.h"

namespace mithep 
{
  class Vect4M
  {
    public:
      Vect4M() : 
        fPt(0), fEta(0), fPhi(0), fM(0) {}
      Vect4M(Double_t pt, Double_t eta, Double_t phi, Double_t m) : 
        fPt(pt), fEta(eta), fPhi(phi), fM(m) {}
      Vect4M(const FourVector &mom) : 
        fPt(mom.Pt()), fEta(mom.Eta()), fPhi(mom.Phi()), fM(mom.M()) {}
      Vect4M(const FourVectorM &mom) : 
        fPt(mom.Pt()), fEta(mom.Eta()), fPhi(mom.Phi()), fM(mom.M()) {}

      Double_t            Eta()        const { return fEta; }
      Double_t            Phi()        const { return fPhi; }
      Double_t            Pt()         const { return fPt; }
      Double_t            M()          const { return fM; }
      void                Set(Double_t pt, Double_t eta, Double_t phi, Double_t m);
      void                SetXYZT(Double_t px, Double_t py, Double_t pz, Double_t e);
      FourVectorM         V()          const { return FourVectorM(fPt,fEta,fPhi,fM); }

    protected:
      Double32_t          fPt;  //[0,0,14]pt-component
      Double32_t          fEta; //[0,0,12]eta-component
      Double32_t          fPhi; //[0,0,12]phi-component
      Double32_t          fM;   //[0,0,14]mass-component

    ClassDef(Vect4M, 1) // Implementation of our own FourVectorM32
  };
}

//--------------------------------------------------------------------------------------------------
inline void mithep::Vect4M::Set(Double_t pt, Double_t eta, Double_t phi, Double_t m)
{ 
  // Set four vector.

  fPt  = pt;
  fEta = eta;
  fPhi = phi;
  fM   = m;
}

//--------------------------------------------------------------------------------------------------
inline void mithep::Vect4M::SetXYZT(Double_t px, Double_t py, Double_t pz, Double_t e)
{ 
  // Set four vector.

  FourVector tmp(px, py, pz, e);
  fPt  = tmp.Pt();
  fEta = tmp.Eta();
  fPhi = tmp.Phi();
  fM   = tmp.M();
}
#endif
