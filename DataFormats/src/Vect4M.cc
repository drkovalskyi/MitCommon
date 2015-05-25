// $Id: Vect4M.cc,v 1.1 2009/03/03 17:03:54 loizides Exp $

#include "MitCommon/DataFormats/interface/Vect4M.h"


void mithep::Vect4M::Set(Double_t pt, Double_t eta, Double_t phi, Double_t m)
{ 
  // Set four vector.

  fPt  = pt;
  fEta = eta;
  fM   = m;
  fPhi = phi;
}
//--------------------------------------------------------------------------------------------------
void mithep::Vect4M::SetXYZT(Double_t px, Double_t py, Double_t pz, Double_t e)
{ 
  // Set four vector.

  FourVector tmp(px, py, pz, e);
  fPt  = tmp.Pt();
  fEta = tmp.Eta();
  fPhi = tmp.Phi();
  fM   = tmp.M();
}


using namespace mithep;
