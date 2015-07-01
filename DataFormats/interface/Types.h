//--------------------------------------------------------------------------------------------------
// Types
//
// Here we define common types. For help on the GenVector lib, see
// http://root.cern.ch/root/html/MATH_GENVECTOR_Index.html
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_DATAFORMATS_TYPES_H
#define MITCOMMON_DATAFORMATS_TYPES_H
 
#include <vector>
#include <Rtypes.h>
#include <Math/GenVector/LorentzVector.h>
#include <Math/SMatrix.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/Vector3Dfwd.h>
#include <Math/Vector4Dfwd.h>

namespace mithep
{
  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<Double_t> >   FourVector;
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double_t> > FourVectorM;
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<Double_t> > FourVectorE;

  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double_t>,
                                           ROOT::Math::DefaultCoordinateSystemTag> 
                                                                        ThreeVector;
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::CylindricalEta3D<Double_t>,
                                           ROOT::Math::DefaultCoordinateSystemTag> 
                                                                        ThreeVectorC;

  typedef ROOT::Math::SMatrix<Double_t,3,3,ROOT::Math::MatRepSym<Double_t,3> >   ThreeSymMatrix;
  typedef ROOT::Math::SMatrix<Double_t,7,7,ROOT::Math::MatRepSym<Double_t,7> >   SevenSymMatrix;
  typedef ROOT::Math::SMatrix<Double_t,3,3,ROOT::Math::MatRepStd<Double_t,3,3> >    ThreeMatrix;
  typedef ROOT::Math::SMatrix<Double_t,7,7,ROOT::Math::MatRepStd<Double_t,7,7> >    SevenMatrix;

  typedef ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<Double32_t> >   FourVector32;
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<Double32_t> > FourVectorM32;
  typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiE4D<Double32_t> > FourVectorE32;

  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::Cartesian3D<Double32_t>,
                                           ROOT::Math::DefaultCoordinateSystemTag> 
                                                                          ThreeVector32;
  typedef ROOT::Math::DisplacementVector3D<ROOT::Math::CylindricalEta3D<Double32_t>,
                                           ROOT::Math::DefaultCoordinateSystemTag> 
                                                                          ThreeVectorC32;

  typedef ROOT::Math::SMatrix<Double32_t,3,3,ROOT::Math::MatRepSym<Double32_t,3> > ThreeSymMatrix32;
  typedef ROOT::Math::SMatrix<Double32_t,7,7,ROOT::Math::MatRepSym<Double32_t,7> > SevenSymMatrix32;
  typedef ROOT::Math::SMatrix<Double32_t,3,3,ROOT::Math::MatRepStd<Double32_t,3,3> >  ThreeMatrix32;
  typedef ROOT::Math::SMatrix<Double32_t,7,7,ROOT::Math::MatRepStd<Double32_t,7,7> >  SevenMatrix32;
  
  typedef std::pair<UInt_t,UInt_t>  UIntPair;
  typedef std::vector<UIntPair> UIntBounds;
}
#endif
