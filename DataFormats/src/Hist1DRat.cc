// $Id: Vect3.cc,v 1.1 2009/03/08 12:00:54 loizides Exp $

#include "MitCommon/DataFormats/interface/Hist1DRat.h"
#include <TList.h>

ClassImp(mithep::Hist1DRat)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
Hist1DRat::Hist1DRat(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup)
  : fNum(Form("%s_num",name),title,nbinsx,xlow,xup),
    fDen(Form("%s_den",name),title,nbinsx,xlow,xup),
    fRat(0)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
void Hist1DRat::Browse(TBrowser *b)
{
  // Browse ratio histogram.

  Create();
  fRat->Browse(b);
}

//--------------------------------------------------------------------------------------------------
void Hist1DRat::Create()
{
  // Create ratio histogram if needed

  if (fRat)
    return;

  fRat = dynamic_cast<TH1D*>(fNum.Clone(GetName()));
  fRat->Divide(&fDen);
}

//--------------------------------------------------------------------------------------------------
void Hist1DRat::Draw(Option_t *option)
{
  // Create ratio histogram if needed and draw it afterwards.

  Create();
  fRat->Draw(option);
}

//--------------------------------------------------------------------------------------------------
void Hist1DRat::Fill(Double_t x, Double_t weight)
{
  // Fill histograms.

  fNum.Fill(x,weight);
  fDen.Fill(x,weight);
  Reset();
}

//--------------------------------------------------------------------------------------------------
void Hist1DRat::FillDen(Double_t x, Double_t weight)
{
  // Fill denominator.

  fDen.Fill(x,weight);
  Reset();
}


//--------------------------------------------------------------------------------------------------
void Hist1DRat::FillNum(Double_t x, Double_t weight)
{
  // Fill numerator.

  fNum.Fill(x,weight);
  Reset();
}

Long64_t Hist1DRat::Merge(TCollection *list)
{
  // Merge numerator and denominator histograms.

  TList nhs;
  TList dhs;

  TIter next(list);
  while (TObject *o = next()) {
    Hist1DRat *rat = dynamic_cast<Hist1DRat*>(o);
    if (rat) {
      nhs.Add(const_cast<TH1D*>(rat->GetNum()));
      dhs.Add(const_cast<TH1D*>(rat->GetDen()));
    }
  }

  Long64_t nres = fNum.Merge(&nhs);
  Long64_t dres = fDen.Merge(&dhs);

  Reset();
  return nres+dres;
}
