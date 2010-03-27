// $Id: TH3DAsymErr.cc,v 1.3 2009/10/21 15:05:22 sixie Exp $

#include "MitCommon/DataFormats/interface/TH3DAsymErr.h"
#include <TList.h>
#include <iostream>

ClassImp(mithep::TH3DAsymErr)

using namespace mithep;
using namespace std;

//--------------------------------------------------------------------------------------------------
TH3DAsymErr::TH3DAsymErr(const char *name, const char *title, 
                         Int_t nbinsx, Double_t xlow, Double_t xup,
                         Int_t nbinsy, Double_t ylow, Double_t yup,
                         Int_t nbinsz, Double_t zlow, Double_t zup) : 
  TH3D(name, title, nbinsx, xlow, xup , nbinsy, ylow, yup, nbinsz, zlow, zup),
  fStatErrorLow(fNcells),
  fStatErrorHigh(fNcells),
  fSysErrorLow(fNcells),
  fSysErrorHigh(fNcells)
{
  // Constructor.
}

TH3DAsymErr::TH3DAsymErr(const TH3D &h3d) : 
  TH3D(h3d),
  fStatErrorLow(fNcells),
  fStatErrorHigh(fNcells),
  fSysErrorLow(fNcells),
  fSysErrorHigh(fNcells)
{
  // Constructor.
}

//--------------------------------------------------------------------------------------------------
Double_t TH3DAsymErr::GetError(Double_t x, Double_t y, Double_t z, TH2DAsymErr::EErrType t)
{
  // Get error corresponding to x,y,z for given error type.

  Int_t nx = fXaxis.GetNbins()+2;
  Int_t ny = fYaxis.GetNbins()+2;
  Int_t binx = fXaxis.FindFixBin(x);
  Int_t biny = fYaxis.FindFixBin(y);
  Int_t binz = fZaxis.FindFixBin(z);
  Int_t bin  =  binx + nx*(biny + ny*binz);
  
  switch (t) {
    case TH2DAsymErr::kStatErrLow:
      return fStatErrorLow.fArray[bin];
      break;
    case TH2DAsymErr::kStatErrHigh:
      return fStatErrorHigh.fArray[bin];
      break;
    case TH2DAsymErr::kSysErrLow:
      return fSysErrorLow.fArray[bin];
      break;
    case TH2DAsymErr::kSysErrHigh:
      return fSysErrorHigh.fArray[bin];
      break;
  }
  return 0;
}

//--------------------------------------------------------------------------------------------------
Double_t TH3DAsymErr::GetBinError(Int_t b, Int_t c, Int_t d, TH2DAsymErr::EErrType t)
{
  // Get error corresponding to bin (b,c) for given error type.

  Int_t nx = fXaxis.GetNbins()+2;
  Int_t ny = fYaxis.GetNbins()+2;
  Int_t bin = b + nx*(c + ny*d);

  switch (t) {
    case TH2DAsymErr::kStatErrLow:      
      return fStatErrorLow.fArray[bin];
      break;
    case TH2DAsymErr::kStatErrHigh:
      return fStatErrorHigh.fArray[bin];
      break;
    case TH2DAsymErr::kSysErrLow:
      return fSysErrorLow.fArray[bin];
      break;
    case TH2DAsymErr::kSysErrHigh:
      return fSysErrorHigh.fArray[bin];
      break;
  }
  return 0;
}


//--------------------------------------------------------------------------------------------------
void TH3DAsymErr::SetBinContent(Int_t binx, Int_t biny, Int_t binz, Double_t value,
                                Double_t statErrorLow, Double_t statErrorHigh,
                                Double_t sysErrorLow, Double_t sysErrorHigh)
{
  // Set bin content for given bin and value.

  Int_t bin = GetBin(binx,biny,binz);
  SetBinContent(bin,value);

  fStatErrorLow.fArray[bin]  = statErrorLow;
  fStatErrorHigh.fArray[bin] = statErrorHigh;
  fSysErrorLow.fArray[bin]   = sysErrorLow;
  fSysErrorHigh.fArray[bin]  = sysErrorHigh;
  return;
}

//--------------------------------------------------------------------------------------------------
void TH3DAsymErr::SetBinError(Int_t binx, Int_t biny, Int_t binz,
                              Double_t statErrorLow, Double_t statErrorHigh,
                              Double_t sysErrorLow, Double_t sysErrorHigh)
{
  // Set bin error for asymetric errors.

  Int_t bin = GetBin(binx,biny,binz);

  fStatErrorLow.fArray[bin] = statErrorLow;
  fStatErrorHigh.fArray[bin] = statErrorHigh;
  fSysErrorLow.fArray[bin] = sysErrorLow;
  fSysErrorHigh.fArray[bin] = sysErrorHigh;

  return;
}
