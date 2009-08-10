// $Id: $

#include "MitCommon/DataFormats/interface/TH2DAsymErr.h"
#include <TList.h>
#include <iostream>


ClassImp(mithep::TH2DAsymErr)

using namespace mithep;

//--------------------------------------------------------------------------------------------------
TH2DAsymErr::TH2DAsymErr(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup
                                 ,Int_t nbinsy, Double_t ylow, Double_t yup) : 
  TH2D(name, title, nbinsx, xlow, xup , nbinsy, ylow, yup)
{
  // Constructor.
  fStatErrorLow = TArrayD(fNcells);
  fStatErrorHigh = TArrayD(fNcells);
  fSysErrorLow = TArrayD(fNcells);
  fSysErrorHigh = TArrayD(fNcells);
}

TH2DAsymErr::TH2DAsymErr(const TH2D &h2d) : TH2D(h2d)
{
  // Constructor.
  fStatErrorLow = TArrayD(fNcells);
  fStatErrorHigh = TArrayD(fNcells);
  fSysErrorLow = TArrayD(fNcells);
  fSysErrorHigh = TArrayD(fNcells);
}

//--------------------------------------------------------------------------------------------------
void TH2DAsymErr::SetBinContent(Int_t binx, Int_t biny, Double_t value,
                              Double_t statErrorLow, Double_t statErrorHigh,
                              Double_t sysErrorLow, Double_t sysErrorHigh)
{
  // Fill histograms.

  Int_t bin;
  fEntries++;
  bin  = biny*(fXaxis.GetNbins()+2) + binx;

  SetBinContent(bin,value);
  if (fStatErrorLow.fN) fStatErrorLow.fArray[bin] = statErrorLow;
  if (fStatErrorHigh.fN) fStatErrorHigh.fArray[bin] = statErrorHigh;
  if (fSysErrorLow.fN) fSysErrorLow.fArray[bin] = sysErrorLow;
  if (fSysErrorHigh.fN) fSysErrorHigh.fArray[bin] = sysErrorHigh;
}

//--------------------------------------------------------------------------------------------------
void TH2DAsymErr::SetBinError(Int_t binx, Int_t biny,
                              Double_t statErrorLow, Double_t statErrorHigh,
                              Double_t sysErrorLow, Double_t sysErrorHigh)
{
  // Fill histograms.

  Int_t bin;
  fEntries++;
  bin  = biny*(fXaxis.GetNbins()+2) + binx;

  if (fStatErrorLow.fN) fStatErrorLow.fArray[bin] = statErrorLow;
  if (fStatErrorHigh.fN) fStatErrorHigh.fArray[bin] = statErrorHigh;
  if (fSysErrorLow.fN) fSysErrorLow.fArray[bin] = sysErrorLow;
  if (fSysErrorHigh.fN) fSysErrorHigh.fArray[bin] = sysErrorHigh;
}

//--------------------------------------------------------------------------------------------------
Double_t TH2DAsymErr::GetStatErrorLow(Double_t x, Double_t y)
{
  // Fill histograms.
  Int_t binx, biny, bin;

  binx = fXaxis.FindBin(x);
  biny = fYaxis.FindBin(y);
  bin  = biny*(fXaxis.GetNbins()+2) + binx;
  return fStatErrorLow.fArray[bin];
}


//--------------------------------------------------------------------------------------------------
Double_t TH2DAsymErr::GetStatErrorHigh(Double_t x, Double_t y)
{
  // Fill histograms.
  Int_t binx, biny, bin;
  binx = fXaxis.FindBin(x);
  biny = fYaxis.FindBin(y);
  bin  = biny*(fXaxis.GetNbins()+2) + binx;
  return fStatErrorHigh.fArray[bin];
}


//--------------------------------------------------------------------------------------------------
Double_t TH2DAsymErr::GetSysErrorLow(Double_t x, Double_t y)
{
  // Fill histograms.
  Int_t binx, biny, bin;
  binx = fXaxis.FindBin(x);
  biny = fYaxis.FindBin(y);
  bin  = biny*(fXaxis.GetNbins()+2) + binx;
  return fSysErrorLow.fArray[bin];
}


//--------------------------------------------------------------------------------------------------
Double_t TH2DAsymErr::GetSysErrorHigh(Double_t x, Double_t y)
{
  // Fill histograms.
  Int_t binx, biny, bin;
  binx = fXaxis.FindBin(x);
  biny = fYaxis.FindBin(y);
  bin  = biny*(fXaxis.GetNbins()+2) + binx;
  return fSysErrorHigh.fArray[bin];
}
