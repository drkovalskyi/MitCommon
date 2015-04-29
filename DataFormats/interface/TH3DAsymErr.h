//--------------------------------------------------------------------------------------------------
// TH3DAsymErr
//
// Histogram that stores separate asymmetric statistical and systematic errors. It is 
// derived from TH3D.
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_DATAFORMATS_HIST3DASYMERROR_H
#define MITCOMMON_DATAFORMATS_HIST3DASYMERROR_H

#include <TH3D.h>
#include "TH2DAsymErr.h"

namespace mithep 
{
  class TH3DAsymErr : public TH3D
  {
    public:

      TH3DAsymErr() {}
      TH3DAsymErr(const char *name, const char *title, 
                  Int_t nbinsx, Double_t xlow, Double_t xup,
                  Int_t nbinsy, Double_t ylow, Double_t yup,
                  Int_t nbinsz, Double_t zlow, Double_t zup);
      TH3DAsymErr(const TH3D &h3d);
      ~TH3DAsymErr() {}

      using TH3D::GetBinContent;
      Double_t  GetError(Double_t x, Double_t y, Double_t z, TH2DAsymErr::EErrType t);
      using TH3D::GetBinError;
      Double_t  GetBinError(Int_t b, Int_t c, Int_t d, TH2DAsymErr::EErrType t);
      Double_t  GetStatErrorLow(Double_t x, Double_t y, Double_t z);
      Double_t  GetBinStatErrorLow(Int_t b, Int_t c, Int_t d);
      Double_t  GetStatErrorHigh(Double_t x, Double_t y, Double_t z);
      Double_t  GetBinStatErrorHigh(Int_t b, Int_t c, Int_t d);
      Double_t  GetSysErrorLow(Double_t x, Double_t y, Double_t z);
      Double_t  GetBinSysErrorLow(Int_t b, Int_t c, Int_t d);      
      Double_t  GetSysErrorHigh(Double_t x, Double_t y, Double_t z);
      Double_t  GetBinSysErrorHigh(Int_t b, Int_t c, Int_t d);
      using TH3D::SetBinContent;
      void      SetBinContent(Int_t binx, Int_t biny, Int_t binz, Double_t value, 
                              Double_t statErrorLow, Double_t statErrorHigh,
                              Double_t sysErrorLow, Double_t sysErrorHigh);
      using TH3D::SetBinError;
      void      SetBinError(Int_t binx, Int_t biny, Int_t binz,
                            Double_t statErrorLow, Double_t statErrorHigh, 
                            Double_t sysErrorLow, Double_t sysErrorHigh);
      
    protected:

      TArrayD       fStatErrorLow;  //array to store statistical low errors
      TArrayD       fStatErrorHigh; //array to store statistical high errors
      TArrayD       fSysErrorLow;   //array to store systematic low errors
      TArrayD       fSysErrorHigh;  //array to store systematic high errors

    ClassDef(TH3DAsymErr, 1) // Histogram for storage of seperate asymmetric errors
  };
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::TH3DAsymErr::GetStatErrorLow(Double_t x, Double_t y, Double_t z)
{
  // Get Upper Statistical Uncertainty
  return GetError(x,y,z,TH2DAsymErr::kStatErrLow);
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::TH3DAsymErr::GetStatErrorHigh(Double_t x, Double_t y, Double_t z)
{
  // Get Upper Statistical Uncertainty
  return GetError(x,y,z,TH2DAsymErr::kStatErrHigh);
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::TH3DAsymErr::GetBinStatErrorLow(Int_t b, Int_t c, Int_t d)
{
  // Get Upper Statistical Uncertainty
  return GetBinError(b,c,d,TH2DAsymErr::kStatErrLow);
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::TH3DAsymErr::GetBinStatErrorHigh(Int_t b, Int_t c, Int_t d)
{
  // Get Upper Statistical Uncertainty
  return GetBinError(b,c,d,TH2DAsymErr::kStatErrHigh);
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::TH3DAsymErr::GetSysErrorLow(Double_t x, Double_t y, Double_t z)
{
  // Get Upper Sysistical Uncertainty
  return GetError(x,y,z,TH2DAsymErr::kSysErrLow);
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::TH3DAsymErr::GetSysErrorHigh(Double_t x, Double_t y, Double_t z)
{
  // Get Upper Sysistical Uncertainty
  return GetError(x,y,z,TH2DAsymErr::kSysErrHigh);
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::TH3DAsymErr::GetBinSysErrorLow(Int_t b, Int_t c, Int_t d)
{
  // Get Upper Sysistical Uncertainty
  return GetBinError(b,c,d,TH2DAsymErr::kSysErrLow);
}

//--------------------------------------------------------------------------------------------------
inline Double_t mithep::TH3DAsymErr::GetBinSysErrorHigh(Int_t b, Int_t c, Int_t d)
{
  // Get Upper Sysistical Uncertainty
  return GetBinError(b,c,d,TH2DAsymErr::kSysErrHigh);
}

#endif
