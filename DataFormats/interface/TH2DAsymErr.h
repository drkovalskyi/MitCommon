//--------------------------------------------------------------------------------------------------
// $Id $
//
// TH2DAsymErr
//
// Histogram that stores separate asymmetric statistical and systematic errors
//
// Authors: S.Xie
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_DATAFORMATS_HIST2DASYMERROR_H
#define MITCOMMON_DATAFORMATS_HIST2DASYMERROR_H

#include <TH2D.h>
class TCollection;

namespace mithep 
{
  class TH2DAsymErr : public TH2D
  {
    public:
      enum ErrorType { 
        kStatErrorLow = 0, 
        kStatErrorHigh, 
        kSysErrorLow, 
        kSysErrorHigh
      };

      TH2DAsymErr() {}
      TH2DAsymErr(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup,
                      Int_t nbinsy, Double_t ylow, Double_t yup);
      TH2DAsymErr(const TH2D &h2d);
      ~TH2DAsymErr() {}

      void                SetBinContent(Int_t binx, Int_t biny, Double_t value, 
                                        Double_t statErrorLow, Double_t statErrorHigh,
                                        Double_t sysErrorLow, Double_t sysErrorHigh);
      void                SetBinError(Int_t binx, Int_t biny, Double_t statErrorLow,                                      
                                      Double_t statErrorHigh, Double_t sysErrorLow,
                                      Double_t sysErrorHigh);

      Double_t            GetStatErrorLow(Double_t x, Double_t y);
      Double_t            GetStatErrorHigh(Double_t x, Double_t y);
      Double_t            GetSysErrorLow(Double_t x, Double_t y);
      Double_t            GetSysErrorHigh(Double_t x, Double_t y);

      using TH2D::SetBinContent;
      using TH2D::GetBinContent;

    protected:

      TArrayD             fStatErrorLow;  //Array to store statistical low errors
      TArrayD             fStatErrorHigh; //Array to store statistical high errors
      TArrayD             fSysErrorLow;   //Array to store systematic low errors
      TArrayD             fSysErrorHigh;  //Array to store systematic high errors

    ClassDef(TH2DAsymErr, 1) // Histogram ratio class in 1D
  };
}


#endif
