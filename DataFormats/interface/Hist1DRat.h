//--------------------------------------------------------------------------------------------------
// $Id: Vect3.h,v 1.3 2009/03/20 12:54:20 loizides Exp $
//
// Hist1DRat
//
// Class to store numerator and denominator histograms in a mergeable fashion. 
// Not yet functional with out-of-the-box hadd command, and therefore not yet useful.
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_DATAFORMATS_HIST1DRAT_H
#define MITCOMMON_DATAFORMATS_HIST1DRAT_H

#include <TH1D.h>
class TBrowser;
class TCollection;

namespace mithep 
{
  class Hist1DRat : public TObject
  {
    public:
      Hist1DRat() : fRat(0) {}
      Hist1DRat(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup);
      ~Hist1DRat() { delete fRat; }

      void                Browse(TBrowser *b);
      void                Draw(Option_t *option = "");
      void                Fill(Double_t x, Double_t weight=1);
      void                FillDen(Double_t x, Double_t weight=1);
      void                FillNum(Double_t x, Double_t weight=1);
      const TH1D         *GetDen()   const { return &fDen; }
      const char         *GetName()  const;
      const TH1D         *GetNum()   const { return &fNum; }
      const TH1D         *GetRatio() const { return fRat;  }
      Long64_t            Merge(TCollection *list);
      void                Sumw2()          { fNum.Sumw2(); fDen.Sumw2();          }

    protected:
      void                Create();
      void                Reset()          { if (fRat) { delete fRat; fRat = 0; } }

      TH1D                fNum;      //numerator histogram
      TH1D                fDen;      //denumerator histogram
      TH1D               *fRat;      //!ratio histogram

    ClassDef(Hist1DRat, 1) // Histogram ratio class in 1D
  };
}

//--------------------------------------------------------------------------------------------------
inline const char *mithep::Hist1DRat::GetName() const
{
  // Get name.

  TString nname(fNum.GetName());
  TString rname(nname(0,nname.Length()-4));
  return Form("%s",rname.Data());
}
#endif
