// $Id: OptLoad.cc,v 1.2 2009/02/24 20:13:57 loizides Exp $

#include "MitCommon/OptIO/interface/OptInt.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TError.h>
#include <TEnv.h>
#include <Riostream.h>

using namespace mithep;

//--------------------------------------------------------------------------------------------------
namespace {
  class OptIntLoad {
    public:
      OptIntLoad() {
        Int_t verbose = gEnv->GetValue("Root.OptIO.Verbose",0);
        if(verbose)
          ::Info("OptIntLoad", "Loading libMitCommonOptIO.so for the optimized IO interface.");
        Int_t zipmode = gEnv->GetValue("Root.ZipMode",-1);
        if (zipmode>0) {
          if(verbose)
            ::Warning("OptIntLoad", "Setting ZipMode to %d as specified in rootrc", zipmode);
          OptInt::SetZipMode(zipmode);
        }
        Double_t lzf = gEnv->GetValue("Root.OptIO.LzoFraction",1.);
        Double_t gzf = gEnv->GetValue("Root.OptIO.GzipFraction",1.);
        Double_t bzf = gEnv->GetValue("Root.OptIO.BzipFraction",1.);
        Double_t lzm = gEnv->GetValue("Root.OptIO.LzmaFraction",1.);
        if (lzf!=1 || gzf!= 1 || bzf!=1 || lzm!=1) {
          if(verbose)
            ::Warning("OptIntLoad", "Setting algo fractions to %.2f, %.2f, %.2f, %.2f"
                      " as specified in rootrc", lzf, gzf, bzf, lzm);
          OptInt::SetAlgoFractions(lzf, gzf, bzf, lzm);
        }
        Int_t smalloc = gEnv->GetValue("Root.OptIO.SMalloc",0);
        if (smalloc) {
          if(verbose)
            ::Warning("OptIntLoad", "Enabling static memory allocation as specified in rootrc");
          OptInt::SetStaticMalloc(smalloc==1);
        }
        if (verbose) {
          ::Warning("OptIntLoad", "Setting verbosity to %d as specified in rootrc", verbose);
          OptInt::SetVerbose(verbose);
        }
      }
  };
  OptIntLoad dummy;
}
