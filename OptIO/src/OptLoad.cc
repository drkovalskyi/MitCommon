// $Id: OptLoad.cc,v 1.3 2009/02/26 16:18:07 loizides Exp $

#include "MitCommon/OptIO/interface/OptInt.h"
#include <Riostream.h>
#include <TEnv.h>
#include <TError.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>


using namespace mithep;
extern int activated;

//--------------------------------------------------------------------------------------------------
namespace {
  class OptIntLoad {
    public:
      OptIntLoad() {
        TString mylib(gSystem->Getenv("LD_PRELOAD"));
        if (mylib.Contains("MitCommonOptIO")) {
          activated = 123;
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
      }
  };
  OptIntLoad dummy;
}
