// $Id: Types.cc,v 1.1 2008/09/27 05:44:11 loizides Exp $

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
        ::Info("OptIntLoad", "Loading libMitCommonOptIO.so for the optimized IO interface.");
        Int_t zipmode = gEnv->GetValue("Root.ZipMode",-1);
        if (zipmode>0) {
          ::Warning("OptIntLoad", "Setting ZipMode to %d as specified in rootrc", zipmode);
          OptInt::SetZipMode(zipmode);
        }
        Double_t lzf = gEnv->GetValue("Root.OptIO.LzoFraction",1.);
        Double_t gzf = gEnv->GetValue("Root.OptIO.GzipFraction",1.);
        Double_t bzf = gEnv->GetValue("Root.OptIO.BzipFraction",1.);
        if (lzf!=1 || gzf!= 1 || bzf!=1 ) {
          ::Warning("OptIntLoad", "Setting algo fractions to %f %f %f as specified in rootrc", 
                    lzf, gzf, bzf);
          OptInt::SetAlgoFractions(lzf, gzf, bzf);
        }
        Int_t smalloc = gEnv->GetValue("Root.OptIO.SMalloc",0);
        if (smalloc) {
          ::Warning("OptIntLoad", "Enabling static memory allocation as specified in rootrc");
          OptInt::SetStaticMalloc(smalloc==1);
        }
        Int_t verbose = gEnv->GetValue("Root.OptIO.Verbose",0);
        if (verbose) {
          ::Warning("OptIntLoad", "Setting verbosity to %d as specified in rootrc", verbose);
          OptInt::SetVerbose(verbose);
        }
      }
  };
  OptIntLoad dummy;
}
