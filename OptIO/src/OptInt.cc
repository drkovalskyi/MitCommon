// $Id: Types.cc,v 1.1 2008/09/27 05:44:11 loizides Exp $

#include "MitCommon/OptIO/interface/OptInt.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TError.h>
#include <TEnv.h>
#include <Riostream.h>

ClassImp(mithep::OptInt)

using namespace mithep;

extern "C" void R__SetZipMode(int zipmode);

extern double lzipfrac;
extern double gzipfrac;
extern double bzipfrac;
extern int myverbose;
extern int mystaticm;

//--------------------------------------------------------------------------------------------------
void OptInt::SetAlgoFractions(Double_t lzo, Double_t gz, Double_t bz)
{
  // Fraction of compression that must be reached to accept result of the "heavy"
  // compression algorithm. Negative values turn of usage of the specific compression 
  // algorithm in ZipMode=5.

  lzipfrac = lzo;
  gzipfrac = gz;
  bzipfrac = bz;
}

//--------------------------------------------------------------------------------------------------
void OptInt::SetStaticMalloc(Bool_t b)
{
  // Set flag to enable the usage of static rather than heap memory in the compression
  // algorithm. No real increase in speed has been observed, therefore off by default.

  mystaticm = b;
}

//--------------------------------------------------------------------------------------------------
void OptInt::SetVerbose(Int_t vl)
{
  // Set verbosity level in the R__zip and R__unzip functions. Use 1 to only enable for
  // writing, 2 only for reading and 10 for both.

  myverbose = vl;
}

//--------------------------------------------------------------------------------------------------
void OptInt::SetZipMode(Int_t zipmode)
{
  // Set ZipMode to be used in R__zip. Supported values range from 1 to 5, where
  //  1 == gzip  (http://www.gzip.org/, standard in ROOT)
  //  2 == bzip2 (http://www.bzip.org/)
  //  3 == lzo   (http://www.oberhumer.com/opensource/lzo/)
  //  4 == rle   (http://bcl.comli.eu/)
  //  5 == combination of 1 to 4.

  R__SetZipMode(zipmode);
}
