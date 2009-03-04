// $Id: OptInt.cc,v 1.3 2009/02/26 16:18:07 loizides Exp $

#include "MitCommon/OptIO/interface/OptInt.h"
#include <TROOT.h>
#include <TSystem.h>
#include <TError.h>
#include <TEnv.h>
#include <Riostream.h>

ClassImp(mithep::OptInt)

using namespace mithep;

extern "C" void R__SetZipMode(int zipmode);

extern int activated;
extern double lzipfrac;
extern double gzipfrac;
extern double bzipfrac;
extern double lzmafrac;
extern int myverbose;
extern int mystaticm;

//--------------------------------------------------------------------------------------------------
Bool_t OptInt::IsActivated()
{
  // Return true if code was properly pre-loaded.

  return (activated==123);
}

//--------------------------------------------------------------------------------------------------
void OptInt::SetAlgoFractions(Double_t lzo, Double_t gz, Double_t bz, Double_t lzma)
{
  // Fraction of compression that must be reached to accept result of the "heavy"
  // compression algorithms. Negative values turn off usage of the specific compression 
  // algorithm in ZipMode=99.

  lzipfrac = lzo;
  gzipfrac = gz;
  bzipfrac = bz;
  lzmafrac = lzma;
}

//--------------------------------------------------------------------------------------------------
void OptInt::SetBzipFraction(Double_t f)
{
  // Set fraction for compression algorithm, see description of SetAlgoFractions.

  bzipfrac = f;
}

//--------------------------------------------------------------------------------------------------
void OptInt::SetGzipFraction(Double_t f)
{
  // Set fraction for compression algorithm, see description of SetAlgoFractions.

  gzipfrac = f;
}

//--------------------------------------------------------------------------------------------------
void OptInt::SetLzmaFraction(Double_t f)
{
  // Set fraction for compression algorithm, see description of SetAlgoFractions.

  lzmafrac = f;
}

//--------------------------------------------------------------------------------------------------
void OptInt::SetLzoFraction(Double_t f)
{
  // Set fraction for compression algorithm, see description of SetAlgoFractions.

  lzipfrac = f;
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
  //  5 == lzma  (http://www.7-zip.org/sdk.html)
  //  99 == combination of 1 to 5.

  R__SetZipMode(zipmode);
}
