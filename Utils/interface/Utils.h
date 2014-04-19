//--------------------------------------------------------------------------------------------------
// Utils
//
// Utility functions useful when workign with root. Basically helping some of the interfaces to
// be more convenient of produce better debug output.
//
// Authors: C.Paus
//--------------------------------------------------------------------------------------------------
#ifndef MITCOMMON_UTILS_UTILS_H
#define MITCOMMON_UTILS_UTILS_H

#include <TString.h>
#include <TSystem.h>

namespace mithep
{
  class Utils
  {

    public:
      static TString     GetEnv(const char* name);    // get environment variable with proper check
      static TString     GetEnv(TString name);        // "
      static TString     DomainName();                // get domain name (uses HOSTNAME)
      static TString     GetCatalogDir(const char* name); // get the catalog directory
      static TString     GetJsonFile(const char* name);   // get the json file

    ClassDef(Utils, 0) // Utitily functions
  };
}

using namespace mithep;
//--------------------------------------------------------------------------------------------------
inline TString Utils::GetEnv(TString name)
{
  return GetEnv(name.Data());
}

inline TString Utils::GetEnv(const char* name)
{
  if (! gSystem->Getenv(name)) {
    printf(" Environment variable: %s  not defined. EXIT!\n",name);
    return TString("");
  } 
  return TString(gSystem->Getenv(name));  
}
inline TString Utils::DomainName()
{
  if (! gSystem->Getenv("HOSTNAME")) {
    printf(" Environment variable: HOSTNAME not defined. EXIT!\n");
    return TString("");
  } 
  TString domainName = gSystem->Getenv("HOSTNAME");
  domainName.Replace(0,domainName.First('.')+1,0,0);
  return domainName;
}

//--------------------------------------------------------------------------------------------------
inline TString Utils::GetCatalogDir(const char* dir)
{
  TString cataDir = TString("./catalog");
  Long_t *id=0,*size=0,*flags=0,*mt=0;

  printf(" Try local catalog first: %s\n",cataDir.Data());
  if (gSystem->GetPathInfo(cataDir.Data(),id,size,flags,mt) != 0) {
    cataDir = TString(dir);
    if (gSystem->GetPathInfo(cataDir.Data(),id,size,flags,mt) != 0) {
      printf(" Requested local (./catalog) and specified catalog do not exist. EXIT!\n");
      return TString("");
    }
  }
  else {
    printf(" Local catalog exists: %s using this one.\n",cataDir.Data()); 
  }

  return cataDir;
}

//--------------------------------------------------------------------------------------------------
inline TString Utils::GetJsonFile(const char* dir)
{
  TString jsonDir  = TString("./json");
  TString json     = Utils::GetEnv("MIT_PROD_JSON");
  Long_t *id=0,*size=0,*flags=0,*mt=0;

  printf(" Try local json first: %s\n",jsonDir.Data());
  if (gSystem->GetPathInfo(jsonDir.Data(),id,size,flags,mt) != 0) {
    jsonDir = TString(dir);
    if (gSystem->GetPathInfo(jsonDir.Data(),id,size,flags,mt) != 0) {
      printf(" Requested local (./json) and specified json directory do not exist. EXIT!\n");
      return TString("");
    }
  }
  else {
    printf(" Local json directory exists: %s using this one.\n",jsonDir.Data()); 
  }

  // Construct the full file name
  TString jsonFile = jsonDir + TString("/") + json;
  if (gSystem->GetPathInfo(jsonFile.Data(),id,size,flags,mt) != 0) {
    printf(" Requested jsonfile (%s) does not exist. EXIT!\n",jsonFile.Data());
    return TString("");
  }
  else {
    printf(" Requested jsonfile (%s) exist. Moving on now!\n",jsonFile.Data());
  }

  return jsonFile;
}
#endif
