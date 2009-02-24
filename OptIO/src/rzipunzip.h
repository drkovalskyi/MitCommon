//--------------------------------------------------------------------------------------------------
// $Id: Types.h,v 1.4 2008/11/05 10:46:48 bendavid Exp $
//
// ZipInt
//
//
// Authors: C.Loizides
//--------------------------------------------------------------------------------------------------

#ifndef MITCOMMON_OPTIO_ZIPINT_H
#define MITCOMMON_OPTIO_ZIPINT_H

void R__zip(int cxlevel, int *srcsize, char *src, int *tgtsize, char *tgt, int *irep);
void R__unzip(int *srcsize, unsigned char *src, int *tgtsize, unsigned char *tgt, int *irep);
#endif
