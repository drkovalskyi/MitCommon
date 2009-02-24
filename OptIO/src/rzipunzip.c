// $Id: Types.cc,v 1.1 2008/09/27 05:44:11 loizides Exp $

#include "MitCommon/OptIO/src/rzipunzip.h"
#include "MitCommon/OptIO/src/zlib.h"
#include "MitCommon/OptIO/src/bzlib.h"
#include "MitCommon/OptIO/src/lzo/lzo1x.h"
#include "MitCommon/OptIO/src/rle.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef HDRSIZE
#define HDRSIZE 9
#else
#error "HDRSIZE already defined."
#endif

#ifndef MYSIZE
#define MYSIZE 5*1024*1024
#else
#error "MYSIZE already defined."
#endif

typedef unsigned char uch;  /* code assumes unsigned bytes; these type- */

extern int R__ZipMode;
char mymem[MYSIZE];
char *myptr = mymem;

// the following can be changed using the OptInt interface
double lzipfrac = 1;
double gzipfrac = 1;
double bzipfrac = 1;
int myverbose = 0;
int mystaticm = 0;

//--------------------------------------------------------------------------------------------------
void *mymalloc(size_t size)
{
  if(!mystaticm || size>MYSIZE)
    return malloc(size);

  if(myptr-mymem>MYSIZE-size) {
    myptr=mymem;
  }

  void *v=(void*)myptr;
  myptr+=size;
  return v;
}
 
//--------------------------------------------------------------------------------------------------
void myfree(void *ptr)
{
  if (mystaticm && (char*)ptr>=mymem && (char*)ptr<mymem+MYSIZE)
    return;
  free(ptr);
}

//--------------------------------------------------------------------------------------------------
void R__zip(int cxlevel, int *srcsize, char *src, int *tgtsize, char *tgt, int *irep)
     /* int cxlevel;                      compression level */
     /* int  *srcsize, *tgtsize, *irep;   source and target sizes, replay */
     /* char *tgt, *src;                  source and target buffers */
{
  int err               = 0;
  int method            = 0; 
  unsigned int in_size  = 0;
  unsigned int out_size = 0;

  if (*tgtsize <= HDRSIZE) {
    R__error("target buffer too small");
    return;
  }
  if (*srcsize > 0xffffff) {
    R__error("source buffer too big");
    return;
  }

  *irep        = 0;
  in_size      = (unsigned)(*srcsize); /* decompressed size */
  char *tgtptr = tgt+HDRSIZE;          /* compress data     */

  if (R__ZipMode == 5) { /*determine best of all methods*/
    int cont = 1;
    int msize  = -1;
    int zmode  = -1;
    int tgtlen = 2*in_size+64*1024;
    char *mptr1 = mymalloc(tgtlen);
    char *mptr2 = mymalloc(tgtlen);
    char *mptr  = 0;
    char *mdum  = 0;

    if (cont) {
      zmode = 4;
      msize = RLE_Compress(src, mptr1, in_size);
      mptr = mptr1;
      mdum = mptr2;
    }

    if (cont && lzipfrac>0) {
      err = lzo_init(); 
      if ( err == LZO_E_OK) {
        lzo_uint outlen = tgtlen;
        char *lwmem = mymalloc(LZO1X_999_MEM_COMPRESS);
        err = lzo1x_999_compress(src,in_size,mdum,&outlen,lwmem);
        if (err == LZO_E_OK) {
          if (outlen<lzipfrac*msize) {
            msize = outlen;
            zmode = 3;
            char *tmp = mptr;
            mptr = mdum;
            mdum = tmp;
          }
        }
        myfree(lwmem);
      }
    }

    if (cont && gzipfrac>0) {
      z_stream stream;
      stream.next_in   = (Bytef*)src;
      stream.avail_in  = (uInt)(in_size);
      stream.next_out  = (Bytef*)(mdum);
      stream.avail_out = (uInt)(*tgtsize);
      stream.zalloc    = (alloc_func)0;
      stream.zfree     = (free_func)0;
      stream.opaque    = (voidpf)0;

      err = deflateInit(&stream, cxlevel);
      if (err == Z_OK) {
        err = deflate(&stream, Z_FINISH);
        if (err == Z_STREAM_END) {
          err = deflateEnd(&stream);
          int outlen = stream.total_out;
          if (outlen<gzipfrac*msize) {
            msize = outlen;
            zmode = 1;
            char *tmp = mptr;
            mptr = mdum;
            mdum = tmp;
          }
        }
      }
    }

    if (cont && bzipfrac>0) {
      bz_stream stream;
      stream.next_in   = src;
      stream.avail_in  = (uInt)(in_size);
      stream.next_out  = mdum;
      stream.avail_out = (uInt)(*tgtsize);
      stream.bzalloc = 0;
      stream.bzfree  = 0;
      stream.opaque  = 0;

      err = BZ2_bzCompressInit(&stream, 9, 0, 1);
      if (err == BZ_OK) {
        err = BZ2_bzCompress(&stream, BZ_FINISH);
        if (err == BZ_STREAM_END) {
          BZ2_bzCompressEnd(&stream);
        }
        int outlen = stream.total_out_lo32;
        if (outlen<bzipfrac*msize) {
          msize = outlen;
          zmode = 2;
          char *tmp = mptr;
          mptr = mdum;
          mdum = tmp;
        }
      }
    }

    if (msize>=in_size) {
      out_size = in_size;
      memcpy(tgtptr,src,out_size);
      tgt[0] = 'X'; /* signature xx */
      tgt[1] = 'X';
      method = 0;
    } else {
      switch (zmode) {
        case 1:
          tgt[0] = 'Z'; /* signature zlib */
          tgt[1] = 'L';
          method = Z_DEFLATED; //==8
          break;
        case 2:
          tgt[0] = 'B'; /* signature bzlib */
          tgt[1] = 'Z';
          method = 2;
          break;
        case 3:
          tgt[0] = 'L'; /* signature lzolib */
          tgt[1] = 'O';
          method = 19;
          break;
        case 4:
          tgt[0] = 'R'; /* signature rlelib */
          tgt[1] = 'E';
          method = 4;
          break;
      }
      out_size = msize;
      memcpy(tgtptr,mptr,out_size);
    }
    myfree(mptr1);
    myfree(mptr2);
  } else if (R__ZipMode == 4) { /*rle*/
    int tgtlen = 2*in_size+1;
    char *lmem1 = mymalloc(tgtlen);

    out_size = RLE_Compress(src, lmem1, in_size);
    if (out_size>=in_size) {
      out_size = in_size;
      memcpy(tgtptr,src,out_size);
      tgt[0] = 'X'; /* signature xx */
      tgt[1] = 'X';
      method = 0;
    } else {
      memcpy(tgtptr,lmem1,out_size);
      tgt[0] = 'R'; /* signature rlelib */
      tgt[1] = 'E';
      method = 4;
    }
    myfree(lmem1);
  } else if (R__ZipMode == 3) { /*lzo*/
    err = lzo_init(); 
    if ( err != LZO_E_OK) {
      printf("error %d - lzo_init()\n", err);
        return;
    }

    lzo_uint tgtlen = 2*in_size+64*1024;
    char *lmem1 = mymalloc(tgtlen);
    char *lmem2 = 0;
    if (cxlevel<=1) {
      lmem2 = mymalloc(LZO1X_1_11_MEM_COMPRESS);
      err = lzo1x_1_11_compress(src,in_size,lmem1,&tgtlen,lmem2);
      method = 11;
    } else if (cxlevel==2) {
      lmem2 = mymalloc(LZO1X_1_12_MEM_COMPRESS);
      err = lzo1x_1_12_compress(src,in_size,lmem1,&tgtlen,lmem2);
      method = 12;
    } else if (cxlevel<=5) {
      lmem2 = mymalloc(LZO1X_1_15_MEM_COMPRESS);
      err = lzo1x_1_15_compress(src,in_size,lmem1,&tgtlen,lmem2);
      method = 15;
    } else {
      lmem2 = mymalloc(LZO1X_999_MEM_COMPRESS);
      err = lzo1x_999_compress(src,in_size,lmem1,&tgtlen,lmem2);
      method = 19;
    } 
    if (err != LZO_E_OK) {
      myfree(lmem1);
      myfree(lmem2);
      return;
    }

    if (tgtlen>=lzipfrac*in_size) {
      out_size = in_size;
      memcpy(tgtptr,src,out_size);
      tgt[0] = 'X'; /* signature xx */
      tgt[1] = 'X';
      method = 0;
    } else  {
      out_size  = tgtlen; /* compressed size */
      memcpy(tgtptr,lmem1,out_size);
      tgt[0] = 'L'; /* signature lzolib */
      tgt[1] = 'O';
    }
    myfree(lmem1);
    myfree(lmem2);
  } else if (R__ZipMode == 2) { /*bzip*/
    bz_stream stream;
    stream.next_in   = src;
    stream.avail_in  = (uInt)(in_size);
    stream.next_out  = tgtptr;
    stream.avail_out = (uInt)(*tgtsize);
    stream.bzalloc = 0;
    stream.bzfree  = 0;
    stream.opaque  = 0;

    err = BZ2_bzCompressInit(&stream, 9, 0, 1+(9-cxlevel)*25);
    if (err != BZ_OK) {
       printf("error %d in BZ2_bzCompressInit (bzlib)\n",err);
       return;
    }
    err = BZ2_bzCompress(&stream, BZ_FINISH);
    if (err != BZ_STREAM_END) {
      BZ2_bzCompressEnd(&stream);
      //printf("error %d in BZ2_bzCompress (bzlib) is not = %d\n",err,BZ_STREAM_END);
      return;
    }
    err = BZ2_bzCompressEnd(&stream);

    out_size  = stream.total_out_lo32; /* compressed size */
    if (out_size>=bzipfrac*in_size) {
      out_size = in_size;
      memcpy(tgtptr,src,out_size);
      tgt[0] = 'X'; /* signature xx */
      tgt[1] = 'X';
      method = 0;
    } else {
      tgt[0] = 'B'; /* signature bzlib */
      tgt[1] = 'Z';
      method = 2;
    }
  } else if (R__ZipMode == 1) { /*zip*/
    z_stream stream;
    stream.next_in   = (Bytef*)src;
    stream.avail_in  = (uInt)(in_size);
    stream.next_out  = (Bytef*)(tgtptr);
    stream.avail_out = (uInt)(*tgtsize);
    stream.zalloc    = (alloc_func)0;
    stream.zfree     = (free_func)0;
    stream.opaque    = (voidpf)0;

    err = deflateInit(&stream, cxlevel);
    if (err != Z_OK) {
       printf("error %d in deflateInit (zlib)\n",err);
       return;
    }
    err = deflate(&stream, Z_FINISH);
    if (err != Z_STREAM_END) {
       deflateEnd(&stream);
       //printf("error %d in deflate (zlib) is not = %d\n",err,Z_STREAM_END);
       return;
    }
    err = deflateEnd(&stream);

    out_size  = stream.total_out; /* compressed size */
    if (out_size>=gzipfrac*in_size) {
      out_size = in_size;
      memcpy(tgtptr,src,out_size);
      tgt[0] = 'X'; /* signature xx */
      tgt[1] = 'X';
      method = 0;
    } else {
      tgt[0] = 'Z'; /* signature zlib */
      tgt[1] = 'L';
      method = Z_DEFLATED;
    }
  } else if (R__ZipMode == 0) {
    printf("error: Old zip method not supported in this patch.");
    return;
  }

  // fill rest of header
  tgt[2] = (char)method;
  tgt[3] = (char)(out_size & 0xff);
  tgt[4] = (char)((out_size >> 8) & 0xff);
  tgt[5] = (char)((out_size >> 16) & 0xff);
  tgt[6] = (char)(in_size & 0xff);         
  tgt[7] = (char)((in_size >> 8) & 0xff);
  tgt[8] = (char)((in_size >> 16) & 0xff);
  *irep = out_size + HDRSIZE;

  if (myverbose==1||myverbose>9) {
    printf("R__zip:: zm=%d m=%d cl=%d: %c%c compressed %lu bytes into %lu bytes -> %.3f%%\n", 
           R__ZipMode, method, cxlevel, tgt[0], tgt[1], 
           (unsigned long)in_size, (unsigned long)out_size, (double)out_size/in_size*100.);
  }
}

//--------------------------------------------------------------------------------------------------
void R__unzip(int *srcsize, uch *src, int *tgtsize, uch *tgt, int *irep)
    /* Input: scrsize - size of input buffer                               */
    /*        src     - input buffer                                       */
    /*        tgtsize - size of target buffer                              */
    /* Output: tgt - target buffer (decompressed)                          */
    /*         irep - size of decompressed data                            */
    /*                0 - if error                                         */
{
  unsigned int osize   = 0;
  unsigned int ibufcnt = 0, obufcnt = 0;
  uch  *ibufptr = 0, *obufptr = 0;
  *irep = 0L;

  /*   C H E C K   H E A D E R   */
  if (*srcsize < HDRSIZE) {
    fprintf(stderr,"R__unzip: too small source\n");
    return;
  }

  char method = src[2];
  if (method!=0 && method!=2 && method!=4 && method!=Z_DEFLATED && method!=11
      && method!=12 && method!=15 && method!=19) {
    fprintf(stderr,"R__unzip: error in header -> unknown method %d\n", method);
    return;
  }

  if ((method==0          &&  (src[0] != 'X' || src[1] != 'X'))    ||
      (method==2          &&  (src[0] != 'B' || src[1] != 'Z'))    ||
      (method==4          &&  (src[0] != 'R' || src[1] != 'E'))    ||
      (method==Z_DEFLATED && ((src[0] != 'C' || src[1] != 'S') && 
                              (src[0] != 'Z' || src[1] != 'L')))   ||
      ((method==11 || method==12 || method==15 || method==19)  && 
       (src[0] != 'L' || src[1] != 'O'))) {
    fprintf(stderr,"R__unzip: error in header -> m=%d with %c%c\n",
            method, src[0], src[1]);
    return;
  }

  ibufptr = src + HDRSIZE;
  ibufcnt = (unsigned int)src[3] | ((unsigned int)src[4] << 8) | ((unsigned int)src[5] << 16);
  osize   = (unsigned int)src[6] | ((unsigned int)src[7] << 8) | ((unsigned int)src[8] << 16);
  obufptr = tgt;
  obufcnt = *tgtsize;

  if (obufcnt < osize) {
    fprintf(stderr,"R__unzip: too small target\n");
    return;
  }

  if (ibufcnt + HDRSIZE != *srcsize) {
    fprintf(stderr,"R__unzip: discrepancy in source length\n");
    return;
  }

  if (myverbose==2 || myverbose>9) {
    printf("R__unzip:: zm=%d m=%d: %c%c uncompressed %lu bytes from %lu bytes (%.3f%%)\n", 
           R__ZipMode, method, src[0], src[1], osize, ibufcnt, (double)ibufcnt/osize*100.);
  }

  if (method==0) { // apparently this is not reached since underlying ROOT code catches this
    if (ibufcnt!=osize) {
      fprintf(stderr,"R__unzip: error in header -> input should be output %d!=%d\n", 
              ibufcnt, osize);
      return;
    }
    memcpy(obufptr,ibufptr,osize);
    *irep = obufcnt;
    return;
  }

  /*   D E C O M P R E S S   D A T A  */
  if (src[0] == 'R' && src[1] == 'E') { /* rle format */
    RLE_Uncompress(ibufptr, obufptr, ibufcnt);
    *irep = obufcnt;
    return;
  } else if (src[0] == 'L' && src[1] == 'O') { /* lolib format */

    int err = lzo_init(); 
    if ( err != LZO_E_OK) {
      fprintf(stderr,"R__unzip: error %d in lzo_init (lzolib)\n",err);
      return;
    }
    lzo_uint new_len = obufcnt;
    err = lzo1x_decompress(ibufptr,ibufcnt,obufptr,&new_len, NULL);
    if ((err != LZO_E_OK) || (new_len!=osize)) {
      fprintf(stderr,"R__unzip: error %d (%d,%d) in lzo1x_decompress (lzolib)\n",
              err, new_len, osize);
      return;
    }
    *irep = obufcnt;
    return;
  } else if (src[0] == 'B' && src[1] == 'Z') { /* bzlib format */
    bz_stream stream; /* decompression stream */
    stream.next_in   = ibufptr;
    stream.avail_in  = ibufcnt;
    stream.next_out  = obufptr;
    stream.avail_out = obufcnt;
    stream.bzalloc   = 0;//my_balloc;
    stream.bzfree    = 0;//my_bfree;
    stream.opaque    = 0;

    int err = BZ2_bzDecompressInit(&stream,0,0);
    if (err != BZ_OK) {
      fprintf(stderr,"R__unzip: error %d in BZ2_bzDecompressInit (bzlib)\n",err);
      return;
    }
    err = BZ2_bzDecompress(&stream);
    if (err != BZ_STREAM_END) {
      fprintf(stderr,"R__unzip: error %d inBZ2_bzDecompress (bzlib)\n",err);
      BZ2_bzDecompressEnd(&stream);
      return;
    }
    BZ2_bzDecompressEnd(&stream);
    *irep = stream.total_out_lo32;
    return;
  } else if (src[0] == 'Z' && src[1] == 'L') { /* zlib format */
    z_stream stream; /* decompression stream */
    stream.next_in   = (Bytef*)(&src[HDRSIZE]);
    stream.avail_in  = (uInt)(*srcsize);
    stream.next_out  = (Bytef*)tgt;
    stream.avail_out = (uInt)(*tgtsize);
    stream.zalloc = 0; // my_zalloc;
    stream.zfree  = 0; //my_zfree;
    stream.opaque    = (voidpf)0;

    int err = inflateInit(&stream);
    if (err != Z_OK) {
      fprintf(stderr,"R__unzip: error %d in inflateInit (zlib)\n",err);
      return;
    }
    err = inflate(&stream, Z_FINISH);
    if (err != Z_STREAM_END) {
      inflateEnd(&stream);
      fprintf(stderr,"R__unzip: error %d in inflate (zlib)\n",err);
      return;
    }
    inflateEnd(&stream);
    *irep = stream.total_out;
    return;
  } else if (src[0] == 'C' && src[1] == 'S') { /* old zlib format */
    if (R__Inflate(&ibufptr, &ibufcnt, &obufptr, &obufcnt)) {
      fprintf(stderr,"R__unzip: error during decompression\n");
      return;
    }

    /* if (obufptr - tgt != osize) {
       There are some rare cases when a few more bytes are required */
    if (obufptr - tgt > *tgtsize) {
      fprintf(stderr,"R__unzip: discrepancy (%ld) with initial size: %ld, tgtsize=%d\n",
              (long)(obufptr - tgt),osize,*tgtsize);
      *irep = obufptr - tgt;
      return;
    }
    *irep = osize;
  } else {
    fprintf(stderr,"R__unzip: Format not supported -> m=%d with %d%d", method, src[0], src[1]);
    return;
  }
}
