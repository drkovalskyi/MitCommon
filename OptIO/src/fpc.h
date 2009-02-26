// $Id:$

//
//Adapted from http://users.ices.utexas.edu/~burtscher/index.html
//

#ifndef FPC_DEFS_H
#define FPC_DEFS_H
void FLTCompress(long predsizem1, 
                 unsigned long long *inbuf, unsigned int inlen,
                 unsigned char *outbuf, unsigned int *outlen);

void FLTDecompress(long predsizem1,
                   unsigned char *inbuf, unsigned int inlen, 
                   unsigned long long *outbuf, unsigned int *outlen);

void FLTCompressBS(long predsizem1, 
                   char *inbuf, unsigned int inlen,
                   unsigned char *outbuf, unsigned int *outlen);
#endif
