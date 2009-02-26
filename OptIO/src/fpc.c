// $Id:$

const long long mask[8] =
{0x0000000000000000LL,
 0x00000000000000ffLL,
 0x000000000000ffffLL,
 0x0000000000ffffffLL,
 0x000000ffffffffffLL,
 0x0000ffffffffffffLL,
 0x00ffffffffffffffLL,
 0xffffffffffffffffLL};

void Compress(long predsizem1, 
              unsigned long long *inbuf, unsigned int inlen,
              unsigned char *outbuf, unsigned int *outlen)
{
  register long i, out, intot, hash, dhash, code, bcode, ioc;
  register long long val, lastval, stride, pred1, pred2, xor1, xor2;
  register long long *fcm, *dfcm;

  predsizem1 = (1L << predsizem1) - 1;

  hash = 0;
  dhash = 0;
  lastval = 0;
  pred1 = 0;
  pred2 = 0;
  fcm = (long long *)calloc(predsizem1 + 1, 8);
  dfcm = (long long *)calloc(predsizem1 + 1, 8);
  if (!fcm || !dfcm) {
    return;
  }

  intot = inlen;
  if (intot>0) {
    val = inbuf[0];
    out = 0 + ((intot + 1) >> 1);
    *((long long *)&outbuf[(out >> 3) << 3]) = 0;
    for (i = 0; i < intot; i += 2) {
      xor1 = val ^ pred1;
      fcm[hash] = val;
      hash = ((hash << 6) ^ ((unsigned long long)val >> 48)) & predsizem1;
      pred1 = fcm[hash];

      stride = val - lastval;
      xor2 = val ^ (lastval + pred2);
      lastval = val;
      val = inbuf[i + 1];
      dfcm[dhash] = stride;
      dhash = ((dhash << 2) ^ ((unsigned long long)stride >> 40)) & predsizem1;
      pred2 = dfcm[dhash];

      code = 0;
      if ((unsigned long long)xor1 > (unsigned long long)xor2) {
        code = 0x80;
        xor1 = xor2;
      }
      bcode = 7;                // 8 bytes
      if (0 == (xor1 >> 56))
        bcode = 6;              // 7 bytes
      if (0 == (xor1 >> 48))
        bcode = 5;              // 6 bytes
      if (0 == (xor1 >> 40))
        bcode = 4;              // 5 bytes
      if (0 == (xor1 >> 24))
        bcode = 3;              // 3 bytes
      if (0 == (xor1 >> 16))
        bcode = 2;              // 2 bytes
      if (0 == (xor1 >> 8))
        bcode = 1;              // 1 byte
      if (0 == xor1)
        bcode = 0;              // 0 bytes

      *((long long *)&outbuf[(out >> 3) << 3]) |= xor1 << ((out & 0x7) << 3);
      if (0 == (out & 0x7))
        xor1 = 0;
      *((long long *)&outbuf[((out >> 3) << 3) + 8]) = 
        (unsigned long long)xor1 >> (64 - ((out & 0x7) << 3));

      out += bcode + (bcode >> 2);
      code |= bcode << 4;

      xor1 = val ^ pred1;
      fcm[hash] = val;
      hash = ((hash << 6) ^ ((unsigned long long)val >> 48)) & predsizem1;
      pred1 = fcm[hash];

      stride = val - lastval;
      xor2 = val ^ (lastval + pred2);
      lastval = val;
      val = inbuf[i + 2];
      dfcm[dhash] = stride;
      dhash = ((dhash << 2) ^ ((unsigned long long)stride >> 40)) & predsizem1;
      pred2 = dfcm[dhash];

      bcode = code | 0x8;
      if ((unsigned long long)xor1 > (unsigned long long)xor2) {
        code = bcode;
        xor1 = xor2;
      }
      bcode = 7;                // 8 bytes
      if (0 == (xor1 >> 56))
        bcode = 6;              // 7 bytes
      if (0 == (xor1 >> 48))
        bcode = 5;              // 6 bytes
      if (0 == (xor1 >> 40))
        bcode = 4;              // 5 bytes
      if (0 == (xor1 >> 24))
        bcode = 3;              // 3 bytes
      if (0 == (xor1 >> 16))
        bcode = 2;              // 2 bytes
      if (0 == (xor1 >> 8))
        bcode = 1;              // 1 byte
      if (0 == xor1)
        bcode = 0;              // 0 bytes

      *((long long *)&outbuf[(out >> 3) << 3]) |= xor1 << ((out & 0x7) << 3);
      if (0 == (out & 0x7))
        xor1 = 0;
      *((long long *)&outbuf[((out >> 3) << 3) + 8]) = 
        (unsigned long long)xor1 >> (64 - ((out & 0x7) << 3));

      out += bcode + (bcode >> 2);
      outbuf[0 + (i >> 1)] = code | bcode;
    }
    if (0 != (intot & 1)) {
      out -= bcode + (bcode >> 2);
    }
  }
  *outlen=out;
  free(fcm);
  free(dfcm);
}

void Decompress(long predsizem1,
                unsigned char *inbuf, unsigned int inlen, 
                unsigned long long *outbuf, unsigned int *outlen)
{
  register long i, in, intot, hash, dhash, code, bcode, end, tmp, ioc;
  register long long val, lastval, stride, pred1, pred2, next;
  register long long *fcm, *dfcm;

  predsizem1 = (1L << predsizem1) - 1;

  hash = 0;
  dhash = 0;
  lastval = 0;
  pred1 = 0;
  pred2 = 0;
  fcm = (long long *)calloc(predsizem1 + 1, 8);
  dfcm = (long long *)calloc(predsizem1 + 1, 8);
  if (!fcm || !dfcm) {
    return;
  }

  intot = *outlen;
  in = inlen;

  if (intot>0) {
    in = (intot + 1) >> 1;
    for (i = 0; i < intot; i += 2) {
      code = inbuf[i >> 1];

      val = *((long long *)&inbuf[(in >> 3) << 3]);
      next = *((long long *)&inbuf[((in >> 3) << 3) + 8]);
      tmp = (in & 0x7) << 3;
      val = (unsigned long long)val >> tmp;
      next <<= 64 - tmp;
      if (0 == tmp)
        next = 0;
      val |= next;

      bcode = (code >> 4) & 0x7;
      val &= mask[bcode];
      in += bcode + (bcode >> 2);

      if (0 != (code & 0x80))
        pred1 = pred2;
      val ^= pred1;

      fcm[hash] = val;
      hash = ((hash << 6) ^ ((unsigned long long)val >> 48)) & predsizem1;
      pred1 = fcm[hash];

      stride = val - lastval;
      dfcm[dhash] = stride;
      dhash = ((dhash << 2) ^ ((unsigned long long)stride >> 40)) & predsizem1;
      pred2 = val + dfcm[dhash];
      lastval = val;

      outbuf[i] = val;

      val = *((long long *)&inbuf[(in >> 3) << 3]);
      next = *((long long *)&inbuf[((in >> 3) << 3) + 8]);
      tmp = (in & 0x7) << 3;
      val = (unsigned long long)val >> tmp;
      next <<= 64 - tmp;
      if (0 == tmp)
        next = 0;
      val |= next;

      bcode = code & 0x7;
      val &= mask[bcode];
      in += bcode + (bcode >> 2);

      if (0 != (code & 0x8))
        pred1 = pred2;
      val ^= pred1;

      fcm[hash] = val;
      hash = ((hash << 6) ^ ((unsigned long long)val >> 48)) & predsizem1;
      pred1 = fcm[hash];

      stride = val - lastval;
      dfcm[dhash] = stride;
      dhash = ((dhash << 2) ^ ((unsigned long long)stride >> 40)) & predsizem1;
      pred2 = val + dfcm[dhash];
      lastval = val;

      outbuf[i + 1] = val;
    }
  }
  free(fcm);
  free(dfcm);
}

void FLTCompressBS(long predsizem1, 
                   char *inbuf, unsigned int inlen,
                   unsigned char *outbuf, unsigned int *outlen)
{
  // test double
  int memsize1 = inlen / sizeof(double) + inlen % sizeof(double);
  char *memblock1 = (char*)calloc(memsize1*sizeof(double),1);
  memcpy(memblock1, inbuf, inlen);
  unsigned int outlen1 = *outlen*2+1024;
  char *outbuf1 = (char*)calloc(outlen1,1);
  Compress(predsizem1, (unsigned long long*)memblock1, memsize1, outbuf1, &outlen1);

  if (0) 
    printf("DecompressionD %d from %i to %i, yields %.2f%%\n", 
           predsizem1, inlen, outlen1, 100.*outlen1/inlen);

  if (outlen1<inlen) {
    *outlen=outlen1;
    memcpy(outbuf,memblock1,outlen1);
  } else {
    *outlen=inlen;
    memcpy(outbuf,inbuf,inlen);
  }

  free(memblock1);
  free(outbuf1);

#if 0
  // test float
  int memsize2 = inlen / sizeof(float) + inlen % sizeof(float);
  double *memblock2 = calloc(memsize2,sizeof(double));
  int i;
  char *iptr=(char*)memblock1;
  for(i=0;i<memsize2;++i) {
    float *fptr=(float*)iptr;
    double val=*fptr;
    memblock2[i]=val;
    iptr+=4;
  }
  unsigned int outlen2 = *outlen*2+1024;
  char *outbuf2 = (char*)calloc(outlen2,1);
  Compress(predsizem1, (unsigned long long*)memblock2, memsize2, outbuf2, &outlen2);

  printf("DecompressionF %d from %i to %i, yields %.2f%%\n", 
         predsizem1, inlen, outlen2, 100.*outlen2/inlen);

  free(memblock2);
  free(outbuf2);
#endif
}
