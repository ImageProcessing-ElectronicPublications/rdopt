
/*
** Copyright 1995,1996 by Viresh Ratnakar, Miron Livny
**
** Permission to use and copy this software and its documentation
** for any non-commercial purpose and without fee is hereby granted,
** provided that the above copyright notice appear in all copies and that
** both that copyright notice and this permission notice appear in
** supporting documentation.
**
**
** The University of Wisconsin and the copyright holders make no
** representations about the suitability of this software for any
** purpose.  It is provided "as is" without express or implied warranty.
**
**
** THE UNIVERSITY OF WISCONSIN AND THE COPYRIGHT HOLDERS DISCLAIM ALL
** WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED WARRANTIES
** OF MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE UNIVERSITY OF
** WISCONSIN OR THE COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
** OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
** OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
** OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
** OR PERFORMANCE OF THIS SOFTWARE.
**
** Author:  Viresh Ratnakar
** Email:   ratnakar@cs.wisc.edu
**
*/

#include "rdopt.h"
#include "qentry.h"

#define H(n) job->Histogram[(n)]

extern void PrepareForErrBpp(OptimJob *job, int unum)
{
    int i,n, offset, totalentries, max;
    int *iptr;
    FFLOAT *fptr;

    /*** set MinTable and MaxTable ***/
    totalentries = 64;
    if (job->VerboseLevel > 1)
    {
        fprintf(errfile,"\t\tQuantTable Min/Max entries:\n");
        fflush(errfile);
    }
    for (n=0; n<64; n++)
    {
        max = H(n).PlusSize;
        if (max < H(n).MinusSize) max = H(n).MinusSize;
        if (max < job->MaxTable[unum][n]) job->MaxTable[unum][n] = max;
        if ((!n) && (job->DCclamp < job->MaxTable[unum][0]))
            job->MaxTable[unum][0] = job->DCclamp;
        if (job->MaxTable[unum][n] < job->MinTable[unum][n])
            job->MaxTable[unum][n] = job->MinTable[unum][n];

        if (job->VerboseLevel > 1)
        {
            if ((n%8)==0) fprintf(errfile,"\t\t");
#if (QTABBITS==8)
            fprintf(errfile,"%3d/%3d ",
#elif (QTABBITS==10)
            fprintf(errfile,"%4d/%4d ",
#elif (QTABBITS==12)
            fprintf(errfile,"%4d/%4d ",
#elif (QTABBITS==16)
            fprintf(errfile,"%5d/%5d ",
#else
            /** force syntax error **/
            QTABBITS must be 8,12, or 16
#endif
                    job->MinTable[unum][n],job->MaxTable[unum][n]);
            if ((n%8)==7)
            {
                fprintf(errfile,"\n");
                fflush(errfile);
            }
        }

        if (job->MapQ)
        {
            job->MaxTable[unum][n] = UnMapQentry(n, job->MaxTable[unum][n]);
            job->MinTable[unum][n] = UnMapQentry(n, job->MinTable[unum][n]);
        }
        totalentries += job->MaxTable[unum][n];
    }

    totalentries *= (job->ThreshSpan+1);

    if ((iptr = (int *) calloc(1,totalentries*sizeof(int)))==NULL)
    {
        FatalError("PrepareForErrBpp out of memory");
    }
    if ((fptr = (FFLOAT *) calloc(1,totalentries*sizeof(FFLOAT)))==NULL)
    {
        FatalError("PrepareForErrBpp out of memory");
    }

    offset = 0;
    for (n=0; n<64; n++)
    {
        job->opt_method.dp.Err[unum][n] = fptr + offset;
        job->opt_method.dp.Bpp[unum][n] = iptr + offset;
        offset += ((job->MaxTable[unum][n]+1)*(job->ThreshSpan+1));
    }
}



/*** MACROS for Err and Bpp ***/

#define SetErrPrepareThresh(E) \
{ \
  int i,intquant, ilim, idx ; \
  FFLOAT orig, err,  ans; \
      ans = ((FFLOAT) 0.0); \
      if (q < H(n).PlusSize) { \
        if ((ilim = (q + job->ThreshSpan)) >= H(n).PlusSize) \
	   ilim = H(n).PlusSize - 1; \
        for (i=H(n).PlusSize-1;i>ilim;i--) { \
          orig = UnDiscretizePlus(i); \
          intquant = QuantizeDis(i,q); \
          err = orig - RealUnQuantize(intquant,q); \
          ans += (err*err*((FFLOAT) H(n).PlusCount[i])); \
        } \
        idx = ilim - q; \
        for (i=ilim;i>=q;i--,idx--) { \
          orig = UnDiscretizePlus(i); \
          intquant = QuantizeDis(i,q); \
          err = orig - RealUnQuantize(intquant,q); \
	  ErrSavedPlus[idx] = (err*err*((FFLOAT) H(n).PlusCount[i])); \
          ans += ErrSavedPlus[idx]; \
        } \
        for (i=q-1;i>=0;i--) { \
          err = UnDiscretizePlus(i); \
          ans += (err*err*((FFLOAT) H(n).PlusCount[i])); \
        } \
      } \
      else { \
        for (i=H(n).PlusSize-1;i>=0;i--) { \
          err = UnDiscretizePlus(i); \
          ans += (err*err*((FFLOAT) H(n).PlusCount[i])); \
        } \
      } \
      if (q < H(n).MinusSize) { \
        if ((ilim = (q + job->ThreshSpan)) >= H(n).MinusSize) \
	   ilim = 1 - H(n).MinusSize ; \
        else ilim = 0 - ilim; \
        for (i=1-H(n).MinusSize;i<ilim;i++) { \
          orig = UnDiscretizeMinus(i); \
          intquant = QuantizeDis(i,q); \
          err = orig - RealUnQuantize(intquant,q); \
          ans += (err*err*((FFLOAT) H(n).MinusCount[0-i])); \
        } \
        idx = 0 - ilim - q; \
        for (i=ilim;i<=(0-q);i++, idx--) { \
          orig = UnDiscretizeMinus(i); \
          intquant = QuantizeDis(i,q); \
          err = orig - RealUnQuantize(intquant,q); \
	  ErrSavedMinus[idx] = (err*err*((FFLOAT) H(n).MinusCount[0-i])); \
          ans += ErrSavedMinus[idx]; \
        } \
        for (i=(1-q);i<=0;i++) { \
          err = UnDiscretizeMinus(i); \
          ans += (err*err*((FFLOAT) H(n).MinusCount[0-i])); \
        } \
      } \
      else { \
        for (i=1-H(n).MinusSize;i<=0;i++) { \
          err = UnDiscretizeMinus(i); \
          ans += (err*err*((FFLOAT) H(n).MinusCount[0-i])); \
        } \
      } \
      if (job->WeightedCoefs[unum]) { \
	ans *= job->CoefWeights[unum][n]; \
	for (i = 0; i<= job->ThreshSpan; i++) { \
	  ErrSavedPlus[i] *= job->CoefWeights[unum][n]; \
	  ErrSavedMinus[i] *= job->CoefWeights[unum][n]; \
	} \
      } \
      LastErr = ans; \
      E = ans*job->CompWeights[unum]; \
}

#define SetErr(E) \
{ \
  /* t2/2 is the threshold */  \
  FFLOAT err2, err; \
  int t2less1; \
      E = LastErr; \
      t2less1 = t2 - 1; \
      err2 = ((FFLOAT) 0.0); \
      if (t2less1 < H(n).PlusSize) { \
	E -= ErrSavedPlus[t2less1 - q]; \
	err = UnDiscretizePlus(t2less1); \
        err2 = (err*err*((FFLOAT) H(n).PlusCount[t2less1])); \
      } \
      if (t2less1 < H(n).MinusSize) { \
	E -= ErrSavedMinus[t2less1 - q]; \
	err = UnDiscretizeMinus((0-t2less1)); \
        err2 += (err*err*((FFLOAT) H(n).MinusCount[t2less1])); \
      } \
      if (job->WeightedCoefs[unum]) err2 *= job->CoefWeights[unum][n]; \
      E += err2; \
      LastErr = E; \
      E *= job->CompWeights[unum]; \
}

#define SetErrNoThresh(E) \
{ \
  int i,intquant; \
  FFLOAT orig, err,  ans; \
      ans = ((FFLOAT) 0.0); \
      for (i=H(n).PlusSize-1;i>=0;i--) { \
        orig = UnDiscretizePlus(i); \
        intquant = QuantizeDis(i,q); \
        err = orig - RealUnQuantize(intquant,q); \
        ans += (err*err*((FFLOAT) H(n).PlusCount[i])); \
      } \
      for (i=1-H(n).MinusSize;i<=0;i++) { \
        orig = UnDiscretizeMinus(i); \
        intquant = QuantizeDis(i,q); \
        err = orig - RealUnQuantize(intquant,q); \
        ans += (err*err*((FFLOAT) H(n).MinusCount[0-i])); \
      } \
      if (job->WeightedCoefs[unum]) ans *= job->CoefWeights[unum][n]; \
      E = ans*job->CompWeights[unum]; \
}

#define SetBppPrepareThresh(B) \
{\
  int v, lastVal,lastCount,quant;\
  int vh,vl,qvh,qvl,qv; \
  FFLOAT Entropy,entr_temp;\
  vh = q + job->ThreshSpan; \
  qvh = QuantizeDis(vh, q); \
  vl = 0 - q - job->ThreshSpan; \
  qvl = QuantizeDis(vl, q); \
  Entropy = 0.0;\
  lastVal = -5000; lastCount = 0;\
  for (v=1-H(n).MinusSize;v<=0;v++) {\
      quant = QuantizeDis(v, q);\
      if (quant > lastVal ) {\
        if (lastCount > 0) {\
	  entr_temp = \
                   ((((logtotal - ((FFLOAT) log10((double) lastCount)))\
                         * ((FFLOAT) lastCount) *\
                         job->bppWeight [unum]\
                         / job->NumBlocks[unum])\
                        )/19.26592);\
        } \
        else { \
	  entr_temp = ((FFLOAT) 0.0); \
        } \
        Entropy += entr_temp; \
        if (lastVal >= qvl) { \
	    SavedCounts[lastVal] = lastCount; \
	    SavedEntr[lastVal] = entr_temp; \
        } \
        lastCount = 0;\
        lastVal = quant;\
      }\
      lastCount += H(n).MinusCount[0-v];\
  }\
  for (v=0;v<H(n).PlusSize;v++) {\
      quant = QuantizeDis(v, q);\
      if (quant > lastVal ) {\
        if (lastCount > 0) {\
          entr_temp = \
                   ((((logtotal - ((FFLOAT) log10((double) lastCount)))\
                         * ((FFLOAT) lastCount) *\
                         job->bppWeight[unum]\
                         / job->NumBlocks[unum])\
                        )/19.26592);\
        } \
        else { \
	  entr_temp = ((FFLOAT) 0.0); \
        } \
	Entropy += entr_temp; \
	if (lastVal <= qvh) { \
	    SavedCounts[lastVal] = lastCount; \
	    SavedEntr[lastVal] = entr_temp; \
	} \
        lastCount = 0;\
        lastVal = quant;\
      }\
      lastCount += H(n).PlusCount[v];\
  }\
  if (lastCount > 0) {\
      entr_temp = \
           ((((logtotal - ((FFLOAT) log10((double) lastCount)))\
                 * ((FFLOAT) lastCount) *\
                 job->bppWeight[unum]\
                 / job->NumBlocks[unum])\
                )/19.26592);\
  } \
  else { \
      entr_temp = ((FFLOAT) 0.0); \
  } \
  Entropy += entr_temp; \
  if (lastVal <= qvh) { \
      SavedCounts[lastVal] = lastCount; \
      SavedEntr[lastVal] = entr_temp; \
  } \
  LastEntr = Entropy; \
  B = Entropy;\
}

#define SetBpp(B) \
{\
  int v; \
  FFLOAT Entropy;\
  int qv,qvneg,change; \
  Entropy = LastEntr; \
  change = 0; \
  v = t2-1; \
  qv = QuantizeDis(v, q); \
  qvneg = 0 - qv; \
  if ((v < H(n).MinusSize) && (H(n).MinusCount[v] > 0)) { \
    Entropy -= SavedEntr[qvneg]; \
    SavedCounts[qvneg] -= H(n).MinusCount[v]; \
    if (SavedCounts[qvneg] > 0) { \
      SavedEntr[qvneg] = \
         ((((logtotal - ((FFLOAT) log10((double) SavedCounts[qvneg])))\
                 * ((FFLOAT) SavedCounts[qvneg]) *\
                   job->bppWeight[unum]\
                 / job->NumBlocks[unum])\
                   )/19.26592);\
    } \
    else { \
      SavedEntr[qvneg] = ((FFLOAT) 0.0); \
    } \
    Entropy += SavedEntr[qvneg]; \
    SavedCounts[0] += H(n).MinusCount[v]; \
    change = 1; \
  } \
  if ((v < H(n).PlusSize) && (H(n).PlusCount[v] > 0)) { \
    Entropy -= SavedEntr[qv]; \
    SavedCounts[qv] -= H(n).PlusCount[v]; \
    if (SavedCounts[qv] > 0) { \
      SavedEntr[qv] = \
         ((((logtotal - ((FFLOAT) log10((double) SavedCounts[qv])))\
                 * ((FFLOAT) SavedCounts[qv]) *\
                   job->bppWeight [unum]\
                 / job->NumBlocks[unum])\
                   )/19.26592);\
    } \
    else { \
      SavedEntr[qv] = ((FFLOAT) 0.0); \
    } \
    Entropy += SavedEntr[qv]; \
    SavedCounts[0] += H(n).PlusCount[v]; \
    change = 1; \
  } \
  if (change) { \
    Entropy -= SavedEntr[0]; \
    SavedEntr[0] = \
         ((((logtotal - ((FFLOAT) log10((double) SavedCounts[0])))\
                 * ((FFLOAT) SavedCounts[0]) *\
                   job->bppWeight [unum]\
                 / job->NumBlocks[unum])\
                   )/19.26592);\
    Entropy += SavedEntr[0]; \
  } \
  LastEntr = Entropy; \
  B = Entropy; \
}

#define SetBppNoThresh(B) \
{\
  int v, lastVal,lastCount,quant;\
  FFLOAT Entropy;\
  Entropy = 0.0;\
  lastVal = -5000; lastCount = 0;\
  for (v=1-H(n).MinusSize;v<=0;v++) {\
        quant = QuantizeDis(v, q);\
        if (quant > lastVal ) {\
        if (lastCount > 0) {\
          Entropy = Entropy +\
                   ((((logtotal - ((FFLOAT) log10((double) lastCount)))\
                         * ((FFLOAT) lastCount) *\
                         job->bppWeight [unum]\
                         / job->NumBlocks[unum])\
                        )/19.26592);\
        }\
        lastCount = 0;\
        lastVal = quant;\
      }\
        lastCount += H(n).MinusCount[0-v];\
  }\
  for (v=0;v<H(n).PlusSize;v++) {\
        quant = QuantizeDis(v, q);\
        if (quant > lastVal ) {\
        if (lastCount > 0) {\
          Entropy = Entropy +\
                   ((((logtotal - ((FFLOAT) log10((double) lastCount)))\
                         * ((FFLOAT) lastCount) *\
                         job->bppWeight [unum]\
                         / job->NumBlocks[unum])\
                        )/19.26592);\
        }\
        lastCount = 0;\
        lastVal = quant;\
      }\
        lastCount += H(n).PlusCount[v];\
  }\
  if (lastCount > 0) {\
        Entropy = Entropy +\
           ((((logtotal - ((FFLOAT) log10((double) lastCount)))\
                 * ((FFLOAT) lastCount) *\
                 job->bppWeight [unum]\
                 / job->NumBlocks[unum])\
                )/19.26592);\
  }\
  B = Entropy;\
}

#define ToIntBpp(b) ((int) ((FFLOAT) ((b)*multiplier)+0.5))

extern void SetErrBpp(OptimJob *job, int unum)
{
    int n, q, q1, t2, t2lim;
    int  bppspan, leastintbpp,  indx, intbpp;
    int tspanplus1 = job->ThreshSpan + 1;
    FFLOAT leastbpp, mostbpp,  multiplier, bpp, error;
    FFLOAT ErrSavedPlus[QTABENTRYMAX+1], ErrSavedMinus[QTABENTRYMAX+1];
    FFLOAT LastErr;
    FFLOAT logtotal;
    int *SavedCounts, SavedCountBuff[2*QTABENTRYMAX];
    FFLOAT *SavedEntr, SavedEntrBuff[2*QTABENTRYMAX];
    FFLOAT LastEntr;



    multiplier = ((FFLOAT) job->opt_method.dp.BppScale);
    logtotal = ((FFLOAT) log10((double) job->NumBlocks[unum]));

    if (job->VerboseLevel > 2) fprintf(errfile,"\t\t");

    if (job->ThreshSpan > 0)
    {

        SavedCounts = &SavedCountBuff[QTABENTRYMAX];
        SavedEntr = &SavedEntrBuff[QTABENTRYMAX];

        for (n=0; n<64; n++)
        {

            q1 = job->MinTable[unum][n];

            if (job->MapQ) q = MapQentry(n,q1);
            else q = q1;

            SetBppNoThresh(mostbpp);
            bppspan = ToIntBpp(mostbpp) + 5;

            if (bppspan <
                    ((job->MaxTable[unum][n]-job->MinTable[unum][n]+1)*tspanplus1))
            {

                job->opt_method.dp.ErrEncodesRow[unum][n] = TRUE;

                /* allocate QforBpp and TforBpp */
                job->opt_method.dp.QforBpp[n] =
                    (Qentry *) calloc(1, sizeof(Qentry) * bppspan);
                job->opt_method.dp.TforBpp[n] =
                    (Qentry *) calloc(1, sizeof(Qentry) * bppspan);
                if ((!job->opt_method.dp.TforBpp[n]) ||
                        (!job->opt_method.dp.QforBpp[n]))
                {
                    FatalError("Could not allocate memory: SetErrBpp");
                }

                job->opt_method.dp.BppOffsetInEncoding[unum][n] = 0;
                job->opt_method.dp.BppMaxInEncoding[unum][n] = bppspan -1;
                for (intbpp = 0; intbpp < bppspan; intbpp++)
                {
                    job->opt_method.dp.Err[unum][n][intbpp] = INFINITY;
                }

                indx = job->MinTable[unum][n]*tspanplus1;
                for (q1=job->MinTable[unum][n]; q1<=job->MaxTable[unum][n]; q1++)
                {

                    if (job->MapQ) q = MapQentry(n,q1);
                    else q = q1;

                    SetBppPrepareThresh(bpp);
                    intbpp = ToIntBpp(bpp);
                    if (intbpp > job->opt_method.dp.BppMaxInEncoding[unum][n])
                        intbpp = job->opt_method.dp.BppMaxInEncoding[unum][n];
                    job->opt_method.dp.Bpp[unum][n][indx++] = intbpp;
                    SetErrPrepareThresh(error);
                    if (error < job->opt_method.dp.Err[unum][n][intbpp])
                    {
                        job->opt_method.dp.Err[unum][n][intbpp] = error;
                        job->opt_method.dp.QforBpp[n][intbpp] = q1;
                        job->opt_method.dp.TforBpp[n][intbpp] = 0;
                    }

                    t2lim = q + job->ThreshSpan;


                    for (t2 = q+1; t2 <= t2lim; t2++, indx++)
                    {

                        SetBpp(bpp);

                        intbpp = ToIntBpp(bpp);

                        if (intbpp > job->opt_method.dp.BppMaxInEncoding[unum][n])
                            intbpp = job->opt_method.dp.BppMaxInEncoding[unum][n];

                        job->opt_method.dp.Bpp[unum][n][indx] = intbpp;

                        SetErr(error);

                        if (error < job->opt_method.dp.Err[unum][n][intbpp])
                        {
                            job->opt_method.dp.Err[unum][n][intbpp] = error;
                            job->opt_method.dp.QforBpp[n][intbpp] = q1;
                            job->opt_method.dp.TforBpp[n][intbpp] = t2 - q;
                        }
                    }
                }

            }
            else
            {
                job->opt_method.dp.ErrEncodesRow[unum][n] = FALSE;

                indx = job->MinTable[unum][n]*tspanplus1;
                for (q1=job->MinTable[unum][n]; q1<=job->MaxTable[unum][n]; q1++)
                {

                    if (job->MapQ) q = MapQentry(n,q1);
                    else q = q1;

                    SetErrPrepareThresh(job->opt_method.dp.Err[unum][n][indx]);
                    SetBppPrepareThresh(bpp);


                    job->opt_method.dp.Bpp[unum][n][indx] = ToIntBpp(bpp);
                    indx++;

                    t2lim = q + job->ThreshSpan;

                    for (t2 = q+1; t2 <= t2lim; t2++,indx++)
                    {

                        SetErr(job->opt_method.dp.Err[unum][n][indx]);
                        SetBpp(bpp);

                        job->opt_method.dp.Bpp[unum][n][indx] = ToIntBpp(bpp);
                    }
                }
            }


            if (job->VerboseLevel > 2) fprintf(errfile,".");
        }
    }
    else
    {
        /* No thresholding */
        for (n=0; n<64; n++)
        {

            q1 = job->MaxTable[unum][n];

            if (job->MapQ) q = MapQentry(n,q1);
            else q = q1;

            SetBppNoThresh(leastbpp);

            leastintbpp = ToIntBpp(leastbpp);

            q1 = job->MinTable[unum][n];
            if (job->MapQ) q = MapQentry(n,q1);
            else q = q1;

            SetBppNoThresh(mostbpp);

            bppspan = ToIntBpp(mostbpp) - leastintbpp + 5;


            if (bppspan <
                    (job->MaxTable[unum][n]-job->MinTable[unum][n]+1))
            {

                job->opt_method.dp.ErrEncodesRow[unum][n] = TRUE;

                /* allocate QforBpp */
                job->opt_method.dp.QforBpp[n] =
                    (Qentry *) calloc(1, sizeof(Qentry) * bppspan);

                if (!job->opt_method.dp.QforBpp[n])
                {
                    FatalError("Could not allocate memory: SetErrBpp");
                }

                job->opt_method.dp.BppOffsetInEncoding[unum][n] = leastintbpp;
                job->opt_method.dp.BppMaxInEncoding[unum][n] = leastintbpp + bppspan -1;
                for (intbpp = 0; intbpp < bppspan; intbpp++)
                {
                    job->opt_method.dp.Err[unum][n][intbpp] = INFINITY;
                }

                for (q1=job->MinTable[unum][n]; q1<=job->MaxTable[unum][n]; q1++)
                {

                    if (job->MapQ) q = MapQentry(n,q1);
                    else q = q1;

                    SetBppNoThresh(bpp);

                    intbpp = ToIntBpp(bpp);

                    if (intbpp < leastintbpp) intbpp = leastintbpp;
                    if (intbpp > job->opt_method.dp.BppMaxInEncoding[unum][n])
                        intbpp = job->opt_method.dp.BppMaxInEncoding[unum][n];

                    job->opt_method.dp.Bpp[unum][n][q1] = intbpp;

                    intbpp -= leastintbpp;


                    SetErrNoThresh(error);


                    if (error < job->opt_method.dp.Err[unum][n][intbpp])
                    {
                        job->opt_method.dp.Err[unum][n][intbpp] = error;
                        job->opt_method.dp.QforBpp[n][intbpp] = q1;
                    }
                }

            }
            else
            {
                job->opt_method.dp.ErrEncodesRow[unum][n] = FALSE;

                for (q1=job->MinTable[unum][n]; q1<=job->MaxTable[unum][n]; q1++)
                {

                    if (job->MapQ) q = MapQentry(n,q1);
                    else q = q1;


                    SetErrNoThresh(job->opt_method.dp.Err[unum][n][q1]);
                    SetBppNoThresh(bpp);

                    job->opt_method.dp.Bpp[unum][n][q1] = ToIntBpp(bpp);
                }
            }


            if (job->VerboseLevel > 2) fprintf(errfile,".");
        }
    }

    if (job->VerboseLevel > 2)
    {
        fprintf(errfile,"\n");
        fflush(errfile);
    }

}




