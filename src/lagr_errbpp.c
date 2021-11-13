
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

extern void lagrPrepareForErrBpp(OptimJob *job, int unum)
{
    int i,n, offset, totalentries, max;
    FFLOAT *fptr1, *fptr2, *fptr3;
    int *iptr;

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

    if ((fptr1 = (FFLOAT *) calloc(1,totalentries*sizeof(FFLOAT)))==NULL)
    {
        FatalError("PrepareForErrBpp out of memory");
    }
    if ((fptr2 = (FFLOAT *) calloc(1,totalentries*sizeof(FFLOAT)))==NULL)
    {
        FatalError("PrepareForErrBpp out of memory");
    }
    if ((fptr3 = (FFLOAT *) calloc(1,totalentries*sizeof(FFLOAT)))==NULL)
    {
        FatalError("PrepareForErrBpp out of memory");
    }
    if ((iptr = (int *) calloc(1,totalentries*sizeof(int)))==NULL)
    {
        FatalError("PrepareForErrBpp out of memory");
    }

    offset = 0;
    for (n=0; n<64; n++)
    {
        job->opt_method.lagr.Err[unum][n] = fptr1 + offset;
        job->opt_method.lagr.Bpp[unum][n] = fptr2 + offset;
        job->opt_method.lagr.Lambda[unum][n] = fptr3 + offset;
        job->opt_method.lagr.SortedIndex[unum][n] = iptr + offset;

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
      E = ans*job->CompWeights[unum]/job->opt_method.lagr.ErrScale; \
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
      E *= (job->CompWeights[unum]/job->opt_method.lagr.ErrScale); \
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
      E = (ans*job->CompWeights[unum])/job->opt_method.lagr.ErrScale; \
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


extern void lagrSetErrBpp(OptimJob *job, int unum)
{
    int n, q, q1, t2, t2lim;
    int  indx;
    int tspanplus1 = job->ThreshSpan + 1;
    FFLOAT ErrSavedPlus[QTABENTRYMAX+1], ErrSavedMinus[QTABENTRYMAX+1];
    FFLOAT LastErr;
    FFLOAT logtotal;
    int *SavedCounts, SavedCountBuff[2*QTABENTRYMAX];
    FFLOAT *SavedEntr, SavedEntrBuff[2*QTABENTRYMAX];
    FFLOAT LastEntr;



    logtotal = ((FFLOAT) log10((double) job->NumBlocks[unum]));

    if (job->VerboseLevel > 2) fprintf(errfile,"\t\t");

    if (job->ThreshSpan > 0)
    {

        SavedCounts = &SavedCountBuff[QTABENTRYMAX];
        SavedEntr = &SavedEntrBuff[QTABENTRYMAX];

        for (n=0; n<64; n++)
        {

            indx = job->MinTable[unum][n]*tspanplus1;
            for (q1=job->MinTable[unum][n]; q1<=job->MaxTable[unum][n]; q1++)
            {

                if (job->MapQ) q = MapQentry(n,q1);
                else q = q1;

                SetErrPrepareThresh(job->opt_method.lagr.Err[unum][n][indx]);
                SetBppPrepareThresh(job->opt_method.lagr.Bpp[unum][n][indx]);


                indx++;

                t2lim = q + job->ThreshSpan;

                for (t2 = q+1; t2 <= t2lim; t2++,indx++)
                {

                    SetErr(job->opt_method.lagr.Err[unum][n][indx]);
                    SetBpp(job->opt_method.lagr.Bpp[unum][n][indx]);
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

            for (q1=job->MinTable[unum][n]; q1<=job->MaxTable[unum][n]; q1++)
            {

                if (job->MapQ) q = MapQentry(n,q1);
                else q = q1;


                SetErrNoThresh(job->opt_method.lagr.Err[unum][n][q1]);
                SetBppNoThresh(job->opt_method.lagr.Bpp[unum][n][q1]);
            }

            if (job->VerboseLevel > 2) fprintf(errfile,".");
        }
    }

#ifdef TESTING
#include "tests.c"
#endif


    if (job->VerboseLevel > 2)
    {
        fprintf(errfile,"\n");
        fflush(errfile);
    }

}


#define EQUALITY_DELTA ((FFLOAT) 1e-30)

#define IS_RIGHT_TURN(x0,y0,x1,y1,x2,y2) \
   (((((x2)-(x0))*((y1)-(y0))) - (((x1)-(x0))*((y2)-(y0)))) >= 0)


#define PUSH(p,S) \
   stacktop++; \
   (S)[stacktop] = (p)

#define POP(S) stacktop--

#define TOP(S) ((S)[stacktop])

#define NEXT_TO_TOP(S) ((S)[stacktop-1])

/* comparison function for qsort */
FFLOAT *tosort;
static int comp_func(int *i, int *j)
{
    if (tosort[*i] >= tosort[*j]) return(-1);
    else return(1);
}

extern void lagrSortErrBpp(OptimJob *job, int unum)
{
    int tspanplus1 = job->ThreshSpan + 1;
    int indx, i, j, n, tot, temp, done;
    int *idxptr,*stackptr;
    int *tempIdx;
    FFLOAT rdiff, rlastval, ddiff, dlastval;
    int stacktop, orig_tot, lastwhere, nextpos;


    if ((tempIdx = (int *)
                   calloc(1,sizeof(int)*(QTABENTRYMAX+1)*tspanplus1)) == NULL)
    {
        FatalError("lagrSortErrBpp out of memory");
    }

#define optstr job->opt_method.lagr

    for (n=0; n<64; n++)
    {

        tot = (job->MaxTable[unum][n] - job->MinTable[unum][n] + 1)*
              tspanplus1;

        indx = job->MinTable[unum][n]*tspanplus1;
        idxptr = tempIdx;

        for (i=0,indx=job->MinTable[unum][n]*tspanplus1; i<tot; i++,indx++)
        {
            idxptr[i] = indx;
        }

        if (job->VerboseLevel > 2)
        {
            fprintf(errfile,"\t\tcoefficient #%d ..",n);
            fflush(errfile);
        }
        if (job->ThreshSpan == 0)
        {
            /** bubble sort to ensure decreasing rate **/
            /** "bubble sort" 'cos the array must already be nearly sorted **/
            done = 0;
            while (!done)
            {
                done = 1;
                for (i=1; i<tot; i++)
                {
                    if (optstr.Bpp[unum][n][idxptr[i-1]] <
                            optstr.Bpp[unum][n][idxptr[i]])
                    {
                        temp = idxptr[i-1];
                        idxptr[i-1] = idxptr[i];
                        idxptr[i] = temp;
                        done = 0;
                    }
                }
            }
        }
        else
        {
            /* quick sort */
            tosort = optstr.Bpp[unum][n];
            qsort((char *) idxptr, tot, sizeof(int), comp_func);
        }

        orig_tot = tot;

        /** ensure "strictly" decreasing rate and "strictly" increasing distortion **/
        rlastval = optstr.Bpp[unum][n][idxptr[0]];
        dlastval = optstr.Err[unum][n][idxptr[0]];
        lastwhere = idxptr[0];
        nextpos = 0;
        for (i=1; i<tot; i++)
        {
            if (lastwhere != -1)
            {
                rdiff = rlastval - optstr.Bpp[unum][n][idxptr[i]];
                if (rdiff <= EQUALITY_DELTA)
                {
                    if (optstr.Err[unum][n][idxptr[i]] < dlastval )
                    {
                        lastwhere = idxptr[i];
                        rlastval = optstr.Bpp[unum][n][idxptr[i]];
                        dlastval = optstr.Err[unum][n][idxptr[i]];
                    }
                }
                else
                {
                    idxptr[nextpos++] = lastwhere;
                    lastwhere = -1;
                }
            }
            else
            {
                lastwhere = idxptr[i];
                rlastval = optstr.Bpp[unum][n][idxptr[i]];
                dlastval = optstr.Err[unum][n][idxptr[i]];
            }
        }
        if (lastwhere != -1) idxptr[nextpos++] = lastwhere;

        tot = nextpos;


        rlastval = optstr.Bpp[unum][n][idxptr[0]];
        dlastval = optstr.Err[unum][n][idxptr[0]];
        lastwhere = idxptr[0];
        nextpos = 0;
        for (i=1; i<tot; i++)
        {
            if (lastwhere != -1)
            {
                ddiff = optstr.Err[unum][n][idxptr[i]] - dlastval;
                if (ddiff <= EQUALITY_DELTA)
                {
                    lastwhere = idxptr[i];
                    rlastval = optstr.Bpp[unum][n][idxptr[i]];
                    dlastval = optstr.Err[unum][n][idxptr[i]];
                }
                else
                {
                    idxptr[nextpos++] = lastwhere;
                    lastwhere = -1;
                }
            }
            else
            {
                lastwhere = idxptr[i];
                rlastval = optstr.Bpp[unum][n][idxptr[i]];
                dlastval = optstr.Err[unum][n][idxptr[i]];
            }
        }
        if (lastwhere != -1) idxptr[nextpos++] = lastwhere;

        tot = nextpos;

        if (job->VerboseLevel > 2)
        {
            fprintf(errfile,".. %d of %d R-D points\n",tot, orig_tot);
            fflush(errfile);
        }

#define Xval(i) job->opt_method.lagr.Err[unum][n][(i)]
#define Yval(i) job->opt_method.lagr.Bpp[unum][n][(i)]

        /* remove points not on the "convex hull" */
        stackptr = optstr.SortedIndex[unum][n];
        if (tot >= 3)
        {
            stacktop = -1;
            PUSH(idxptr[0],stackptr);
            PUSH(idxptr[1],stackptr);
            PUSH(idxptr[2],stackptr);
            for (i=3; i<tot; i++)
            {
                while ((stacktop > 0) &&
                        (IS_RIGHT_TURN(Xval(NEXT_TO_TOP(stackptr)),
                                       Yval(NEXT_TO_TOP(stackptr)),
                                       Xval(TOP(stackptr)),
                                       Yval(TOP(stackptr)),
                                       Xval(idxptr[i]),
                                       Yval(idxptr[i]))))
                {
                    POP(stackptr);
                }
                PUSH(idxptr[i],stackptr);
            }
            tot = optstr.IndexEntries[unum][n] = stacktop + 1;
        }
        else
        {
            for (i=0; i< tot; i++)
            {
                stackptr[i] = idxptr[i];
            }
            optstr.IndexEntries[unum][n] = tot;
        }

        idxptr = stackptr;


        tot--;

        /* calculate slopes */
        for (i=0; i<tot; i++)
        {
            optstr.Lambda[unum][n][i] =
                (optstr.Bpp[unum][n][idxptr[i]] -
                 optstr.Bpp[unum][n][idxptr[i+1]]) /
                (optstr.Err[unum][n][idxptr[i+1]] -
                 optstr.Err[unum][n][idxptr[i]]);
        }
        optstr.Lambda[unum][n][tot] = ((FFLOAT) 0.0);

        if (job->VerboseLevel > 3)
        {
            fprintf(errfile,"\t\t\tBPP range %6.5lf to %6.5lf over %d R-D points\n", optstr.Bpp[unum][n][idxptr[tot]], optstr.Bpp[unum][n][idxptr[0]],tot+1);
            fflush(errfile);
        }


    }

    free(tempIdx);


}

extern void lagrEpilogue(OptimJob *job)
{
    int unum, i, j, n;
    FFLOAT err, psnr, bpp, bincr;
    FFLOAT lambda, askingbpp, delta;
    FILE *fp;
    int tab[MAXCOMPONENTS][64];
    FFLOAT actualBpp, TargetBpp, TargetErr;
    int Q[MAXCOMPONENTS][64], T[MAXCOMPONENTS][64],
        QandT[MAXCOMPONENTS][64];
    FFLOAT ErrDist[MAXCOMPONENTS][64], BppDist[MAXCOMPONENTS][64];
    char qfname[STRLENMAX];
    FFLOAT lambdahigh,lambdalow;

#define optstr job->opt_method.lagr

    if (job->VerboseLevel && (job->NumTables > 1))
    {
        /* for compatibility of output # lines with dyn prog */
        fprintf(errfile,"Postprocessing...\n");
        fflush(errfile);
        if (job->VerboseLevel > 1)
        {
            n = 2; /* # of blank lines */
        }
        else
        {
            n = 1;
        }
        for (unum=0; unum < job->NumTables; unum++)
        {
            for (i=0; i<n; i++) fprintf(errfile,"\t\n");
        }
        fflush(errfile);
    }

    /**** compute LogSignalSq ****/
    job->LogTotalSignalSq = 0.0;
    for (unum=0; unum < job->NumTables; unum++)
    {
        job->LogTotalSignalSq += job->LogSignalSq[unum];
        job->LogSignalSq[unum] = ((FFLOAT) log10((double) job->LogSignalSq[unum]));
    }
    job->LogTotalSignalSq = ((FFLOAT) log10((double) job->LogTotalSignalSq));
    optstr.LambdaMax = ((FFLOAT) 0.0);
    optstr.BppMax = ((FFLOAT) 0.0);
    optstr.BppMin = ((FFLOAT) 0.0);
    optstr.ErrMax = ((FFLOAT) 0.0);
    optstr.ErrMin = ((FFLOAT) 0.0);
    for (unum = 0; unum < job->NumTables; unum++)
    {
        for (n = 0; n < 64; n++)
        {
            if (optstr.Lambda[unum][n][0] > optstr.LambdaMax)
            {
                optstr.LambdaMax = optstr.Lambda[unum][n][0];
            }
            optstr.BppMax += optstr.Bpp[unum][n][optstr.SortedIndex[unum][n][0]];
            optstr.ErrMin += optstr.Err[unum][n][optstr.SortedIndex[unum][n][0]];
            optstr.BppMin += optstr.Bpp[unum][n][optstr.SortedIndex[unum][n][optstr.IndexEntries[unum][n]-1]];
            optstr.ErrMax += optstr.Err[unum][n][optstr.SortedIndex[unum][n][optstr.IndexEntries[unum][n]-1]];
        }
    }

    optstr.LambdaMax += ((FFLOAT) 0.1);

    optstr.BppHook[0] = optstr.BppMax;
    optstr.LambdaHook[0] = optstr.LambdaMax;
    optstr.ErrHook[0] = optstr.ErrMin;

    optstr.BppHook[optstr.NumHooks-1] = optstr.BppMin;
    optstr.LambdaHook[optstr.NumHooks-1] = ((FFLOAT) 0.0);
    optstr.ErrHook[optstr.NumHooks-1] = optstr.ErrMax;


    bincr = (optstr.BppMax - optstr.BppMin)/((FFLOAT) optstr.NumHooks);
    delta = bincr/((FFLOAT) 10.0);
    askingbpp = optstr.BppMax - bincr;
    for (i = 1; i < (optstr.NumHooks-1); i++, askingbpp -= bincr)
    {
        bpp = askingbpp;
        lagrFindRate(job, &bpp, &err, delta, 0.0, optstr.LambdaHook[i-1], &lambda, tab);
        optstr.BppHook[i] = bpp;
        optstr.LambdaHook[i] = lambda;
        optstr.ErrHook[i] = err;
    }

    if (job->useCorrection)
    {
        if (job->TheImage.ImFileName[0] == '\0')
        {
            fprintf(errfile,"Ignoring -correct: image read off stdin\n");
            job->useCorrection = FALSE;
            goto afterCorrection;
        }
        TargetBpp = job->correctionBpp;
        i = 0;
        while ((i < optstr.NumHooks) && (optstr.BppHook[i] >= TargetBpp))
            i++;

        if (i == 0)
        {
            lambdalow = lambdahigh = optstr.LambdaMax;
        }
        else if (i == optstr.NumHooks)
        {
            lambdalow = lambdahigh = ((FFLOAT) 0.0);
        }
        else
        {
            lambdahigh = optstr.LambdaHook[i-1];
            lambdalow = optstr.LambdaHook[i];
        }

        lagrFindRate(job, &TargetBpp, &TargetErr, 0.01,
                     lambdalow, lambdahigh, &lambda, QandT);

        sprintf(qfname,"/tmp/rdopt.%s.correctionQ\0",
                FileNameTail(job->TheImage.ImFileName));
        if ((fp = fopen(qfname,"w")) == NULL)
        {
            fprintf(errfile,"Could not open %s\n",qfname);
            job->useCorrection = FALSE;
            goto afterCorrection;
        }
        if (job->NumTables==1)
        {
            lagrRecoverQandT(job,QandT,Q,T,BppDist,ErrDist);
            if (job->ThreshSpan > 0)
            {
                WriteIntTable(fp, Q[0], "");
                WriteThreshTable(fp, T[0], "#T ");
            }
            else
            {
                WriteIntTable(fp, Q[0], "");
            }
        }
        else
        {
            lagrRecoverQandT(job,QandT,Q,T,BppDist,ErrDist);
            for (unum=0; unum < job->NumTables; unum++)
            {
                if (job->ThreshSpan > 0)
                {
                    WriteIntTable(fp, Q[unum], "");
                    WriteThreshTable(fp, T[unum], "#T ");
                }
                else
                {
                    WriteIntTable(fp, Q[unum], "");
                }
            }
        }
        fclose(fp);
        if (GetActualBpp(job, qfname, &actualBpp))
        {
            job->correctionBpp = TargetBpp;
            job->addToTableBpp = actualBpp - TargetBpp;
        }
        else
        {
            job->useCorrection = FALSE;
        }
        remove(qfname);
    }


afterCorrection:
    if (!job->Silent)
    {
        /** print range **/
        psnr = (job->LogTotalPeakSq -
                ((FFLOAT) log10((double) (optstr.ErrMax*optstr.ErrScale)))) *10.0;
        fprintf(errfile," Range: bpp (psnr) %lf (%lf) to ",
                (double) optstr.BppMin +
                ((double) (job->useCorrection?job->addToTableBpp:0)),
                (double) psnr);
        psnr = (job->LogTotalPeakSq -
                ((FFLOAT) log10((double) (optstr.ErrMin*optstr.ErrScale)))) *10.0;
        fprintf(errfile,"  %lf (%lf)\n",
                (double) optstr.BppMax +
                ((double) (job->useCorrection?job->addToTableBpp:0)),
                (double) psnr);
        fflush(errfile);
    }

    if (job->DumpPlot)
    {
        if (!strcmp(job->PlotFileName,"-")) fp = stdout;
        else if ((fp = fopen(job->PlotFileName,"w")) == NULL)
        {
            fprintf(errfile,"Could not open plot file.. ignoring\n");
            return;
            fflush(errfile);
        }

        fprintf(fp,"# BPP -- PSNR plot generated by RDOPT\n");
        if (strcmp(job->PlotFileName,"-")) DumpJobChars(job,fp);

        bpp = job->PlotBppMax;
        if (bpp > optstr.BppMax) bpp = optstr.BppMax;
        bincr = (bpp - optstr.BppMin)/((FFLOAT) job->PlotPoints);
        bpp = optstr.BppMin;
        for (i=0; i<= job->PlotPoints; i++, bpp += bincr)
        {
            TargetBpp = bpp;

            j = 0;
            while ((j < optstr.NumHooks) && (optstr.BppHook[j] >= TargetBpp))
                j++;

            if (j == 0)
            {
                lambdalow = lambdahigh = optstr.LambdaMax;
            }
            else if (j == optstr.NumHooks)
            {
                lambdalow = lambdahigh = ((FFLOAT) 0.0);
            }
            else
            {
                lambdahigh = optstr.LambdaHook[j-1];
                lambdalow = optstr.LambdaHook[j];
            }

            lagrFindRate(job, &TargetBpp, &TargetErr, 0.01,
                         lambdalow, lambdahigh, &lambda, QandT);

            psnr = (job->LogTotalPeakSq -
                    ((FFLOAT) log10((double) (TargetErr*optstr.ErrScale)))) *10.0;
            fprintf(fp,"%lf   %lf\n", (double) TargetBpp +
                    ((double) (job->useCorrection?job->addToTableBpp:0)),
                    (double) psnr);
        }

        if (strcmp(job->PlotFileName,"-")) fclose(fp);
    }
}

