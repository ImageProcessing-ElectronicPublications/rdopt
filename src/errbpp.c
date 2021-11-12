/*
** Copyright 1995 by Viresh Ratnakar, Miron Livny
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
        fprintf(stderr,"\t\tQuantTable Min/Max entries:\n");
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
            if ((n%8)==0) fprintf(stderr,"\t\t");
#if (QTABBITS==8)
            fprintf(stderr,"%3d/%3d ",
#elif (QTABBITS==10)
            fprintf(stderr,"%4d/%4d ",
#elif (QTABBITS==12)
            fprintf(stderr,"%4d/%4d ",
#elif (QTABBITS==16)
            fprintf(stderr,"%5d/%5d ",
#else
            /** force syntax error **/
            QTABBITS must be 8,12, or 16
#endif
                    job->MinTable[unum][n],job->MaxTable[unum][n]);
            if ((n%8)==7) fprintf(stderr,"\n");
        }

        if (job->MapQ)
        {
            job->MaxTable[unum][n] = UnMapQentry(n, job->MaxTable[unum][n]);
            job->MinTable[unum][n] = UnMapQentry(n, job->MinTable[unum][n]);
        }
        totalentries += job->MaxTable[unum][n];
    }

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
        job->Err[n] = fptr + offset;
        job->Bpp[unum][n] = iptr + offset;
        offset += (job->MaxTable[unum][n]+1);
    }
}


extern void SetErr(OptimJob *job, int unum)
{
    int i,n,q,intquant, q1;
    FFLOAT orig, err, ssq;

    if (job->VerboseLevel > 2) fprintf(stderr,"\t\t");
    for (n=0; n<64; n++)
    {
        for (q1=job->MinTable[unum][n]; q1<=job->MaxTable[unum][n]; q1++)
        {
            if (job->MapQ) q = MapQentry(n,q1);
            else q = q1;
            job->Err[n][q1] = 0.0;
            for (i=H(n).PlusSize-1; i>=0; i--)
            {
                orig = UnDiscretize(i,TRUE);
                intquant = QuantizeDis(i,q);
                err = orig - RealUnQuantize(intquant,q);
                job->Err[n][q1] += (err*err*((FFLOAT) H(n).PlusCount[i]));
            }
            for (i=1-H(n).MinusSize; i<=0; i++)
            {
                orig = UnDiscretize(i,FALSE);
                intquant = QuantizeDis(i,q);
                err = orig - RealUnQuantize(intquant,q);
                job->Err[n][q1] += (err*err*((FFLOAT) H(n).MinusCount[0-i]));
            }
            if (job->WeightedCoefs) job->Err[n][q1] *= job->CoefWeights[n];
        }
        if (job->VerboseLevel > 2) fprintf(stderr,".");
    }
    if (job->VerboseLevel > 2) fprintf(stderr,"\n");
}

#define ToIntBpp(b) ((int) ((FFLOAT) ((b)*multiplier)+0.5))

extern FFLOAT GetBpp(OptimJob *job, int unum, int n, int q)
{
    int v, lastVal,lastCount,quant, q1;
    FFLOAT Entropy;
    static int lastunum = -1;
    static FFLOAT logtotal;

    if (job->MapQ) q1 = MapQentry(n,q);
    else q1 = q;
    if (unum > lastunum)
    {
        logtotal = ((FFLOAT) log10((double) job->NumBlocks[unum]));
        lastunum = unum;
    }

    Entropy = 0.0;
    lastVal = -5000;
    lastCount = 0;
    for (v=1-H(n).MinusSize; v<=0; v++)
    {
        quant = QuantizeDis(v, q1);
        if (quant > lastVal )
        {
            if (lastCount > 0)
            {
                Entropy = Entropy +
                          ((((logtotal - ((FFLOAT) log10((double) lastCount)))
                             * ((FFLOAT) lastCount) *
                             job->bppWeight [unum]
                             / job->NumBlocks[unum])
                           )/0.30103);

            }
            lastCount = 0;
            lastVal = quant;
        }
        lastCount += H(n).MinusCount[0-v];
    }
    for (v=0; v<H(n).PlusSize; v++)
    {
        quant = QuantizeDis(v, q1);
        if (quant > lastVal )
        {
            if (lastCount > 0)
            {
                Entropy = Entropy +
                          ((((logtotal - ((FFLOAT) log10((double) lastCount)))
                             * ((FFLOAT) lastCount) *
                             job->bppWeight [unum]
                             / job->NumBlocks[unum])
                           )/0.30103);

            }
            lastCount = 0;
            lastVal = quant;
        }
        lastCount += H(n).PlusCount[v];
    }
    if (lastCount > 0)
    {
        Entropy = Entropy +
                  ((((logtotal - ((FFLOAT) log10((double) lastCount)))
                     * ((FFLOAT) lastCount) *
                     job->bppWeight [unum]
                     / job->NumBlocks[unum])
                   )/0.30103);
    }
    Entropy /= 64.0;
    return(Entropy);
}


extern void SetBpp(OptimJob *job, int unum)
{
    int n,v,q,lastVal,lastCount,quant, q1;
    FFLOAT multiplier,Entropy,logtotal;

    multiplier = ((FFLOAT) job->BppScale)/64.0;
    logtotal = ((FFLOAT) log10((double) job->NumBlocks[unum]));

    if (job->VerboseLevel > 2) fprintf(stderr,"\t\t");
    for (n=0; n<64; n++)
    {
        for (q1=job->MinTable[unum][n]; q1<=job->MaxTable[unum][n]; q1++)
        {
            if (job->MapQ) q = MapQentry(n,q1);
            else q = q1;
            Entropy = 0.0;
            lastVal = -5000;
            lastCount = 0;
            for (v=1-H(n).MinusSize; v<=0; v++)
            {
                quant = QuantizeDis(v, q);
                if (quant > lastVal )
                {
                    if (lastCount > 0)
                    {
                        Entropy = Entropy +
                                  ((((logtotal - ((FFLOAT) log10((double) lastCount)))
                                     * ((FFLOAT) lastCount) *
                                     job->bppWeight [unum]
                                     / job->NumBlocks[unum])
                                   )/0.30103);

                    }
                    lastCount = 0;
                    lastVal = quant;
                }
                lastCount += H(n).MinusCount[0-v];
            }
            for (v=0; v<H(n).PlusSize; v++)
            {
                quant = QuantizeDis(v, q);
                if (quant > lastVal )
                {
                    if (lastCount > 0)
                    {
                        Entropy = Entropy +
                                  ((((logtotal - ((FFLOAT) log10((double) lastCount)))
                                     * ((FFLOAT) lastCount) *
                                     job->bppWeight [unum]
                                     / job->NumBlocks[unum])
                                   )/0.30103);

                    }
                    lastCount = 0;
                    lastVal = quant;
                }
                lastCount += H(n).PlusCount[v];
            }
            if (lastCount > 0)
            {
                Entropy = Entropy +
                          ((((logtotal - ((FFLOAT) log10((double) lastCount)))
                             * ((FFLOAT) lastCount) *
                             job->bppWeight [unum]
                             / job->NumBlocks[unum])
                           )/0.30103);
            }
            job->Bpp[unum][n][q1] = ToIntBpp(Entropy);

        }
        if (job->VerboseLevel > 2) fprintf(stderr,".");
    }
    if (job->VerboseLevel > 2) fprintf(stderr,"\n");
}

extern void TranslateBppRange(OptimJob *job, int unum)
{
    int n,l,u,curr,ltemp,utemp, u1, l1, ltemp1, utemp1;
    FFLOAT bl,bu,m,bltemp,butemp,bcurr;
    boolean done;
    FFLOAT x,y;

#define BPPZERO ((FFLOAT) 1.0/job->NumBlocks[unum])

    if (job->VerboseLevel > 3)
    {
        fprintf(stderr,"\t\tLeast:Most bits per pixel specified:\n");
        for (n=0; n<64; n++)
        {
            if ((n%8)==0) fprintf(stderr,"\t\t  ");
            fprintf(stderr,"%.3lf:%.3lf ",job->BppRange[unum][n][0],
                    job->BppRange[unum][n][1]);
            if ((n%8)==7) fprintf(stderr,"\n");
        }
    }
    if (job->VerboseLevel > 2) fprintf(stderr,"\t\t");


    if (!job->MapQ)
    {

        for (n=0; n<64; n++)
        {

            x = job->BppRange[unum][n][0];
            y = job->BppRange[unum][n][1];


            if (y <= BPPZERO)
            {
                job->MinTable[unum][n] = QTABENTRYMAX;
                job->MaxTable[unum][n] = QTABENTRYMAX;
                if (job->VerboseLevel > 2) fprintf(stderr,".");
                continue;
            }

            l = 1;

            bl = GetBpp(job, unum, n, l);

            if (bl < x)
            {
                job->MinTable[unum][n] = l;
                job->MaxTable[unum][n] = l;
                if (job->VerboseLevel > 2) fprintf(stderr,".");
                continue;
            }

            u = QTABENTRYMAX;

            bu = GetBpp(job, unum, n, u);

            if (bu > y)
            {
                job->MinTable[unum][n] = u;
                job->MaxTable[unum][n] = u;
                if (job->VerboseLevel > 2) fprintf(stderr,".");
                continue;
            }

            if ((bu >= x) && (bl <= y))
            {
                job->MinTable[unum][n] = l;
                job->MaxTable[unum][n] = u;
                if (job->VerboseLevel > 2) fprintf(stderr,".");
                continue;
            }

            if (bu >= x)
            {
                job->MaxTable[unum][n] = u;
            }
            else
            {
                done = FALSE;
                utemp = u;
                butemp = bu;
                ltemp = l;
                bltemp = bl;
                while (!done)
                {
                    while (butemp < x)
                    {
                        bu = butemp;
                        u = utemp;
                        curr = (utemp+ltemp)/2;
                        while (curr >= utemp) curr--;
                        while (curr < ltemp) curr++;
                        utemp = curr;
                        butemp = GetBpp(job, unum, n, utemp);
                    }
                    if (GetBpp(job, unum, n, utemp+1) <= x)
                    {
                        u = utemp;
                        bu = butemp;
                        done = TRUE;
                    }
                    else
                    {
                        bltemp = butemp;
                        ltemp = utemp;
                        butemp = bu;
                        utemp = u;
                    }
                }
                job->MaxTable[unum][n] = u;
            }


            if (bu >= y)
            {
                job->MinTable[unum][n] = u;
            }
            else if (bl <= y)
            {
                job->MinTable[unum][n] = l;
            }
            else
            {
                done = FALSE;
                utemp = u;
                butemp = bu;
                ltemp = l;
                bltemp = bl;
                while (!done)
                {
                    while (butemp < y)
                    {
                        bu = butemp;
                        u = utemp;
                        curr = (utemp+ltemp)/2;
                        while (curr >= utemp) curr--;
                        while (curr < ltemp) curr++;
                        utemp = curr;
                        butemp = GetBpp(job, unum, n, utemp);
                    }
                    if (GetBpp(job, unum, n, utemp+1) <= y)
                    {
                        u = utemp;
                        bu = butemp;
                        done = TRUE;
                    }
                    else
                    {
                        bltemp = butemp;
                        ltemp = utemp;
                        butemp = bu;
                        utemp = u;
                    }
                }
                job->MinTable[unum][n] = u;
            }
            if (job->VerboseLevel > 2) fprintf(stderr,".");
        }
        if (job->VerboseLevel > 2) fprintf(stderr,"\n");
    }

    else
    {

        for (n=0; n<64; n++)
        {

            x = job->BppRange[unum][n][0];
            y = job->BppRange[unum][n][1];


            if (y <= BPPZERO)
            {
                job->MinTable[unum][n] = QTABENTRYMAX;
                job->MaxTable[unum][n] = QTABENTRYMAX;
                if (job->VerboseLevel > 2) fprintf(stderr,".");
                continue;
            }

            l = 1;
            l1 = MapQentry(n,l);

            bl = GetBpp(job, unum, n, l);

            if (bl < x)
            {
                job->MinTable[unum][n] = l1;
                job->MaxTable[unum][n] = l1;
                if (job->VerboseLevel > 2) fprintf(stderr,".");
                continue;
            }

            u = 255;
            u1 = MapQentry(n,u);

            bu = GetBpp(job, unum, n, u);

            if (bu > y)
            {
                job->MinTable[unum][n] = u1;
                job->MaxTable[unum][n] = u1;
                if (job->VerboseLevel > 2) fprintf(stderr,".");
                continue;
            }

            if ((bu >= x) && (bl <= y))
            {
                job->MinTable[unum][n] = l1;
                job->MaxTable[unum][n] = u1;
                if (job->VerboseLevel > 2) fprintf(stderr,".");
                continue;
            }

            if (bu >= x)
            {
                job->MaxTable[unum][n] = u1;
            }
            else
            {
                done = FALSE;
                utemp = u;
                butemp = bu;
                utemp1 = u1;
                ltemp = l;
                bltemp = bl;
                ltemp1 = l1;
                while (!done)
                {
                    while (butemp < x)
                    {
                        bu = butemp;
                        u = utemp;
                        u1 = utemp1;
                        curr = (utemp+ltemp)/2;
                        while (curr >= utemp) curr--;
                        while (curr < ltemp) curr++;
                        utemp = curr;
                        utemp1 = MapQentry(n,utemp);
                        butemp = GetBpp(job, unum, n, utemp);
                    }
                    if (GetBpp(job, unum, n, utemp+1) <= x)
                    {
                        u = utemp;
                        u1 = utemp1;
                        bu = butemp;
                        done = TRUE;
                    }
                    else
                    {
                        bltemp = butemp;
                        ltemp = utemp;
                        ltemp1 = utemp1;
                        butemp = bu;
                        utemp = u;
                        utemp1 = u1;
                    }
                }
                job->MaxTable[unum][n] = u1;
            }


            if (bu >= y)
            {
                job->MinTable[unum][n] = u1;
            }
            else if (bl <= y)
            {
                job->MinTable[unum][n] = l1;
            }
            else
            {
                done = FALSE;
                utemp = u;
                butemp = bu;
                utemp1 = u1;
                ltemp = l;
                bltemp = bl;
                ltemp1 = l1;
                while (!done)
                {
                    while (butemp < y)
                    {
                        bu = butemp;
                        u = utemp;
                        u1 = utemp1;
                        curr = (utemp+ltemp)/2;
                        while (curr >= utemp) curr--;
                        while (curr < ltemp) curr++;
                        utemp = curr;
                        utemp1 = MapQentry(n,utemp);
                        butemp = GetBpp(job, unum, n, utemp);
                    }
                    if (GetBpp(job, unum, n, utemp+1) <= y)
                    {
                        u = utemp;
                        u1 = utemp1;
                        bu = butemp;
                        done = TRUE;
                    }
                    else
                    {
                        bltemp = butemp;
                        ltemp = utemp;
                        ltemp1 = utemp1;
                        butemp = bu;
                        utemp = u;
                        utemp1 = u1;
                    }
                }
                job->MinTable[unum][n] = u1;
            }
            if (job->VerboseLevel > 2) fprintf(stderr,".");
        }
        if (job->VerboseLevel > 2) fprintf(stderr,"\n");
    }
}
