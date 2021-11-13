
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



#define C(x) (((x) == 0) ? ((FFLOAT) M_SQRT1_2) : ((FFLOAT) 1.0))

#include "CosTable.h"

static int IJtoWhere(int i, int j, int height, int width)
{
    /* used for blocks at the boundaries */
    if ((i<height)&&(j<width)) return((i*width)+j);
    if (i>=height)
    {
        i = i - 8 + height - ((height>>3)<<3);
    }
    if (j>=width)
    {
        j = j - 8 + width - ((width>>3)<<3);
    }
    return((i*width)+j);
}

#define Linearize(i,j,w) (((i)*(w))+(j))


static int LastDC = 0;
static boolean LastDCwasPos = TRUE;

static void BlockStats(Pixel *image, int height, int width, boolean ErrImg,
                       int where_i, int where_j, Hist *H, FFLOAT *mssq,
                       FFLOAT *SigTab, boolean UseDCDPCM)
{
    int nc,ip,jp,ic,jc,lowip,lowjp,loc;
    FFLOAT sum, coef;
    int intcoef,temp;
    boolean positive;
    boolean btemp;


    lowip = where_i<<3;
    lowjp = where_j<<3;
    if (((lowip+8) > height) || ((lowjp+8) > width))
    {

        for (ip=lowip; ip<8+lowip; ip++)
        {
            for (jp=lowjp; jp<8+lowjp; jp++)
            {
                loc =  IJtoWhere(ip,jp,height,width);
                *mssq += ((FFLOAT) ((unsigned long)image[loc]*(unsigned long)image[loc]));
            }
        }

        nc = 0;
        for (ic=0; ic<8; ic++)
        {
            for (jc=0; jc<8; jc++,nc++)
            {
                sum = ((FFLOAT) 0.0);
                for (ip=lowip; ip<8+lowip; ip++)
                {
                    for (jp=lowjp; jp<8+lowjp; jp++)
                    {
                        sum += ((FFLOAT) image[IJtoWhere(ip,jp,height,width)])*
                               CosTable[((ip-lowip)<<3)+ic] *
                               CosTable[((jp-lowjp)<<3)+jc];
                    }
                }
                coef = sum * (0.25) * C(ic) * C(jc);
                SigTab[nc] += (coef*coef);
                if ((nc == 0) && (!ErrImg))
                {
                    coef = coef - ZERO_CORRECTION;
                    if (UseDCDPCM)
                    {
                        btemp = TRUE;
                        if (coef < 0.0) btemp = FALSE;
                        temp = Discretize(coef);
                        coef = coef - UnDiscretize(LastDC, LastDCwasPos);
                        LastDC = temp;
                        LastDCwasPos = btemp;
                    }
                }
                positive = TRUE;
                if (coef < 0.0) positive = FALSE;
                intcoef = Discretize(coef);
                HistIncrCount(H,nc,intcoef,positive);
            }
        }
    }
    else
    {


        for (ip=lowip; ip<8+lowip; ip++)
        {
            for (jp=lowjp; jp<8+lowjp; jp++)
            {
                loc =  Linearize(ip,jp,width);
                *mssq += ((FFLOAT) ((unsigned long)image[loc]*(unsigned long)image[loc]));
            }
        }

        nc = 0;
        for (ic=0; ic<8; ic++)
        {
            for (jc=0; jc<8; jc++,nc++)
            {
                sum = ((FFLOAT) 0.0);
                for (ip=lowip; ip<8+lowip; ip++)
                {
                    for (jp=lowjp; jp<8+lowjp; jp++)
                    {
                        sum += ((FFLOAT) image[Linearize(ip,jp,width)])*
                               CosTable[((ip-lowip)<<3)+ic] *
                               CosTable[((jp-lowjp)<<3)+jc];
                    }
                }
                coef = sum * (0.25) * C(ic) * C(jc);
                SigTab[nc] += (coef*coef);
                if ((nc == 0)  && (!ErrImg))
                {
                    coef = coef - ZERO_CORRECTION;
                    if (UseDCDPCM)
                    {
                        btemp = TRUE;
                        if (coef < 0.0) btemp = FALSE;
                        temp = Discretize(coef);
                        coef = coef - UnDiscretize(LastDC, LastDCwasPos);
                        LastDC = temp;
                        LastDCwasPos = btemp;
                    }
                }
                positive = TRUE;
                if (coef < 0.0) positive = FALSE;
                intcoef = Discretize(coef);
                HistIncrCount(H,nc,intcoef,positive);
            }
        }
    }
}



extern void SetHistogram(Hist *H, Image *Im, int cnum, FFLOAT *mssq,
                         FFLOAT *SigTab, boolean UseDCDPCM)
{
    int i,j,imax,jmax,h,w;

    LastDC = 0;
    LastDCwasPos = TRUE;

    h = Im->NumRows/Im->OutSamplingFactor[cnum][0];
    w = Im->NumCols/Im->OutSamplingFactor[cnum][1];
    imax = h/8;
    if ((imax << 3) != h) imax++;
    jmax = w/8;
    if ((jmax << 3) != w) jmax++;


    for (i=0; i<imax; i++)
    {
        for (j=0; j<jmax; j++)
        {
            BlockStats(Im->Im[cnum], h, w, Im->IsErrImage, i, j, H, mssq,
                       SigTab, UseDCDPCM);
        }
    }
}

