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


extern void HistIncrCount(Hist *H, int n, int val,
                          boolean waspositive)
{
    int absval, newsize, *temp;
    if (waspositive)
    {
        if (H[n].PlusSize <= val)
        {
            newsize = H[n].PlusSize << 1;
            while (newsize <= val) newsize <<= 1;
            if ((temp = (int *) calloc(1,newsize*sizeof(int))) == NULL)
            {
                FatalError("HistIncrCount: out of memory");
            }
            memcpy((char *) temp, (char *) H[n].PlusCount,
                   sizeof(int) * H[n].PlusSize);
            free(H[n].PlusCount);
            H[n].PlusCount = temp;
            H[n].PlusSize = newsize;
        }
        H[n].PlusCount[val]++;
    }
    else
    {
        absval = 0-val;
        if (H[n].MinusSize <= absval)
        {
            newsize = H[n].MinusSize << 1;
            while (newsize <= absval) newsize <<= 1;
            if ((temp = (int *) calloc(1,newsize*sizeof(int))) == NULL)
            {
                FatalError("HistIncrCount: out of memory");
            }
            memcpy((char *) temp, (char *) H[n].MinusCount,
                   sizeof(int) * H[n].MinusSize);
            free(H[n].MinusCount);
            H[n].MinusCount = temp;
            H[n].MinusSize = newsize;
        }
        H[n].MinusCount[absval]++;
    }
}



static void HistSetCount(Hist *H, int n, int val,
                         boolean waspositive, int count)
{
    int absval, newsize, *temp;
    if (waspositive)
    {
        if (H[n].PlusSize <= val)
        {
            newsize = H[n].PlusSize << 1;
            while (newsize <= val) newsize <<= 1;
            if ((temp = (int *) calloc(1,newsize*sizeof(int))) == NULL)
            {
                FatalError("HistSetCount: out of memory");
            }
            memcpy((char *) temp, (char *) H[n].PlusCount,
                   sizeof(int) * H[n].PlusSize);
            free(H[n].PlusCount);
            H[n].PlusCount = temp;
            H[n].PlusSize = newsize;
        }
        H[n].PlusCount[val] = count;
    }
    else
    {
        absval = 0-val;
        if (H[n].MinusSize <= absval)
        {
            newsize = H[n].MinusSize << 1;
            while (newsize <= absval) newsize <<= 1;
            if ((temp = (int *) calloc(1,newsize*sizeof(int))) == NULL)
            {
                FatalError("HistSetCount: out of memory");
            }
            memcpy((char *) temp, (char *) H[n].MinusCount,
                   sizeof(int) * H[n].MinusSize);
            free(H[n].MinusCount);
            H[n].MinusCount = temp;
            H[n].MinusSize = newsize;
        }
        H[n].MinusCount[absval] = count;
    }
}


extern void InitHistogram(Hist *H)
{
    int n;

    for (n=0; n<64; n++)
    {
        H[n].PlusSize = 16;
        if ((H[n].PlusCount = (int *) calloc(1,16*sizeof(int))) == NULL)
        {
            FatalError("InitHistogram: out of memory");
        }
        H[n].MinusSize = 16;
        if ((H[n].MinusCount = (int *) calloc(1,16*sizeof(int))) == NULL)
        {
            FatalError("InitHistogram: out of memory");
        }
    }
}


/***********************************************************
Format of histogram files:
#RDOPT.H<version>
#<Description-String>
#Mean-Signal-Squared
#Coefficient number 0
#Negative values
<value> <count>
<value> <count>
  ..      ..
  ..      ..
#Positive values
<value> <count>
<value> <count>
  ..      ..
  ..      ..
#Coefficient number 1
  ..      ..
  ..      ..
  ..      ..
  ..      ..
#Coefficient number 63
  ..      ..
  ..      ..
#END
**************************************************************/


extern void ReadHistogram(FILE *hfile, Hist *H, FFLOAT *mssq)
{
    char nextline[STRLENMAX];
    int n, v, c;
    double dtemp;

    newline(nextline,hfile);
    if (strncmp(nextline,"#RDOPT.H",8))
    {
        FatalError("ReadHistogram: unknown format");
    }
    newline(nextline,hfile); /* discard description string */

    newline(nextline,hfile);
    sscanf(nextline,"#%lf",&dtemp);
    *mssq = ((FFLOAT) dtemp);


    /* read the "#Coefficient number 0" line */
    newline(nextline,hfile);

    for (n=0; n<64; n++)
    {
        newline(nextline,hfile); /* #Negative values */

        newline(nextline,hfile);
        while (nextline[0] != '#')
        {
            sscanf(nextline,"%d%d",&v,&c);
            HistSetCount(H,n,v,FALSE,c);
            newline(nextline,hfile);
        }

        /* "#Positive values" has been read */
        newline(nextline,hfile);
        while (nextline[0] != '#')
        {
            sscanf(nextline,"%d%d",&v,&c);
            HistSetCount(H,n,v,TRUE,c);
            newline(nextline,hfile);
        }
    }
    if (strncmp(nextline,"#END",4)) FatalError("ReadHistogram: Unknown format");

}



extern void DumpHistogram(FILE *hfile, Hist *H, char *Description, FFLOAT mssq)
{
    int n,v;

    fprintf(hfile,"#%s\n",CQHIST);
    fprintf(hfile,"#%s\n",Description);
    fprintf(hfile,"#%lf\n",(double) mssq);

    for (n=0; n<64; n++)
    {

        fprintf(hfile,"#Coefficient number %d\n",n);

        fprintf(hfile,"#Negative values\n");
        for (v=1-H[n].MinusSize; v<=0; v++)
        {
            if (H[n].MinusCount[0-v])
            {
                fprintf(hfile,"%d %d\n",v,H[n].MinusCount[0-v]);
            }
        }

        fprintf(hfile,"#Positive values\n");
        for (v=H[n].PlusSize-1; v>=0; v--)
        {
            if (H[n].PlusCount[v])
            {
                fprintf(hfile,"%d %d\n",v,H[n].PlusCount[v]);
            }
        }

    }
    fprintf(hfile,"#END\n");
}

extern void FreeHistogram(Hist *H)
{
    int n;
    for (n=0; n<64; n++)
    {
        free(H[n].PlusCount);
        free(H[n].MinusCount);
    }
}

extern FFLOAT GetTotalBlocks(Hist *H)
{
    int n;
    FFLOAT count = 0.0;

    for (n=0; n<H[63].PlusSize; n++)
        count += ((FFLOAT) H[63].PlusCount[n]);
    for (n=0; n<H[63].MinusSize; n++)
        count += ((FFLOAT) H[63].MinusCount[n]);

    return(count);
}
