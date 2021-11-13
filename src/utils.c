
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>

#include "rdopt.h"

static int ZigZagToN[64] =
{
    0, 1, 8, 16, 9, 2, 3, 10,
    17, 24, 32, 25, 18, 11, 4, 5,
    12, 19, 26, 33, 40, 48, 41, 34,
    27, 20, 13, 6, 7, 14, 21, 28,
    35, 42, 49, 56, 57, 50, 43, 36,
    29, 22, 15, 23, 30, 37, 44, 51,
    58, 59, 52, 45, 38, 31, 39, 46,
    53, 60, 61, 54, 47, 55, 62, 63,
};

extern int newline(char s[], FILE *fp)
{
    int i, c;
    i = 0;
    while ((c = getc(fp)) != EOF)
    {
        if (c == '\n')
        {
            s[i] = '\0';
            return (i+1);
        }
        s[i] = c;
        i++;
    }
    return EOF;
}

extern int newlinefd(char s[], int fd)
{
    int i;
    char c;
    i = 0;
    while (read(fd,&c,1) == 1)
    {
        if (c == '\n')
        {
            s[i] = '\0';
            return (i+1);
        }
        s[i] = c;
        i++;
    }
    return EOF;
}

extern void FatalError(char *s)
{
    fprintf(errfile,"%s\n",s);
    exit(1);
}


extern void ReadRealTable(FILE *fp, FFLOAT *tab)
{
    int n;
    double temp;
    for (n=0; n<64; n++)
    {
        if (fscanf(fp,"%lf",&temp) == EOF)
            FatalError("Could not read table off file");
        tab[n] = ((FFLOAT) temp);
    }
}

extern void WriteIntTable(FILE *fp, int *tab, char *prefix)
{
    int i,j;
    for (i=0; i<8; i++)
    {
        fprintf(fp,"%s",prefix);
        for (j=0; j<8; j++)
        {
            fprintf(fp,"%d ",tab[(i<<3)+j]);
        }
        fprintf(fp,"\n");
    }
}

extern void WriteThreshTable(FILE *fp, int *tab, char *prefix)
{
    int i,j;

    fprintf(fp,"# Table of thresholds\n");
    for (i=0; i<8; i++)
    {
        fprintf(fp,"%s",prefix);
        for (j=0; j<8; j++)
        {
            fprintf(fp,"%5.1f ", ((float) tab[(i<<3)+j])/2.0);
        }
        fprintf(fp,"\n");
    }
}

extern void WriteErrFBppDist(FILE *fp, FFLOAT *bpptab,
                             FFLOAT *errtab, char *prefix)
{
    int n, zn;

    fprintf(fp,"# Rate-Distortion distribution\n");
    fprintf(fp,"# ZigZagIndex CumulativeBpp   RMSE       PSNR\n");
    for (n=0; n<64; n++)
    {
        zn = ZigZagToN[n];
        fprintf(fp,"%s",prefix);
        fprintf(fp,"%2d             %7.6lf    %7.3lf %6.2lf\n",
                n, ((double) bpptab[zn]),
                sqrt((double) errtab[zn]),
                ((double) 10.0*log10(((double) 65025.0)/((double) errtab[zn]))));
    }
}
extern void WriteErrBppDist(FILE *fp, int bppscale, int *bpptab,
                            FFLOAT *errtab, char *prefix)
{
    int n, zn;

    fprintf(fp,"# Rate-Distortion distribution\n");
    fprintf(fp,"# ZigZagIndex CumulativeBpp   RMSE       PSNR\n");
    for (n=0; n<64; n++)
    {
        zn = ZigZagToN[n];
        fprintf(fp,"%s",prefix);
        fprintf(fp,"%2d             %7.6lf    %7.3lf %6.2lf\n",
                n, ((double) bpptab[zn])/((double) bppscale),
                sqrt((double) errtab[zn]),
                ((double) 10.0*log10(((double) 65025.0)/((double) errtab[zn]))));
    }
}

#define EQUALITY_DELTA ((FFLOAT) 1.0)

extern boolean NotEqualReal(FFLOAT x, FFLOAT y)
{
    FFLOAT diff;
    diff = x-y;
    if (diff < 0.0) diff = 0.0 - diff;
    if (diff <= EQUALITY_DELTA) return(FALSE);
    else return(TRUE);
}



extern void ReadIntTable(FILE *fp, int *tab)
{
    int j,k,l;
    char nextline[STRLENMAX], lilbuff[20];

    j = 0;

    while (newline(nextline,fp) != EOF)
    {
        if (!strncmp(nextline,"#END",4)) break;
        if (nextline[0] == '#') continue;
        k = 0;
        while (nextline[k] != '\0')
        {
            if ((nextline[k] >= '0') && (nextline[k] <= '9'))
            {
                l = 0;
                while ((nextline[k] >= '0') && (nextline[k] <= '9'))
                {
                    lilbuff[l] = nextline[k];
                    l++;
                    k++;
                }
                lilbuff[l] = '\0';
                tab[j] = atoi(lilbuff);
                j++;
                if (j == 64) break;
            }
            else k++;
        }
    }
    if (j != 64) FatalError("Could not read table of integers");
}

extern void SigToRemainingSig(FFLOAT *Sig, FFLOAT *RemainingSig)
{
    int i;

    RemainingSig[ZigZagToN[63]] = ((FFLOAT) 0.0);
    for (i=62; i>=0; i--)
        RemainingSig[ZigZagToN[i]]
            = Sig[ZigZagToN[i+1]]
              + RemainingSig[ZigZagToN[i+1]];
}

extern char * FileNameTail(char *s)
{
    int l;
    char *ans;

    l = strlen(s);
    ans = s + l;
    while ((l>0) && ((*(ans-1)) != '/'))
    {
        l--;
        ans--;
    }
    return(ans);
}


