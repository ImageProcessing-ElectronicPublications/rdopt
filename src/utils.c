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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>

#include "rdopt.h"

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
    fprintf(stderr,"%s\n",s);
    exit(1);
}


extern void ReadIntTable(FILE *fp, int *tab)
{
    int n;
    for (n=0; n<64; n++)
    {
        if (fscanf(fp,"%d",&tab[n]) == EOF)
            FatalError("Could not read table off file");
    }
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

extern void WriteIntTable(FILE *fp, int *tab)
{
    int i,j;
    for (i=0; i<8; i++)
    {
        for (j=0; j<8; j++)
        {
            fprintf(fp,"%d ",tab[(i<<3)+j]);
        }
        fprintf(fp,"\n");
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


