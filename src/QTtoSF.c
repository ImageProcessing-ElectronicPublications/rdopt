
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

#define MAXCOMPS 5

int getaline(FILE *fp, char *buff, int lim)
{
    int n = 0,c;
    while (((c = getc(fp)) != '\n') && (c != EOF))
    {
        buff[n++] = c;
        if (n >= (lim-1))
        {
            buff[lim-1] = '\0';
            return(1);
        }
    }
    if ((n==0) && (c==EOF)) return(0);
    else
    {
        buff[n] = '\0';
        return(1);
    }
}

int
main (int argc, char *argv[])
{
    FILE *qtfile, *sffile;
    int numACscans;
    int nplanes;
    double bpp[MAXCOMPS][64], psnr[MAXCOMPS][64];
    double bppspan[MAXCOMPS], bppincr[MAXCOMPS], currbpp;
    int lastcoeff[MAXCOMPS];
    char nextline[200], *ptr;
    int n, p, s, tc;
    int ntemp;
    double trmse;

    if (argc < 3)
    {
        fprintf(stderr,"Usage: QTtoSF <numACscans> qtabfile [sffile]\n");
        exit(1);
    }

    numACscans = atoi(argv[1]);
    if ((numACscans < 1) || (numACscans > 63))
    {
        fprintf(stderr,"numACscans must be in range 1-63\n");
        exit(1);
    }

    if ((qtfile = fopen(argv[2],"r")) == NULL)
    {
        fprintf(stderr,"Could not open %s\n",argv[2]);
        exit(1);
    }

    if (argc > 3)
    {
        if ((sffile = fopen(argv[3],"w")) == NULL)
        {
            fprintf(stderr,"Could not open %s\n",argv[3]);
            exit(1);
        }
    }
    else sffile = stdout;

    nplanes = 0;
    n = 0;
    while (getaline(qtfile, nextline, 200) == 1)
    {
        if (strncmp(nextline,"#D ",3)) continue;
        ptr = &nextline[3];
        sscanf(ptr,"%d%lf%lf%lf", &ntemp,&bpp[nplanes][n],
               &trmse, &psnr[nplanes][n]);
        if (ntemp != n)
        {
            fprintf(stderr,"Bad format in %s\n",argv[2]);
            exit(1);
        }
        n++;
        if (n==64)
        {
            n=0;
            nplanes++;
            if (nplanes >= MAXCOMPS)
            {
                fprintf(stderr,"Too many components\n");
                exit(1);
            }
        }
    }

    fclose(qtfile);

    for (p=0; p<(nplanes-1); p++) fprintf(sffile,"%d,",p);
    fprintf(sffile,"%d: 0 0 0 0;\n",nplanes-1);

    for (p=0; p<nplanes; p++)
    {
        bppspan[p] = bpp[p][63] - bpp[p][0];
        bppincr[p] = bppspan[p]/((double) numACscans);
        lastcoeff[p] = 0;
    }

    for (s=1; s<=numACscans; s++)
    {
        for (p=0; p<nplanes; p++)
        {
            tc = lastcoeff[p] + 1;
            if (s==numACscans) lastcoeff[p] = 63;
            else
            {
                lastcoeff[p]++;
                currbpp = bpp[p][0] + (s*bppincr[p]);
                while ((lastcoeff[p] < (64 - numACscans + s)) &&
                        (bpp[p][lastcoeff[p]] < currbpp))
                {
                    lastcoeff[p]++;
                }
            }
            fprintf(sffile,"%d: %d %d 0 0;\n",p,tc,lastcoeff[p]);

        }
    }

    if (sffile != stdout) fclose(sffile);
    return 0;
}

