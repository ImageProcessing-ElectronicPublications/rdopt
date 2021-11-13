
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
#include "qentry.h"

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


static int FindLambda(FFLOAT *llist, FFLOAT lambda, int num)
{
    int left, right, mid;

    if (lambda > llist[0]) return(0);
    if (lambda <= 0.0)
    {
        return(num-1);
    }

    left = 0;
    right = num;
    mid = num/2;
    while (mid > 0)
    {
        if ((lambda > llist[mid]) && (lambda <= llist[mid - 1]))
            return(mid);
        else if (llist[mid] > lambda)
        {
            left = mid;
            mid = (left + right)/2;
        }
        else
        {
            right = mid;
            mid = (left + right)/2;
        }
    }
    return(0);
}

static boolean EqualTabs(int (*tab1)[64], int (*tab2)[64], int numtabs)
{
    int u, n;
    for (u=0; u<numtabs; u++)
    {
        for (n=0; n<64; n++)
        {
            if (tab1[u][n] != tab2[u][n]) return(FALSE);
        }
    }
    return(TRUE);
}

static FFLOAT deltalambda = ((FFLOAT) 1.0e-8);

extern boolean
lagrFindRate(OptimJob *job, FFLOAT *bpp, FFLOAT *err, FFLOAT delta,
             FFLOAT lambdamin, FFLOAT lambdamax, FFLOAT *lambda,
             int (*anstab)[64])
{
    int unum, n;
    int tab1[MAXCOMPONENTS][64];
    int tab2[MAXCOMPONENTS][64];
    int tab3[MAXCOMPONENTS][64];

    int (*lefttab)[64], (*righttab)[64], (*midtab)[64], (*temptab)[64];
    FFLOAT leftbpp, rightbpp, lefterr, righterr, midbpp, miderr;
    FFLOAT lleft, lright, lmid;
    FFLOAT diff, absdiff, tempdiff, tempabsdiff;

    int where;


#define optstr job->opt_method.lagr

    lefttab = tab1;
    righttab = tab2;
    midtab = tab3;


    /* find first tables */

    lleft = lambdamax;
    leftbpp = ((FFLOAT) 0.0);
    lefterr = ((FFLOAT) 0.0);
    for (unum = 0; unum < job->NumTables; unum++)
    {
        for (n = 0; n < 64; n++)
        {
            where = FindLambda(optstr.Lambda[unum][n], lleft, optstr.IndexEntries[unum][n]);
            lefttab[unum][n] = optstr.SortedIndex[unum][n][where];
            leftbpp += optstr.Bpp[unum][n][lefttab[unum][n]];
            lefterr += optstr.Err[unum][n][lefttab[unum][n]];
        }
    }

    diff = leftbpp - (*bpp);
    if (diff < ((FFLOAT) 0.0)) absdiff = ((FFLOAT) 0.0) - diff;
    else absdiff = diff;

    lmid = lleft;
    for (unum = 0; unum < job->NumTables; unum++)
    {
        for (n = 0; n < 64; n++)
        {
            midtab[unum][n] = lefttab[unum][n];
        }
    }
    miderr = lefterr;
    midbpp = leftbpp;


    lright = lambdamin;
    rightbpp = ((FFLOAT) 0.0);
    righterr = ((FFLOAT) 0.0);
    for (unum = 0; unum < job->NumTables; unum++)
    {
        for (n = 0; n < 64; n++)
        {
            where = FindLambda(optstr.Lambda[unum][n], lright,
                               optstr.IndexEntries[unum][n]);
            righttab[unum][n] = optstr.SortedIndex[unum][n][where];
            rightbpp += optstr.Bpp[unum][n][righttab[unum][n]];
            righterr += optstr.Err[unum][n][righttab[unum][n]];
        }
    }

    tempdiff = rightbpp - (*bpp);
    if (tempdiff < ((FFLOAT) 0.0)) tempabsdiff = ((FFLOAT) 0.0) -tempdiff;
    else tempabsdiff = tempdiff;

    if (tempabsdiff < absdiff)
    {
        diff = tempdiff;
        absdiff = tempabsdiff;
        lmid = lright;
        miderr = righterr;
        midbpp = rightbpp;
        for (unum = 0; unum < job->NumTables; unum++)
        {
            for (n = 0; n < 64; n++)
            {
                midtab[unum][n] = righttab[unum][n];
            }
        }
    }

    while (absdiff > delta)
    {
        if (diff > ((FFLOAT) 0.0))
        {
            /* we are greater than target rate.. reduce lambda */
            if ((lmid - lright) < deltalambda) break;
            if (EqualTabs(midtab, righttab, job->NumTables)) break;
            lleft = lmid;
            temptab = lefttab;
            lefttab = midtab;
            midtab = temptab;
            lefterr = miderr;
            leftbpp = midbpp;
        }
        else
        {
            if ((lleft - lmid) < deltalambda) break;
            if (EqualTabs(midtab, lefttab, job->NumTables)) break;
            lright = lmid;
            temptab = righttab;
            righttab = midtab;
            midtab = temptab;
            righterr = miderr;
            rightbpp = midbpp;
        }


        lmid = (lleft + lright)/((FFLOAT) 2.0);
        midbpp = ((FFLOAT) 0.0);
        miderr = ((FFLOAT) 0.0);

        for (unum = 0; unum < job->NumTables; unum++)
        {
            for (n = 0; n < 64; n++)
            {
                where = FindLambda(optstr.Lambda[unum][n], lmid, optstr.IndexEntries[unum][n]);
                midtab[unum][n] = optstr.SortedIndex[unum][n][where];
                midbpp += optstr.Bpp[unum][n][midtab[unum][n]];
                miderr += optstr.Err[unum][n][midtab[unum][n]];
            }
        }


        diff = midbpp - (*bpp);
        if (diff < ((FFLOAT) 0.0)) absdiff = ((FFLOAT) 0.0) - diff;
        else absdiff = diff;
    }

    *bpp = midbpp;
    *err = miderr;
    *lambda = lmid;
    for (unum = 0; unum < job->NumTables; unum++)
    {
        for (n = 0; n < 64; n++)
        {
            anstab[unum][n] = midtab[unum][n];
        }
    }
    return(TRUE);
}


extern boolean
lagrFindErr(OptimJob *job, FFLOAT *bpp, FFLOAT *err, FFLOAT delta,
            FFLOAT lambdamin, FFLOAT lambdamax, FFLOAT *lambda,
            int (*anstab)[64])
{
    int unum, n;
    int tab1[MAXCOMPONENTS][64];
    int tab2[MAXCOMPONENTS][64];
    int tab3[MAXCOMPONENTS][64];

    int (*lefttab)[64], (*righttab)[64], (*midtab)[64], (*temptab)[64];
    FFLOAT leftbpp, rightbpp, lefterr, righterr, midbpp, miderr;
    FFLOAT lleft, lright, lmid;
    FFLOAT diff, absdiff, tempdiff, tempabsdiff;

    int where;


#define optstr job->opt_method.lagr

    lefttab = tab1;
    righttab = tab2;
    midtab = tab3;


    /* find first tables */

    lleft = lambdamax;
    leftbpp = ((FFLOAT) 0.0);
    lefterr = ((FFLOAT) 0.0);
    for (unum = 0; unum < job->NumTables; unum++)
    {
        for (n = 0; n < 64; n++)
        {
            where = FindLambda(optstr.Lambda[unum][n], lleft, optstr.IndexEntries[unum][n]);
            lefttab[unum][n] = optstr.SortedIndex[unum][n][where];
            leftbpp += optstr.Bpp[unum][n][lefttab[unum][n]];
            lefterr += optstr.Err[unum][n][lefttab[unum][n]];
        }
    }

    diff = lefterr - (*err);
    if (diff < ((FFLOAT) 0.0)) absdiff = ((FFLOAT) 0.0) - diff;
    else absdiff = diff;

    lmid = lleft;
    for (unum = 0; unum < job->NumTables; unum++)
    {
        for (n = 0; n < 64; n++)
        {
            midtab[unum][n] = lefttab[unum][n];
        }
    }
    miderr = lefterr;
    midbpp = leftbpp;


    lright = lambdamin;
    rightbpp = ((FFLOAT) 0.0);
    righterr = ((FFLOAT) 0.0);
    for (unum = 0; unum < job->NumTables; unum++)
    {
        for (n = 0; n < 64; n++)
        {
            where = FindLambda(optstr.Lambda[unum][n], lright, optstr.IndexEntries[unum][n]);
            righttab[unum][n] = optstr.SortedIndex[unum][n][where];
            rightbpp += optstr.Bpp[unum][n][righttab[unum][n]];
            righterr += optstr.Err[unum][n][righttab[unum][n]];
        }
    }

    tempdiff = righterr - (*err);
    if (tempdiff < ((FFLOAT) 0.0)) tempabsdiff = ((FFLOAT) 0.0)- tempdiff;
    else tempabsdiff = tempdiff;


    if (tempabsdiff < absdiff)
    {
        diff = tempdiff;
        absdiff = tempabsdiff;
        lmid = lright;
        miderr = righterr;
        midbpp = rightbpp;
        for (unum = 0; unum < job->NumTables; unum++)
        {
            for (n = 0; n < 64; n++)
            {
                midtab[unum][n] = righttab[unum][n];
            }
        }
    }

    while (absdiff > delta)
    {
        if (diff < ((FFLOAT) 0.0))
        {
            /* we are less than target error.. reduce lambda */
            if ((lmid - lright) < deltalambda) break;
            if (EqualTabs(midtab, righttab, job->NumTables)) break;
            lleft = lmid;
            temptab = lefttab;
            lefttab = midtab;
            midtab = temptab;
            lefterr = miderr;
            leftbpp = midbpp;
        }
        else
        {
            if ((lleft - lmid) < deltalambda) break;
            if (EqualTabs(midtab, lefttab, job->NumTables)) break;
            lright = lmid;
            temptab = righttab;
            righttab = midtab;
            midtab = temptab;
            righterr = miderr;
            rightbpp = midbpp;
        }


        lmid = (lleft + lright)/((FFLOAT) 2.0);
        midbpp = ((FFLOAT) 0.0);
        miderr = ((FFLOAT) 0.0);

        for (unum = 0; unum < job->NumTables; unum++)
        {
            for (n = 0; n < 64; n++)
            {
                where = FindLambda(optstr.Lambda[unum][n], lmid, optstr.IndexEntries[unum][n]);
                midtab[unum][n] = optstr.SortedIndex[unum][n][where];
                midbpp += optstr.Bpp[unum][n][midtab[unum][n]];
                miderr += optstr.Err[unum][n][midtab[unum][n]];
            }
        }


        diff = miderr - (*err);
        if (diff < ((FFLOAT) 0.0)) absdiff = ((FFLOAT) 0.0) - diff;
        else absdiff = diff;
    }

    *bpp = midbpp;
    *err = miderr;
    *lambda = lmid;
    for (unum = 0; unum < job->NumTables; unum++)
    {
        for (n = 0; n < 64; n++)
        {
            anstab[unum][n] = midtab[unum][n];
        }
    }
    return(TRUE);
}


extern void lagrRecoverQandT(OptimJob *job, int (*QandT)[64], int (*Q)[64],
                             int (*T)[64], FFLOAT (*Bpp)[64], FFLOAT (*Err)[64])
{
    int unum;
    int n;
    int tspanplus1 = job->ThreshSpan + 1;

    for (unum=0; unum < job->NumTables; unum++)
    {
        for (n=0; n<64; n++)
        {
            Q[unum][n] = QandT[unum][n]/tspanplus1;
            T[unum][n] = QandT[unum][n] - (Q[unum][n] * tspanplus1);

            Bpp[unum][n] = job->opt_method.lagr.Bpp[unum][n][QandT[unum][n]];
            Err[unum][n] = (job->opt_method.lagr.Err[unum][n][QandT[unum][n]]/
                            job->CompWeights[unum]);
        }
    }


    if (job->MapQ)
    {
        for (unum=0; unum < job->NumTables; unum++)
        {
            for (n=0; n<64; n++) Q[unum][n] = MapQentry(n,Q[unum][n]);
        }
    }

    for (unum=0; unum < job->NumTables; unum++)
    {
        for (n=0; n<64; n++) T[unum][n] += Q[unum][n];
    }

    for (unum=0; unum < job->NumTables; unum++)
    {
        /* convert Err, Bpp */
        for (n=1; n<64; n++)
        {
            Err[unum][ZigZagToN[n]] += Err[unum][ZigZagToN[n-1]];
            Bpp[unum][ZigZagToN[n]] += Bpp[unum][ZigZagToN[n-1]];
        }
        for (n=0; n<64; n++)
        {
            Err[unum][n] *= job->opt_method.lagr.ErrScale;
            Err[unum][n] += job->RemainingSignal[unum][n];
            Err[unum][n] /= job->NumPixels[unum];
        }
    }
}


static void ListCommands(void)
{
    fprintf(errfile,"Valid commands are:\n");
    fprintf(errfile,"\t [compress] size <target_size>\n");
    fprintf(errfile,"\t [compress] bpp <target_bpp>\n");
    fprintf(errfile,"\t [compress] psnr <target_psnr>\n");
    fprintf(errfile,"\t [compress] snr <target_snr>\n");
    fprintf(errfile,"\t [compress] rmse <target_rmse>\n");
    fprintf(errfile,"A file RDOPT.Q.<command>.<target> will be generated\n");

    fprintf(errfile,"  (unless a qfile command has already been issued)\n");
    fprintf(errfile,"If the 'compress' suffix is given, then a compressed\n");
    fprintf(errfile,"  file <imagename>.jpg will be generated\n");
    fprintf(errfile,"\t correct <use_bpp>\n");
    fprintf(errfile,"\t nocorrect\n");
    fprintf(errfile,"\t qfile fname (use fname as qtab file)\n");
    fprintf(errfile,"\t cfile fname (use fname as compressed file)\n");
    fprintf(errfile,"\t deltabpp value [def 0.001]\n");
    fprintf(errfile,"\t deltamse value [def 0.01]\n");
    fprintf(errfile,"\t deltalambda value [def 1e-8]\n");
    fprintf(errfile,"\t stats\n");
    fprintf(errfile,"\t nostats\n");
}

static void BriefListCommands(void)
{
    fprintf(errfile,"Valid commands: [compress] size/bpp/psnr/snr/rmse target\n");
    fprintf(errfile,"                correct bpp/nocorrect/help/qfile fname/stats/nostats/quit\n");
    fprintf(errfile,"                deltabpp val/deltamse val/deltalambda val\n");
}

extern void lagrCommand(OptimJob *job)
{
    int ResultSize;
    FFLOAT ResultSNR, ResultPSNR, TargetErr, TargetBpp;
    FFLOAT temp, actualBpp;
    double arg;
    int Q[MAXCOMPONENTS][64], T[MAXCOMPONENTS][64],
        QandT[MAXCOMPONENTS][64];
    FFLOAT ErrDist[MAXCOMPONENTS][64], BppDist[MAXCOMPONENTS][64];
    int PrintStats = 0;

    char tempqfname[STRLENMAX], cfname[STRLENMAX];
    char command[STRLENMAX], nextline[STRLENMAX], qfilename[STRLENMAX];
    boolean TargetIsBpp;
    int unum, numread;
    FILE *qfile;
    int infile;
    FFLOAT deltabpp, deltamse;
    FFLOAT lambdahigh,lambdalow,lambda;
    int i;
    boolean doingCorrection, doingCompression;
    boolean qfnameIsSet, cfnameIsSet;

    deltabpp = ((FFLOAT) 0.001);
    deltamse = ((FFLOAT) 0.01);

#define optstr job->opt_method.lagr

    if ((job->ReadCmdFromFile) && (strcmp(job->CmdFileName,"-")))
    {
        infile = open(job->CmdFileName,O_RDONLY,0);
        if (infile < 0) FatalError("Could not open command file");
    }
    else infile = 0;

    if ((!job->Silent) && (infile==0))
    {
        fprintf(errfile,"\n**** RD-OPT Command Interface ****\n\n");
        fflush(errfile);
    }

    qfnameIsSet = FALSE;
    cfnameIsSet = FALSE;

    while (TRUE)
    {

        doingCorrection = FALSE;
        doingCompression = FALSE;

        if ((!job->Silent) && (infile==0))
        {
            fprintf(errfile,"Command> ");
            fflush(errfile);
        }
        if (newlinefd(nextline,infile)==EOF) break;

        if (!strncmp(nextline,"qfile",5))
        {
            sscanf(nextline,"%s%s",command,qfilename);
            qfnameIsSet = TRUE;
            continue;
        }

        if (!strncmp(nextline,"cfile",5))
        {
            sscanf(nextline,"%s%s",command,cfname);
            cfnameIsSet = TRUE;
            continue;
        }

        if (!strncmp(nextline,"deltabpp",6))
        {
            sscanf(nextline,"%s%lf",command,&arg);
            deltabpp = ((FFLOAT) arg);
            continue;
        }

        if (!strncmp(nextline,"deltamse",6))
        {
            sscanf(nextline,"%s%lf",command,&arg);
            deltamse = ((FFLOAT) arg);
            continue;
        }

        if (!strncmp(nextline,"deltalambda",6))
        {
            sscanf(nextline,"%s%lf",command,&arg);
            deltalambda = ((FFLOAT) arg);
            continue;
        }

        if (!strncmp(nextline,"stats",5))
        {
            PrintStats = 1;
            continue;
        }

        if (!strncmp(nextline,"nostats",5))
        {
            PrintStats = 0;
            continue;
        }

        if (!strncmp(nextline,"nocorrect",5))
        {
            job->useCorrection = FALSE;
            continue;
        }

        if (!strncmp(nextline,"help",4))
        {
            ListCommands();
            continue;
        }

        if (!strncmp(nextline,"quit",4))
        {
            break;
        }

        if (!strncmp(nextline,"compress",4))
        {
            if (!strcmp(job->TheImage.ImFileName,""))
            {
                fprintf(errfile,"Cannot compress: image was read off stdin\n");
                continue;
            }
            doingCompression = TRUE;
            if (!cfnameIsSet) sprintf(cfname,"%s.jpg\0",
                                          FileNameTail(job->TheImage.ImFileName));
            strcpy(tempqfname,qfilename);
            sprintf(qfilename,"/tmp/rdopt.%s.compQ\0",
                    FileNameTail(job->TheImage.ImFileName));
            i = 0;
            while ((nextline[i] != ' ') && (nextline[i] != '\t')
                    && (nextline[i] != '\0')) nextline[i++] = ' ';
            /* fall thru */
        }

        numread = sscanf(nextline,"%s%lf",command,&arg);
        if (numread <= 0) goto nextLoop;



        if (!strncmp(command,"size",4))
        {
            if (numread < 2)
            {
                fprintf(errfile,"Must specify target size in bytes\n");
                goto nextLoop;
            }
            TargetIsBpp = TRUE;
            TargetBpp = ((FFLOAT) (arg*8.0));
            TargetBpp /= job->NumPixelsForBpp;
            if (job->useCorrection) TargetBpp -= job->addToTableBpp;
        }
        else if (!strncmp(command,"bpp",3))
        {
            if (numread < 2)
            {
                fprintf(errfile,"Must specify target bits-per-pixel\n");
                goto nextLoop;
            }
            TargetIsBpp = TRUE;
            TargetBpp = ((FFLOAT) arg);
            if (job->useCorrection) TargetBpp -= job->addToTableBpp;
        }
        else if ((!doingCompression) && (!strncmp(command,"correct",3)))
        {
            if (numread < 2)
            {
                fprintf(errfile,"Must specify correction bits-per-pixel\n");
                goto nextLoop;
            }
            TargetIsBpp = TRUE;
            TargetBpp = ((FFLOAT) arg);
            doingCorrection = TRUE;
            strcpy(tempqfname,qfilename);
            sprintf(qfilename,"/tmp/rdopt.%s.correctionQ\0",
                    FileNameTail(job->TheImage.ImFileName));
        }
        else if (!strncmp(command,"psnr",3))
        {
            if (numread < 2)
            {
                fprintf(errfile,"Must specify target PSNR in dB\n");
                goto nextLoop;
            }
            TargetIsBpp = FALSE;
            TargetErr = job->LogTotalPeakSq;
            TargetErr -= ((FFLOAT) arg/10.0);
            TargetErr = ((FFLOAT) pow(10.0,(double) TargetErr));
            TargetErr /= optstr.ErrScale;
        }
        else if (!strncmp(command,"snr",3))
        {
            if (numread < 2)
            {
                fprintf(errfile,"Must specify target SNR in dB\n");
                goto nextLoop;
            }
            TargetIsBpp = FALSE;
            TargetErr = job->LogTotalSignalSq;
            TargetErr -= ((FFLOAT) arg/10.0);
            TargetErr = ((FFLOAT) pow(10.0,(double) TargetErr));
            TargetErr /= optstr.ErrScale;

        }
        else if (!strncmp(command,"rmse",3))
        {
            if (numread < 2)
            {
                fprintf(errfile,"Must specify target RMSE\n");
                goto nextLoop;
            }
            TargetIsBpp = FALSE;
            TargetErr = ((FFLOAT) (arg*arg));
            TargetErr *= (job->TotalNumPixels/optstr.ErrScale);
        }
        else
        {
            BriefListCommands();
            goto nextLoop;
        }


        /******** extract Qtable(s) **********/
        if (TargetIsBpp)
        {
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

            lagrFindRate(job, &TargetBpp, &TargetErr, deltabpp,
                         lambdalow, lambdahigh, &lambda, QandT);
        }
        else
        {
            i = 0;
            while ((i < optstr.NumHooks) && (optstr.ErrHook[i] <= TargetErr))
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

            lagrFindErr(job, &TargetBpp, &TargetErr, deltamse,
                        lambdalow, lambdahigh, &lambda, QandT);

        }


        /******** open qfile ***********/
        if (!(qfnameIsSet || doingCorrection || doingCompression))
        {
            sprintf(qfilename,"RDOPT.Q.%s.%.3lf\0",command,arg);
        }

        if (!strcmp(qfilename,"-")) qfile = stdout;
        else
        {
            if ((qfile = fopen(qfilename,"w")) == NULL)
            {
                fprintf(errfile,"Could not open qtable file %s for writing\n",
                        qfilename);
                goto nextLoop;
            }
        }

        fprintf(qfile,"#%s\n",CQTAB);
        DumpJobChars(job, qfile);
        fprintf(qfile,"#Bits per pixel for each unit and the whole image\n");
        fprintf(qfile,"#    reported as num-bits/%lf, the denominator\n", ((double) job->NumPixelsForBpp));
        fprintf(qfile,"#    being the ");
        if (job->bppplane==-1)
            fprintf(qfile,"sum of num pixels in all planes\n#\n");
        else
            fprintf(qfile,"num pixels in plane number %d\n#\n",job->bppplane);
        fprintf(qfile,"# R + (lambda*D/%lf) minimized at l = %lf\n#\n",
                optstr.ErrScale, lambda);
        if (job->NumTables==1)
        {
            fprintf(qfile,"#Quantization table for color planes %d thru %d\n",
                    0,job->TheImage.NumComponents-1);
            lagrRecoverQandT(job,QandT,Q,T,BppDist,ErrDist);
            if (job->ThreshSpan > 0)
            {
                WriteIntTable(qfile, Q[0], "");
                WriteThreshTable(qfile, T[0], "#T ");
                if (PrintStats)
                    WriteErrFBppDist(qfile, BppDist[0], ErrDist[0], "#D ");
            }
            else
            {
                WriteIntTable(qfile, Q[0], "");
                if (PrintStats)
                    WriteErrFBppDist(qfile, BppDist[0], ErrDist[0], "#D ");
            }

            temp = TargetBpp * job->NumPixelsForBpp / 8.0;
            ResultSize = RoundOff(temp);

            if (job->useCorrection && !doingCorrection)
            {
                fprintf(qfile,"#RdoptBpp = %lf bits per pixel\n",((double)TargetBpp));
                fprintf(qfile,"#RdoptSize = %d bytes\n", ResultSize);
                fprintf(qfile,"#CorrectionToBpp = %lf bits per pixel\n",
                        (double) job->addToTableBpp);
                TargetBpp += job->addToTableBpp;
                temp = TargetBpp * job->NumPixelsForBpp / 8.0;
                ResultSize = RoundOff(temp);
            }

            ResultSNR = (job->LogTotalSignalSq -
                         ((FFLOAT) log10( (double) (
                                              (TargetErr*optstr.ErrScale)))))*10.0;
            ResultPSNR = (job->LogTotalPeakSq -
                          ((FFLOAT) log10( (double) (
                                               (TargetErr*optstr.ErrScale)))))*10.0;

            fprintf(qfile,"#Bpp = %lf bits per pixel\n", ((double) TargetBpp));
            fprintf(qfile,"#Size = %d bytes\n", ResultSize);
            fprintf(qfile,"#SNR = %lf dB\n", ((double) ResultSNR));
            fprintf(qfile,"#PSNR = %lf dB\n", ((double) ResultPSNR));
            TargetErr *= (optstr.ErrScale/job->TotalNumPixels);
            fprintf(qfile,"#RMSE = %lf\n", sqrt((double) TargetErr));

        }
        else
        {
            lagrRecoverQandT(job,QandT,Q,T,BppDist,ErrDist);
            for (unum=0; unum < job->NumTables; unum++)
            {
                fprintf(qfile,"#Quantization table for color planes %d thru %d\n",unum,((unum==(job->NumTables-1))?(job->TheImage.NumComponents-1):unum));

                if (job->ThreshSpan > 0)
                {
                    WriteIntTable(qfile, Q[unum], "");
                    WriteThreshTable(qfile, T[unum], "#T ");
                    if (PrintStats)
                        WriteErrFBppDist(qfile, BppDist[unum],ErrDist[unum], "#D ");
                }
                else
                {
                    WriteIntTable(qfile, Q[unum], "");
                    if (PrintStats)
                        WriteErrFBppDist(qfile, BppDist[unum], ErrDist[unum],"#D ");
                }

                temp = BppDist[unum][63] * job->NumPixelsForBpp / 8.0;
                ResultSize = RoundOff(temp);
                ResultSNR = (job->LogSignalSq[unum] -
                             ((FFLOAT) log10( (double) (
                                                  (ErrDist[unum][63]*job->NumPixels[unum])))))*10.0;
                ResultPSNR = (job->LogPeakSq[unum] -
                              ((FFLOAT) log10( (double) (
                                                   (ErrDist[unum][63]*job->NumPixels[unum])))))*10.0;

                fprintf(qfile,"#Bpp = %lf bits per pixel\n",
                        ((double) BppDist[unum][63]));
                fprintf(qfile,"#Size = %d bytes\n", ResultSize);
                fprintf(qfile,"#SNR = %lf dB\n", ((double) ResultSNR));
                fprintf(qfile,"#PSNR = %lf dB\n", ((double) ResultPSNR));
                fprintf(qfile,"#RMSE = %lf\n", sqrt(ErrDist[unum][63]));
                fprintf(qfile,"#\n");

            }
            fprintf(qfile,"#For the entire image:\n");
            temp = TargetBpp * job->NumPixelsForBpp/ 8.0;
            ResultSize = RoundOff(temp);
            ResultSNR = (job->LogTotalSignalSq -
                         ((FFLOAT) log10( (double) (
                                              (TargetErr*optstr.ErrScale)))))*10.0;
            ResultPSNR = (job->LogTotalPeakSq -
                          ((FFLOAT) log10( (double) (
                                               (TargetErr*optstr.ErrScale)))))*10.0;

            if (job->useCorrection && !doingCorrection)
            {
                fprintf(qfile,"#RdoptBpp = %lf bits per pixel\n",((double)TargetBpp));
                fprintf(qfile,"#RdoptSize = %d bytes\n", ResultSize);
                fprintf(qfile,"#CorrectionToBpp = %lf bits per pixel\n",
                        (double) job->addToTableBpp);
                TargetBpp += job->addToTableBpp;
                temp = TargetBpp * job->NumPixelsForBpp / 8.0;
                ResultSize = RoundOff(temp);
            }


            fprintf(qfile,"#Bpp = %lf bits per pixel\n", ((double) TargetBpp));
            fprintf(qfile,"#Size = %d bytes\n", ResultSize);
            fprintf(qfile,"#SNR = %lf dB\n", ((double) ResultSNR));
            fprintf(qfile,"#PSNR = %lf dB\n", ((double) ResultPSNR));
            fprintf(qfile,"#RMSE = %lf\n",
                    sqrt((double) (TargetErr*optstr.ErrScale/job->TotalNumPixels)));
        }

        fprintf(qfile,"#END\n");

        if (strcmp(qfilename,"-")) fclose(qfile);

        if (doingCorrection)
        {
            if (GetActualBpp(job, qfilename, &actualBpp))
            {
                job->useCorrection = TRUE;
                job->correctionBpp = TargetBpp;
                job->addToTableBpp = actualBpp - TargetBpp;
                if ((!job->Silent) && (infile==0))
                {
                    fprintf(errfile,"Correction = %lf - %lf = %lf bpp\n",
                            actualBpp, TargetBpp, job->addToTableBpp);
                    fflush(errfile);
                }
            }
            else
            {
                if ((!job->Silent) && (infile==0))
                {
                    fprintf(errfile,"Error executing cjpeg\n");
                    fflush(errfile);
                }
            }
            remove(qfilename);
        }
        else if (doingCompression)
        {
            if (Compress(job, qfilename, cfname, &ResultSize))
            {
                if ((!job->Silent) && (infile==0))
                {
                    fprintf(errfile,"Size = %d  bytes (%lf bpp)\n",
                            ResultSize,
                            ((double) ResultSize*8)/((double) job->NumPixelsForBpp));
                    fflush(errfile);
                }
            }
            else
            {
                if ((!job->Silent) && (infile==0))
                {
                    fprintf(errfile,"Error executing cjpeg\n");
                    fflush(errfile);
                }
            }
            remove(qfilename);
        }


nextLoop:
        if (doingCorrection || doingCompression)
            strcpy(qfilename,tempqfname);




    }
    if ((job->ReadCmdFromFile) && (strcmp(job->CmdFileName,"-")))
        close(infile);
}


