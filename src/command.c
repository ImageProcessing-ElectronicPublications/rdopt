
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



extern int ScanRowForErr(FFLOAT *row, FFLOAT err, int least, int most)
{
    int i;
    for (i=least; i<=most; i++)
    {
        if (row[i] < err) return(i);
    }
    return(most);
}

extern int ScanRowForBpp(FFLOAT *row, int bpp, int least, int most)
{
    int i;
    if (bpp < least) return(least);
    if (bpp > most) return(most);
    i = bpp;
    while (i>=least)
    {
        if (row[i] < INFINITY) return(i);
        i--;
    }
    /* will not come here */
    return(-1);
}

extern void GetStartingPoints(int *ToCombine[], int *StartingPt,
                              int Loc, int NumUnits)
{
    int i, NextLoc;
    NextLoc = Loc;
    for (i=NumUnits-1; i>=1; i--)
    {
        StartingPt[i] = ToCombine[i][NextLoc];
        NextLoc -= StartingPt[i];
    }
    StartingPt[0] = NextLoc;
}


extern void RecoverQ(int Loc, int unum, OptimJob *job, int *Q, int *Bpp,
                     FFLOAT *Err)
{
    /* ThreshSpan was 0 */
    int n, NextLoc;

    NextLoc = Loc;

    for (n=0; n<64; n++)
    {
        Q[n] = (int) job->opt_method.dp.QChoice[unum][n][NextLoc];

        Bpp[n] = job->opt_method.dp.Bpp[unum][n][Q[n]];

        if (job->opt_method.dp.ErrEncodesRow[unum][n])
        {
            Err[n] = job->opt_method.dp.Err[unum][n][Bpp[n] -
                     job->opt_method.dp.BppOffsetInEncoding[unum][n]];
        }
        else
        {
            Err[n] = job->opt_method.dp.Err[unum][n][Q[n]];
        }
        Err[n] /= job->CompWeights[unum];


        NextLoc -= Bpp[n];
    }

    /* convert Err and Bpp to cumulative zigzag values */
    for (n=1; n<64; n++)
    {
        Err[ZigZagToN[n]] += Err[ZigZagToN[n-1]];
        Bpp[ZigZagToN[n]] += Bpp[ZigZagToN[n-1]];
    }
    for (n=0; n<64; n++)
    {
        Err[n] += job->RemainingSignal[unum][n];
        Err[n] /= job->NumPixels[unum];
    }



    if (job->MapQ)
    {
        for (n=0; n<64; n++) Q[n] = MapQentry(n,Q[n]);
    }

}

extern void RecoverQandT(int Loc, int unum, OptimJob *job, int *Q, int *T, int *Bpp, FFLOAT *Err)
{
    int n, NextLoc, indx;
    NextLoc = Loc;
    for (n=0; n<64; n++)
    {
        Q[n] = (int) job->opt_method.dp.QChoice[unum][n][NextLoc];
        T[n] = (int) job->opt_method.dp.TChoice[unum][n][NextLoc];
        indx = Q[n]*(job->ThreshSpan+1) + T[n];

        Bpp[n] = job->opt_method.dp.Bpp[unum][n][indx];

        NextLoc -= Bpp[n];

        if (job->opt_method.dp.ErrEncodesRow[unum][n])
        {
            Err[n] = job->opt_method.dp.Err[unum][n][Bpp[n] -
                     job->opt_method.dp.BppOffsetInEncoding[unum][n]];
        }
        else
        {
            Err[n] = job->opt_method.dp.Err[unum][n][indx];
        }
        Err[n] /= job->CompWeights[unum];
    }

    if (job->MapQ)
    {
        for (n=0; n<64; n++) Q[n] = MapQentry(n,Q[n]);
    }

    for (n=0; n<64; n++) T[n] += Q[n];

    /* convert Err */
    for (n=1; n<64; n++)
    {
        Err[ZigZagToN[n]] += Err[ZigZagToN[n-1]];
        Bpp[ZigZagToN[n]] += Bpp[ZigZagToN[n-1]];
    }
    for (n=0; n<64; n++)
    {
        Err[n] += job->RemainingSignal[unum][n];
        Err[n] /= job->NumPixels[unum];
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
    fprintf(errfile,"\t stats\n");
    fprintf(errfile,"\t nostats\n");
}

static void BriefListCommands(void)
{
    fprintf(errfile,"Valid commands: [compress] size/bpp/psnr/snr/rmse target\n");
    fprintf(errfile,"                correct bpp/nocorrect/help/qfile fname/stats/nostats/quit\n");
}

extern void Command(OptimJob *job)
{
    int IntTargetBpp, ResultSize;
    FFLOAT ResultBpp, ResultSNR, ResultPSNR, TargetErr, RealTargetBpp;
    FFLOAT ResultErr, temp, actualBpp;
    double arg;
    int Q[64], T[64], BppDist[64];
    FFLOAT ErrDist[64];
    int PrintStats = 0;

    char command[STRLENMAX], nextline[STRLENMAX], qfilename[STRLENMAX];
    char cfname[STRLENMAX];
    char tempqfname[STRLENMAX];
    boolean TargetIsBpp;
    int unum, StartingPt[MAXCOMPONENTS], numread, Loc;
    FILE *qfile;
    int infile;
    boolean doingCorrection, doingCompression;
    boolean qfnameIsSet, cfnameIsSet;
    int i;

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
            RealTargetBpp = ((FFLOAT) (arg*8.0));
            RealTargetBpp /= job->NumPixelsForBpp;
            if (job->useCorrection) RealTargetBpp -= job->addToTableBpp;

            RealTargetBpp *= ((FFLOAT) job->opt_method.dp.BppScale);
            IntTargetBpp = RoundOff(RealTargetBpp);

        }
        else if (!strncmp(command,"bpp",3))
        {
            if (numread < 2)
            {
                fprintf(errfile,"Must specify target bits-per-pixel\n");
                goto nextLoop;
            }
            TargetIsBpp = TRUE;
            RealTargetBpp = ((FFLOAT) arg);
            if (job->useCorrection) RealTargetBpp -= job->addToTableBpp;
            RealTargetBpp *= ((FFLOAT) job->opt_method.dp.BppScale);
            IntTargetBpp = RoundOff(RealTargetBpp);
        }
        else if ((!doingCompression) && (!strncmp(command,"correct",3)))
        {
            if (numread < 2)
            {
                fprintf(errfile,"Must specify correction bits-per-pixel\n");
                goto nextLoop;
            }
            if (!strcmp(job->TheImage.ImFileName,""))
            {
                fprintf(errfile,"Cannot compress: image was read off stdin\n");
                goto nextLoop;
            }
            TargetIsBpp = TRUE;
            doingCorrection = TRUE;
            RealTargetBpp = ((FFLOAT) arg);
            RealTargetBpp *= ((FFLOAT) job->opt_method.dp.BppScale);
            IntTargetBpp = RoundOff(RealTargetBpp);
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

        }
        else if (!strncmp(command,"rmse",3))
        {
            if (numread < 2)
            {
                fprintf(errfile,"Must specify target RMSE\n");
                goto nextLoop;
            }
            TargetIsBpp = FALSE;
            TargetErr = ((FFLOAT) log10((double) arg))*2.0;
            TargetErr += ((FFLOAT) log10((double) job->TotalNumPixels));
        }
        else
        {
            BriefListCommands();
            goto nextLoop;
        }


        /******** extract Qtable(s) **********/
        if (job->NumTables==1)
        {
            if (TargetIsBpp)
            {
                Loc = ScanRowForBpp(job->opt_method.dp.LastRow[0],IntTargetBpp,
                                    job->opt_method.dp.LeastLast[0], job->opt_method.dp.MostLast[0]);
            }
            else
            {
                Loc = ScanRowForErr(job->opt_method.dp.LastRow[0],TargetErr,
                                    job->opt_method.dp.LeastLast[0], job->opt_method.dp.MostLast[0]);
            }
            if (Loc < 0)
            {
                /* should not happen */
                fprintf(errfile,"Set target cannot be achieved\n");
                goto nextLoop;
            }
        }
        else
        {
            if (TargetIsBpp)
            {
                Loc = ScanRowForBpp(job->opt_method.dp.CombinedLastRow,IntTargetBpp,
                                    job->opt_method.dp.CombinedLeast, job->opt_method.dp.CombinedMost);
            }
            else
            {
                Loc = ScanRowForErr(job->opt_method.dp.CombinedLastRow,TargetErr,
                                    job->opt_method.dp.CombinedLeast, job->opt_method.dp.CombinedMost);
            }
            if (Loc < 0)
            {
                /* should not happen */
                fprintf(errfile,"Set target cannot be achieved\n");
                goto nextLoop;
            }
            GetStartingPoints(job->opt_method.dp.ToCombine, StartingPt,
                              Loc, job->NumTables);
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

        if (job->NumTables==1)
        {
            fprintf(qfile,"#Quantization table for color planes %d thru %d\n",0,job->TheImage.NumComponents-1);

            if (job->ThreshSpan > 0)
            {
                RecoverQandT(Loc, 0, job, Q, T, BppDist, ErrDist);
                WriteIntTable(qfile, Q, "");
                WriteThreshTable(qfile, T, "#T ");
                if (PrintStats)
                    WriteErrBppDist(qfile, job->opt_method.dp.BppScale,
                                    BppDist, ErrDist, "#D ");
            }
            else
            {
                RecoverQ(Loc, 0, job, Q, BppDist, ErrDist);
                WriteIntTable(qfile, Q, "");
                if (PrintStats)
                    WriteErrBppDist(qfile, job->opt_method.dp.BppScale,
                                    BppDist, ErrDist, "#D ");
            }

            ResultErr = job->opt_method.dp.LastRow[0][Loc];
            ResultBpp = ((FFLOAT) Loc)/((FFLOAT) job->opt_method.dp.BppScale);
            temp = ResultBpp * job->NumPixelsForBpp / 8.0;
            ResultSize = RoundOff(temp);

            if (job->useCorrection && !doingCorrection)
            {
                fprintf(qfile,"#RdoptBpp = %lf bits per pixel\n",((double)ResultBpp));
                fprintf(qfile,"#RdoptSize = %d bytes\n", ResultSize);
                fprintf(qfile,"#CorrectionToBpp = %lf bits per pixel\n",
                        (double) job->addToTableBpp);
                ResultBpp += job->addToTableBpp;
                temp = ResultBpp * job->NumPixelsForBpp / 8.0;
                ResultSize = RoundOff(temp);
            }

            ResultSNR = (job->LogTotalSignalSq - ResultErr)*10.0;
            ResultPSNR = (job->LogTotalPeakSq - ResultErr)*10.0;

            fprintf(qfile,"#Bpp = %lf bits per pixel\n", ((double) ResultBpp));
            fprintf(qfile,"#Size = %d bytes\n", ResultSize);
            fprintf(qfile,"#SNR = %lf dB\n", ((double) ResultSNR));
            fprintf(qfile,"#PSNR = %lf dB\n", ((double) ResultPSNR));
            ResultErr -= ((FFLOAT) log10((double) job->TotalNumPixels));
            fprintf(qfile,"#RMSE = %lf\n",
                    sqrt(pow(10.0,((double) ResultErr))));

        }
        else
        {
            for (unum=0; unum < job->NumTables; unum++)
            {
                fprintf(qfile,"#Quantization table for color planes %d thru %d\n",unum,((unum==(job->NumTables-1))?(job->TheImage.NumComponents-1):unum));

                if (job->ThreshSpan > 0)
                {
                    RecoverQandT(StartingPt[unum], unum, job, Q, T, BppDist, ErrDist);
                    WriteIntTable(qfile, Q, "");
                    WriteThreshTable(qfile, T, "#T ");
                    if (PrintStats)
                        WriteErrBppDist(qfile, job->opt_method.dp.BppScale,
                                        BppDist, ErrDist, "#D ");
                }
                else
                {
                    RecoverQ(StartingPt[unum], unum, job, Q, BppDist, ErrDist);
                    WriteIntTable(qfile, Q, "");
                    if (PrintStats)
                        WriteErrBppDist(qfile, job->opt_method.dp.BppScale,
                                        BppDist, ErrDist, "#D ");
                }
                ResultErr = job->opt_method.dp.LastRow[unum][StartingPt[unum]];
                ResultBpp = ((FFLOAT) StartingPt[unum])/((FFLOAT) job->opt_method.dp.BppScale);
                temp = ResultBpp * job->NumPixelsForBpp / 8.0;
                ResultSize = RoundOff(temp);
                ResultSNR = (job->LogSignalSq[unum] - ResultErr)*10.0;
                ResultPSNR = (job->LogPeakSq[unum] - ResultErr)*10.0;

                fprintf(qfile,"#Bpp = %lf bits per pixel\n", ((double) ResultBpp));
                fprintf(qfile,"#Size = %d bytes\n", ResultSize);
                fprintf(qfile,"#SNR = %lf dB\n", ((double) ResultSNR));
                fprintf(qfile,"#PSNR = %lf dB\n", ((double) ResultPSNR));
                ResultErr -= ((FFLOAT) log10((double) job->NumPixels[unum]));
                fprintf(qfile,"#RMSE = %lf\n", sqrt(pow(10.0,((double) ResultErr))));
                fprintf(qfile,"#\n");

            }
            fprintf(qfile,"#For the entire image:\n");
            ResultErr = job->opt_method.dp.CombinedLastRow[Loc];
            ResultBpp = ((FFLOAT) Loc)/((FFLOAT) job->opt_method.dp.BppScale);
            temp = ResultBpp * job->NumPixelsForBpp/ 8.0;
            ResultSize = RoundOff(temp);
            if (job->useCorrection && !doingCorrection)
            {
                fprintf(qfile,"#RdoptBpp = %lf bits per pixel\n", ((double) ResultBpp));
                fprintf(qfile,"#RdoptSize = %d bytes\n", ResultSize);
                fprintf(qfile,"#CorrectionToBpp = %lf bits per pixel\n",
                        (double) job->addToTableBpp);
                ResultBpp += job->addToTableBpp;
                temp = ResultBpp * job->NumPixelsForBpp / 8.0;
                ResultSize = RoundOff(temp);
            }
            ResultSNR = (job->LogTotalSignalSq - ResultErr)*10.0;
            ResultPSNR = (job->LogTotalPeakSq - ResultErr)*10.0;

            fprintf(qfile,"#Bpp = %lf bits per pixel\n", ((double) ResultBpp));
            fprintf(qfile,"#Size = %d bytes\n", ResultSize);
            fprintf(qfile,"#SNR = %lf dB\n", ((double) ResultSNR));
            fprintf(qfile,"#PSNR = %lf dB\n", ((double) ResultPSNR));
            ResultErr -= ((FFLOAT) log10((double) job->TotalNumPixels));
            fprintf(qfile,"#RMSE = %lf\n", sqrt(pow(10.0,((double) ResultErr))));
        }

        fprintf(qfile,"#END\n");

        if (strcmp(qfilename,"-")) fclose(qfile);

        if (doingCorrection)
        {
            if (GetActualBpp(job, qfilename, &actualBpp))
            {
                job->useCorrection = TRUE;
                job->correctionBpp = ResultBpp;
                job->addToTableBpp = actualBpp - ResultBpp;
                if ((!job->Silent) && (infile==0))
                {
                    fprintf(errfile,"Correction = %lf - %lf = %lf bpp\n",
                            actualBpp, ResultBpp, job->addToTableBpp);
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

