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
#include "qentry.h"


static int ScanRowForErr(FFLOAT *row, FFLOAT err, int least, int most)
{
    int i;
    for (i=least; i<=most; i++)
    {
        if (row[i] < err) return(i);
    }
    return(-1);
}

static int ScanRowForBpp(FFLOAT *row, int bpp, int least, int most)
{
    int i;
    if (((i=bpp) < least) || (i > most)) return(-1);
    while (i>=least)
    {
        if (row[i] < INFINITY) return(i);
        i--;
    }
    return(-1);
}

static void GetStartingPoints(int *ToCombine[], int *StartingPt,
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


static void RecoverQ(int Loc, int unum, OptimJob *job, int *Q)
{
    int n, NextLoc;
    NextLoc = Loc;
    for (n=0; n<64; n++)
    {
        Q[n] = (int) job->QChoice[unum][n][NextLoc];

        NextLoc -= job->Bpp[unum][n][Q[n]];
        if (job->MapQ) Q[n] = MapQentry(n,Q[n]);
    }
}


static void ListCommands(void)
{
    fprintf(stderr,"Valid commands are:\n");
    fprintf(stderr,"\t size <target_size> [<unit_num>]\n");
    fprintf(stderr,"\t bpp <target_bpp> [<unit_num>]\n");
    fprintf(stderr,"\t psnr <target_psnr> [<unit_num>]\n");
    fprintf(stderr,"\t snr <target_snr> [<unit_num>]\n");
    fprintf(stderr,"\t rmse <target_snr> [<unit_num>]\n");
    fprintf(stderr,"\t qfile fname\n");
    fprintf(stderr,"A file RDOPT.Q.<command>.<target>[.<unit_num>] will be generated\n");
    fprintf(stderr,"  (unless a qfile command has been issued after the last command)\n");
}

static void BriefListCommands(void)
{
    fprintf(stderr,"Valid commands: size/bpp/psnr/snr/rmse target [unit]\n");
    fprintf(stderr,"                help/qfile fname/quit\n");
}

extern void Command(OptimJob *job)
{
    int IntTargetBpp, ResultSize;
    FFLOAT ResultBpp, ResultSNR, ResultPSNR, TargetErr, RealTargetBpp;
    FFLOAT ResultErr, temp;
    double arg;
    int Q[64];
    char command[STRLENMAX], nextline[STRLENMAX], qfilename[STRLENMAX];
    boolean AllUnits;
    boolean TargetIsBpp;
    int unum, StartingPt[MAXCOMPONENTS], numread, Loc;
    FILE *qfile;
    int infile;

    if ((job->ReadCmdFromFile) && (strcmp(job->CmdFileName,"-")))
    {
        infile = open(job->CmdFileName,O_RDONLY,0);
        if (infile < 0) FatalError("Could not open command file");
    }
    else infile = 0;

    if ((!job->Silent) && (infile==0))
        fprintf(stderr,"\n**** RD-OPT Command Interface ****\n\n");

    qfilename[0] = '\0';

    while (TRUE)
    {
        if ((!job->Silent) && (infile==0)) fprintf(stderr,"Command> ");
        if (newlinefd(nextline,infile)==EOF) break;

        if (!strncmp(nextline,"qfile",5))
        {
            sscanf(nextline,"%s%s",command,qfilename);
            continue;
        }

        numread = sscanf(nextline,"%s%lf%d",command,&arg,&unum);
        if (numread <= 0) continue;

        AllUnits = TRUE;
        if (numread==3)
        {
            if ((unum < 0) || (unum >= job->NumTables))
            {
                fprintf(stderr,"Unit_num must be in the range %d - %d\n",
                        0,job->NumTables-1);
                continue;
            }
            if (job->NumTables > 1)
            {
                AllUnits = FALSE;
            }
        }


        if (!strncmp(command,"size",4))
        {
            if (numread < 2)
            {
                fprintf(stderr,"Must specify target size in bytes\n");
                continue;
            }
            TargetIsBpp = TRUE;
            RealTargetBpp = ((FFLOAT) (arg*8.0));
            if (AllUnits) RealTargetBpp /= job->TotalNumPixels;
            else
            {
                RealTargetBpp /= job->NumPixels[unum];
                RealTargetBpp *= job->bppWeight[unum];
            }
            RealTargetBpp *= ((FFLOAT) job->BppScale);
            IntTargetBpp = RoundOff(RealTargetBpp);

        }
        else if (!strncmp(command,"bpp",3))
        {
            if (numread < 2)
            {
                fprintf(stderr,"Must specify target bits-per-pixel\n");
                continue;
            }
            TargetIsBpp = TRUE;
            RealTargetBpp = ((FFLOAT) arg);
            if (!AllUnits)
            {
                RealTargetBpp *= job->bppWeight[unum];
            }
            RealTargetBpp *= ((FFLOAT) job->BppScale);
            IntTargetBpp = RoundOff(RealTargetBpp);
        }
        else if (!strncmp(command,"psnr",3))
        {
            if (numread < 2)
            {
                fprintf(stderr,"Must specify target PSNR in dB\n");
                continue;
            }
            TargetIsBpp = FALSE;
            if (AllUnits)
            {
                TargetErr = job->LogTotalPeakSq;
            }
            else
            {
                TargetErr = job->LogPeakSq[unum];
            }
            TargetErr -= ((FFLOAT) arg/10.0);
        }
        else if (!strncmp(command,"snr",3))
        {
            if (numread < 2)
            {
                fprintf(stderr,"Must specify target SNR in dB\n");
                continue;
            }
            TargetIsBpp = FALSE;
            if (AllUnits)
            {
                TargetErr = job->LogTotalSignalSq;
            }
            else
            {
                TargetErr = job->LogSignalSq[unum];
            }
            TargetErr -= ((FFLOAT) arg/10.0);

        }
        else if (!strncmp(command,"rmse",3))
        {
            if (numread < 2)
            {
                fprintf(stderr,"Must specify target RMSE\n");
                continue;
            }
            TargetIsBpp = FALSE;
            TargetErr = ((FFLOAT) log10((double) arg))*2.0;
            if (AllUnits)
            {
                TargetErr += ((FFLOAT) log10((double) job->TotalNumPixels));
            }
            else
            {
                TargetErr += ((FFLOAT) log10((double) job->NumPixels[unum]));
            }
        }
        else if (!strncmp(command,"help",4))
        {
            ListCommands();
            continue;
        }
        else if (!strncmp(command,"quit",4))
        {
            break;
        }
        else
        {
            BriefListCommands();
            continue;
        }


        /******** extract Qtable(s) **********/
        if (AllUnits)
        {
            if (job->NumTables==1)
            {
                if (TargetIsBpp)
                {
                    Loc = ScanRowForBpp(job->LastRow[0],IntTargetBpp,
                                        job->LeastLast[0], job->MostLast[0]);
                }
                else
                {
                    Loc = ScanRowForErr(job->LastRow[0],TargetErr,
                                        job->LeastLast[0], job->MostLast[0]);
                }
                if (Loc < 0)
                {
                    fprintf(stderr,"Set target cannot be achieved\n");
                    continue;
                }
            }
            else
            {
                if (TargetIsBpp)
                {
                    Loc = ScanRowForBpp(job->CombinedLastRow,IntTargetBpp,
                                        job->CombinedLeast, job->CombinedMost);
                }
                else
                {
                    Loc = ScanRowForErr(job->CombinedLastRow,TargetErr,
                                        job->CombinedLeast, job->CombinedMost);
                }
                if (Loc < 0)
                {
                    fprintf(stderr,"Set target cannot be achieved\n");
                    continue;
                }
                GetStartingPoints(job->ToCombine, StartingPt,
                                  Loc, job->NumTables);
            }
        }
        else
        {
            if (TargetIsBpp)
            {
                Loc = ScanRowForBpp(job->LastRow[unum],IntTargetBpp,
                                    job->LeastLast[unum], job->MostLast[unum]);
            }
            else
            {
                Loc = ScanRowForErr(job->LastRow[unum],TargetErr,
                                    job->LeastLast[unum], job->MostLast[unum]);
            }
            if (Loc < 0)
            {
                fprintf(stderr,"Set target cannot be achieved\n");
                continue;
            }
        }

        /******** open qfile ***********/
        if (qfilename[0] == '\0')
        {
            if (AllUnits)
                sprintf(qfilename,"RDOPT.Q.%s.%.3lf\0",command,arg);
            else
                sprintf(qfilename,"RDOPT.Q.%s.%.3lf.%d\0",command,arg,unum);
        }

        if (!strcmp(qfilename,"-")) qfile = stdout;
        else
        {
            if ((qfile = fopen(qfilename,"w")) == NULL)
            {
                fprintf(stderr,"Could not open qtable file %s for writing\n",
                        qfilename);
                continue;
            }
        }

        fprintf(qfile,"#%s\n",CQTAB);
        DumpJobChars(job, qfile);

        if (AllUnits)
        {
            if (job->NumTables==1)
            {
                fprintf(qfile,"#Quantization table for color planes %d thru %d\n",0,job->TheImage.NumComponents-1);
                RecoverQ(Loc, 0, job, Q);
                WriteIntTable(qfile, Q);

                ResultErr = job->LastRow[0][Loc];
                ResultBpp = ((FFLOAT) Loc)/((FFLOAT) job->BppScale);
                temp = ResultBpp * job->TotalNumPixels / 8.0;
                ResultSize = RoundOff(temp);
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
                    RecoverQ(StartingPt[unum], unum, job, Q);
                    WriteIntTable(qfile, Q);

                    ResultErr = job->LastRow[unum][StartingPt[unum]];
                    ResultBpp = ((FFLOAT) StartingPt[unum])/((FFLOAT) job->BppScale);
                    ResultBpp /= job->bppWeight[unum];
                    temp = ResultBpp * job->NumPixels[unum] / 8.0;
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
                ResultErr = job->CombinedLastRow[Loc];
                ResultBpp = ((FFLOAT) Loc)/((FFLOAT) job->BppScale);
                temp = ResultBpp * job->TotalNumPixels / 8.0;
                ResultSize = RoundOff(temp);
                ResultSNR = (job->LogTotalSignalSq - ResultErr)*10.0;
                ResultPSNR = (job->LogTotalPeakSq - ResultErr)*10.0;

                fprintf(qfile,"#Bpp = %lf bits per pixel\n", ((double) ResultBpp));
                fprintf(qfile,"#Size = %d bytes\n", ResultSize);
                fprintf(qfile,"#SNR = %lf dB\n", ((double) ResultSNR));
                fprintf(qfile,"#PSNR = %lf dB\n", ((double) ResultPSNR));
                ResultErr -= ((FFLOAT) log10((double) job->TotalNumPixels));
                fprintf(qfile,"#RMSE = %lf\n", sqrt(pow(10.0,((double) ResultErr))));
            }
        }
        else
        {
            fprintf(qfile,"# THIS FILE DOES NOT HAVE TABLES FOR ALL UNITS\n");
            fprintf(qfile,"#Quantization table for color planes %d thru %d\n",unum,((unum==(job->NumTables-1))?(job->TheImage.NumComponents-1):unum));
            RecoverQ(Loc, unum, job, Q);
            WriteIntTable(qfile, Q);

            ResultErr = job->LastRow[unum][Loc];
            ResultBpp = ((FFLOAT) Loc)/((FFLOAT) job->BppScale);
            ResultBpp /= job->bppWeight[unum];
            temp = ResultBpp * job->NumPixels[unum] / 8.0;
            ResultSize = RoundOff(temp);
            ResultSNR = (job->LogSignalSq[unum] - ResultErr)*10.0;
            ResultPSNR = (job->LogPeakSq[unum] - ResultErr)*10.0;

            fprintf(qfile,"#Bpp = %lf bits per pixel\n", ((double) ResultBpp));
            fprintf(qfile,"#Size = %d bytes\n", ResultSize);
            fprintf(qfile,"#SNR = %lf dB\n", ((double) ResultSNR));
            fprintf(qfile,"#PSNR = %lf dB\n", ((double) ResultPSNR));
            ResultErr -= ((FFLOAT) log10((double) job->NumPixels[unum]));
            fprintf(qfile,"#RMSE = %lf\n", sqrt(pow(10.0,((double) ResultErr))));

        }

        fprintf(qfile,"#END\n");

        if (strcmp(qfilename,"-")) fclose(qfile);
        qfilename[0] = '\0';



    }
    if ((job->ReadCmdFromFile) && (strcmp(job->CmdFileName,"-")))
        close(infile);
}
