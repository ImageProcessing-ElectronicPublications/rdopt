
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



extern void InitOptimJob(OptimJob *job)
{
    int n,i;
    InitImage(&job->TheImage);
    job->NumTables = 1;
    job->opt_method.dp.BppMax = 1.0;
    job->opt_method.dp.BppScale = 5000;
    job->ThreshSpan = 0;
    job->MapQ = FALSE;
    job->VerboseLevel = 0;
    job->HistFilePresent = FALSE;
    job->ImFilePresent = FALSE;
    job->DumpStats = FALSE;
    job->OnlyDumpStats = FALSE;
    job->ReadCmdFromFile = FALSE;
    job->Silent = FALSE;
    job->DumpPlot = FALSE;
    job->PlotPoints = PLOT_POINTS;
    job->PlotBppMax = ((FFLOAT) 1.0);
    job->WeightedComps = FALSE;
    job->UseLagrangian = FALSE;
    job->useCorrection = FALSE;
    job->UseDCDPCM = FALSE;

    job->bppplane = -1;

    job->DCclamp = QTABENTRYMAX;
    for (i=0; i<MAXCOMPONENTS; i++)
    {
        job->WeightedCoefs[i] = FALSE;
        for (n=0; n<64; n++)
        {
            job->MinTable[i][n] = 1;
            job->MaxTable[i][n] = QTABENTRYMAX;
        }
    }

}


extern void GetParams(OptimJob *job, int argc, char *argv[])
{
    int i,j,n;
    int h,v;
    FFLOAT ftemp;
    double dtemp;
    FILE *fp;
    double b1,b2;
    char *tempstr;

    i=1;

    while (i<argc)
    {

        if (!strncmp(argv[i],"-stats",3))
        {
            job->DumpStats = TRUE;
        }
        else if (!strncmp(argv[i],"-mapq",4))
        {
            job->MapQ = TRUE;
        }
        else if (!strncmp(argv[i],"-dcdpcm",3))
        {
            job->UseDCDPCM = TRUE;
        }
        else if (!strncmp(argv[i],"-bppplane",5))
        {
            i++;
            if (i >= argc) BriefUsage();
            if ((job->bppplane = atoi(argv[i])) < 0)
            {
                job->bppplane = -1;
            }
        }
        else if (!strncmp(argv[i],"-silent",4))
        {
            job->Silent = TRUE;
            job->VerboseLevel = 0;
            job->TheImage.Silent = TRUE;
        }
        else if (!strncmp(argv[i],"-onlystats",6))
        {
            job->DumpStats = TRUE;
            job->OnlyDumpStats = TRUE;
        }
        else if (!strncmp(argv[i],"-height",4))
        {
            i++;
            if (i >= argc) BriefUsage();
            if ((job->TheImage.NumRows = atoi(argv[i])) <= 0)
            {
                FatalError("Image height must be > 0!");
            }
        }
        else if (!strncmp(argv[i],"-width",3))
        {
            i++;
            if (i >= argc) BriefUsage();
            if ((job->TheImage.NumCols = atoi(argv[i])) <= 0)
            {
                FatalError("Image width must be > 0!");
            }
        }
        else if (!strncmp(argv[i],"-planes",4))
        {
            i++;
            if (i >= argc) BriefUsage();
            if ((job->TheImage.NumComponents = atoi(argv[i])) <= 0)
            {
                FatalError("Number of color planes must be > 0!");
            }
            else if (job->TheImage.NumComponents > MAXCOMPONENTS)
            {
                FatalError("Number of color planes too high");
            }
        }
        else if (!strncmp(argv[i],"-subsamp",4))
        {
            if ((i+3) >= argc) BriefUsage();
            i++;
            n = atoi(argv[i]);
            if ((n >= MAXCOMPONENTS) || (n < 0))
                FatalError("Bad color plane number in -subsamp");
            i++;
            job->TheImage.OutSamplingFactor[n][0] = atoi(argv[i]);
            i++;
            job->TheImage.OutSamplingFactor[n][1] = atoi(argv[i]);

            if ((job->TheImage.OutSamplingFactor[n][0] <= 0) ||
                    (job->TheImage.OutSamplingFactor[n][1] <= 0))
            {
                FatalError("Sampling factors must be > 0!");
            }
        }
        else if (!strncmp(argv[i],"-insubsamp",4))
        {
            if ((i+3) >= argc) BriefUsage();
            i++;
            n = atoi(argv[i]);
            if ((n >= MAXCOMPONENTS) || (n < 0))
                FatalError("Bad color plane number in -inubsamp");
            i++;
            job->TheImage.InSamplingFactor[n][0] = atoi(argv[i]);
            i++;
            job->TheImage.InSamplingFactor[n][1] = atoi(argv[i]);

            if ((job->TheImage.InSamplingFactor[n][0] <= 0) ||
                    (job->TheImage.InSamplingFactor[n][1] <= 0))
            {
                FatalError("InSampling factors must be > 0!");
            }
        }
        else if (!strncmp(argv[i],"-v",2))
        {
            job->VerboseLevel++;
        }
        else if (!strncmp(argv[i],"-method",4))
        {
            i++;
            if (i>=argc) BriefUsage();
            if (!strncmp(argv[i],"lagrangian",4)) job->UseLagrangian = TRUE;
            else job->UseLagrangian = FALSE;
        }
        else if (!strncmp(argv[i],"-numtables",5))
        {
            i++;
            if (i>=argc) BriefUsage();
            if (((job->NumTables = atoi(argv[i])) <= 0)  ||
                    (job->NumTables > MAXCOMPONENTS))
            {
                FatalError("Bad value in -numtables");
            }
        }
        else if (!strncmp(argv[i],"-clampDC",5))
        {
            i++;
            if (i>=argc) BriefUsage();
            if ((job->DCclamp = atoi(argv[i])) < 0)
                FatalError("Bad value in -clampDC");
        }
        else if (!strncmp(argv[i],"-dontclampDC",6))
        {
            job->DCclamp = QTABENTRYMAX;
        }
        else if (!strncmp(argv[i],"-mintable",5))
        {
            if ((i+2) >= argc) BriefUsage();
            i++;
            n = atoi(argv[i]);
            if ((n >= MAXCOMPONENTS) || (n < 0))
                FatalError("Bad unit number in -mintable");
            i++;
            if (!strcmp(argv[i],"-")) fp = stdin;
            else if ((fp = fopen(argv[i],"r")) == NULL)
                FatalError("Could not open file for -mintable");
            ReadIntTable(fp,job->MinTable[n]);
            if (strcmp(argv[i],"-")) fclose(fp);
        }
        else if (!strncmp(argv[i],"-maxtable",5))
        {
            if ((i+2) >= argc) BriefUsage();
            i++;
            n = atoi(argv[i]);
            if ((n >= MAXCOMPONENTS) || (n < 0))
                FatalError("Bad unit number in -maxtable");
            i++;
            if (!strcmp(argv[i],"-")) fp = stdin;
            else if ((fp = fopen(argv[i],"r")) == NULL)
                FatalError("Could not open file for -maxtable");
            ReadIntTable(fp,job->MaxTable[n]);
            if (strcmp(argv[i],"-")) fclose(fp);
        }
        else if (!strncmp(argv[i],"-bppmax",5))
        {
            i++;
            if (i >= argc) BriefUsage();
            job->opt_method.dp.BppMax = ((FFLOAT) atof(argv[i]));
            if (job->opt_method.dp.BppMax <= 0.0)
                FatalError("Bad value in -bppmax");
        }
        else if (!strncmp(argv[i],"-pbppmax",5))
        {
            i++;
            if (i >= argc) BriefUsage();
            job->PlotBppMax = ((FFLOAT) atof(argv[i]));
            if (job->PlotBppMax <= 0.0)
                job->PlotBppMax = ((FFLOAT) 1.0);
        }
        else if (!strncmp(argv[i],"-correct",4))
        {
            i++;
            if (i >= argc) BriefUsage();
            job->correctionBpp = ((FFLOAT) atof(argv[i]));
            job->useCorrection = TRUE;
            if (job->correctionBpp <= 0.0)
                FatalError("Bad value in -correct");
        }
        else if (!strncmp(argv[i],"-bppscale",5))
        {
            i++;
            if (i >= argc) BriefUsage();
            job->opt_method.dp.BppScale = atoi(argv[i]);
            if (job->opt_method.dp.BppScale <= 1)
                FatalError("Bad value in -bppscale");
        }
        else if (!strncmp(argv[i],"-thresh",3))
        {
            i++;
            if (i >= argc) BriefUsage();
            job->ThreshSpan = atoi(argv[i]);
            if (job->ThreshSpan < 0)
                FatalError("Bad value in -thresh");
            if (job->ThreshSpan > QTABENTRYMAX)
                FatalError("Bad value (too high) in -thresh");
        }
        else if (!strncmp(argv[i],"-weights",4))
        {
            if ((i+2) >= argc) BriefUsage();
            i++;
            n = atoi(argv[i]);
            if ((n >= MAXCOMPONENTS) || (n < 0))
                FatalError("Bad unit number in -weights");
            i++;
            if (!strcmp(argv[i],"-")) fp = stdin;
            else if ((fp = fopen(argv[i],"r")) == NULL)
                FatalError("Could not open file for -weights");
            ReadRealTable(fp,job->CoefWeights[n]);
            if (strcmp(argv[i],"-")) fclose(fp);
            job->WeightedCoefs[n] = TRUE;
            /*** normalize weights ***/
            ftemp = 0.0;
            for (j=0; j<64; j++) ftemp += job->CoefWeights[n][j];
            ftemp /= 64.0;
            if (ftemp <= 0.0) FatalError("Bad weights in -weights");
            for (j=0; j<64; j++) job->CoefWeights[n][j] /= ftemp;
        }
        else if (!strncmp(argv[i],"-pweights",3))
        {
            i++;
            if (i >= argc) BriefUsage();
            job->WeightedComps = TRUE;
            tempstr = argv[i];
            sscanf(tempstr,"%lf",&dtemp);
            if (dtemp <= 0.0) FatalError("Bad weights in -pweights");
            job->CompWeights[0] = ((FFLOAT) dtemp);
            tempstr = strchr(tempstr,',');
            n = 1;
            while ((tempstr) && (n < MAXCOMPONENTS))
            {
                tempstr++;
                sscanf(tempstr,"%lf",&dtemp);
                if (dtemp <= 0.0) FatalError("Bad weights in -pweights");
                job->CompWeights[n] = ((FFLOAT) dtemp);
                n++;
                tempstr = strchr(tempstr,',');
            }
            while (n < MAXCOMPONENTS)
            {
                job->CompWeights[n] = ((FFLOAT) dtemp);
                n++;
            }
        }
        else if (!strncmp(argv[i],"-plot",3))
        {
            i++;
            if (i >= argc) BriefUsage();
            strcpy(job->PlotFileName,argv[i]);
            job->DumpPlot = TRUE;
        }
        else if (!strncmp(argv[i],"-points",3))
        {
            i++;
            if (i >= argc) BriefUsage();
            if ((job->PlotPoints = atoi(argv[i])) < 2)
                job->PlotPoints = PLOT_POINTS;
        }
        else if (!strncmp(argv[i],"-errfile",5))
        {
            i++;
            if (i >= argc) BriefUsage();
            if (!strcmp(argv[i],"-")) errfile = stdout;
            else
            {
                if ((errfile = fopen(argv[i],"w")) == NULL)
                {
                    errfile = stderr;
                }
            }
            job->TheImage.errfile = errfile;
        }
        else if (!strncmp(argv[i],"-help",4))
        {
            Usage();
        }
        else if (!strncmp(argv[i],"-cmdfile",4))
        {
            i++;
            if (i >= argc) BriefUsage();
            strcpy(job->CmdFileName,argv[i]);
            job->ReadCmdFromFile = TRUE;
        }

        else if (!strncmp(argv[i],"-hist",5))
        {
            i++;
            if (i >= argc) BriefUsage();
            job->HistFilePresent = TRUE;
            strcpy(job->HistFile,argv[i]);
        }


        else if (!strncmp(argv[i],"-im",3))
        {
            i++;
            if (i >= argc) BriefUsage();
            job->ImFilePresent = TRUE;
            if (strcmp(argv[i],"-")) strcpy(job->TheImage.ImFileName,argv[i]);
            i++;
            break;
        }

        else
        {
            BriefUsage();
        }

        i++;
    }

    if (job->HistFilePresent && job->ImFilePresent)
    {
        fprintf(errfile,"Cannot give both a histogram and an image as input!\n");
        exit(1);
    }
    if (!job->HistFilePresent && !job->ImFilePresent)
    {
        fprintf(errfile,"Must give a histogram or an image as input!\n");
        exit(1);
    }

    /* get image flags, if any */
    while (i<argc)
    {

        if (!strncmp(argv[i],"-pnm",4))
        {
            job->TheImage.ImKind = IM_PNM;
        }
        else if (!strncmp(argv[i],"-pgm",4))
        {
            job->TheImage.ImKind = IM_PGM;
        }
        else if (!strncmp(argv[i],"-ppm",4))
        {
            job->TheImage.ImKind = IM_PPM;
        }
        else if (!strncmp(argv[i],"-raw",4))
        {
            job->TheImage.ImKind = IM_RAW;
        }
        else if (!strncmp(argv[i],"-rgbtoycc",9))
        {
            job->TheImage.ColorConvKind = RGBtoYCC;
        }
        else if (!strncmp(argv[i],"-rgbto2ycc",10))
        {
            job->TheImage.ColorConvKind = RGBto2YCC;
        }
        else
        {
            BriefUsage();
        }
        i++;
    }

    if (job->ImFilePresent)
    {
        /***** decide image kind and fill derived info
         ***** in job->TheImage
        **********************************************/
        if (!PeekImage(&job->TheImage)) FatalError("Could not peek into image file");
    }

    /**********************************************************
     Fill derived information in job->...
    ***********************************************************/

    if (job->NumTables > job->TheImage.NumComponents)
    {
        if (!job->Silent)
            fprintf(errfile,"Number of units (%d) reset to number of color planes (%d)\n",
                    job->NumTables,job->TheImage.NumComponents);
        job->NumTables = job->TheImage.NumComponents;
    }



    for (i=0; i<job->TheImage.NumComponents; i++)
    {
        if (job->TheImage.OutSamplingFactor[i][0] <
                job->TheImage.InSamplingFactor[i][0])
        {

            if (!job->Silent)
                fprintf(errfile, "Resetting vertical samp factor for plane %d to %d\n",
                        i,job->TheImage.InSamplingFactor[i][0]);
            job->TheImage.OutSamplingFactor[i][0] =
                job->TheImage.InSamplingFactor[i][0];

        }
        if (job->TheImage.OutSamplingFactor[i][1] <
                job->TheImage.InSamplingFactor[i][1])
        {

            if (!job->Silent)
                fprintf(errfile, "Resetting horizontal samp factor for plane %d to %d\n",
                        i,job->TheImage.InSamplingFactor[i][1]);
            job->TheImage.OutSamplingFactor[i][1] =
                job->TheImage.InSamplingFactor[i][1];

        }
    }

    for (i=0; i<job->TheImage.NumComponents; i++)
    {
        v = job->TheImage.NumRows / (8*job->TheImage.OutSamplingFactor[i][0]);
        if ((v*8*job->TheImage.OutSamplingFactor[i][0]) !=
                job->TheImage.NumRows) v++;
        h = job->TheImage.NumCols / (8*job->TheImage.OutSamplingFactor[i][1]);
        if ((h*8*job->TheImage.OutSamplingFactor[i][1]) !=
                job->TheImage.NumCols) h++;
        job->NumBlocks[i] = ((FFLOAT) (h*v));
    }

    if ((job->bppplane >= 0) && (job->bppplane < job->TheImage.NumComponents))
    {
        job->NumPixelsForBpp =
            job->NumBlocks[job->bppplane]*((FFLOAT) 64.0);
    }
    else
    {
        job->bppplane = -1;
    }

    for (i=job->NumTables; i<job->TheImage.NumComponents; i++)
    {
        job->NumBlocks[job->NumTables-1] += job->NumBlocks[i];
    }

    ftemp = ((FFLOAT) log10(((double) MAXSAMPLE)*
                            ((double) MAXSAMPLE)));
    job->TotalNumPixels = 0.0;
    for (i=0; i<job->NumTables; i++)
    {
        job->NumPixels[i] = job->NumBlocks[i] * 64.0;
        job->LogPeakSq[i] = ftemp +
                            ((FFLOAT) log10((double) job->NumPixels[i]));
        job->TotalNumPixels += job->NumPixels[i];
    }
    job->LogTotalPeakSq = ftemp +
                          ((FFLOAT) log10((double) job->TotalNumPixels));


    if (job->bppplane == -1) job->NumPixelsForBpp = job->TotalNumPixels;

    for (i=0; i<job->NumTables; i++)
    {
        job->bppWeight[i] = job->NumPixels[i]/job->NumPixelsForBpp;
    }

    if (job->WeightedComps)
    {
        /* normalize weights */
        ftemp = ((FFLOAT) 0.0);
        for (n=0; n<job->NumTables; n++)
        {
            ftemp += job->CompWeights[n];
        }
        if (ftemp == 0.0) FatalError("Bad weights in -pweights");
        ftemp /= ((FFLOAT) job->NumTables);
        for (n=0; n<job->NumTables; n++)
        {
            job->CompWeights[n] /= ftemp;
        }
    }
    else
    {
        /* avoid a silly test later */
        job->WeightedComps = TRUE;
        for (n=0; n<job->NumTables; n++)
        {
            job->CompWeights[n] = ((FFLOAT) 1.0);
        }
    }

    if (job->UseLagrangian)
    {
        job->opt_method.lagr.ErrScale = job->NumPixelsForBpp;
        job->opt_method.lagr.NumHooks = NUM_HOOKS;
        if ((job->opt_method.lagr.BppHook = (FFLOAT *)
                                            calloc(1,sizeof(FFLOAT)*job->opt_method.lagr.NumHooks*3)) ==NULL)
            FatalError("GetParams out of memory");
        job->opt_method.lagr.LambdaHook = job->opt_method.lagr.BppHook +
                                          job->opt_method.lagr.NumHooks;
        job->opt_method.lagr.ErrHook = job->opt_method.lagr.LambdaHook +
                                       job->opt_method.lagr.NumHooks;
    }
    else
    {
        job->opt_method.dp.TableRowSize = RoundOff(job->opt_method.dp.BppMax * ((FFLOAT) job->opt_method.dp.BppScale)) + 2;
    }

    if (job->HistFilePresent)
    {
        if ((job->DumpStats) || (job->OnlyDumpStats))
        {
            if (!job->Silent)
                fprintf(errfile,"No image specified: -stats/-onlystats ignored\n");
            job->DumpStats = FALSE;
            job->OnlyDumpStats = FALSE;
        }
    }

}

extern void DumpJobChars(OptimJob *job, FILE *fp)
{
    int i;
    Image *Im;

    Im = &(job->TheImage);

    if (job->HistFilePresent)
    {
        fprintf(fp,"#Image statistics read from file HISTOGRAM\n");
    }
    else
    {
        fprintf(fp,"#Image read from %s file %s\n",
                ImKindString(Im->ImKind),Im->ImFileName);
    }

    fprintf(fp,"#Width x Height:%d x %d / ",Im->NumCols,Im->NumRows);
    for (i=0; i<Im->NumComponents-1; i++)
    {
        fprintf(fp, "(%d x %d):",Im->OutSamplingFactor[i][1],
                Im->OutSamplingFactor[i][0]);
    }
    fprintf(fp, "(%d x %d)\n",Im->OutSamplingFactor[Im->NumComponents-1][1],
            Im->OutSamplingFactor[Im->NumComponents-1][0]);


    fprintf(fp,"#Color-conversion applied: ");
    if (Im->ColorConvKind == RGBtoYCC) fprintf(fp,"RGB to YCC\n");
    else if (Im->ColorConvKind == RGBto2YCC)
        fprintf(fp,"RGB to YCC with subsamp for chrominance\n");
    else fprintf(fp,"none\n");

    fprintf(fp,"#Number of qtables used: %d\n",job->NumTables);
    if (job->UseLagrangian)
    {
        fprintf(fp,"#Optimization parameters (method: Lagrangian): \n");
        fprintf(fp,"#\t DC clamp = %d, ThreshSpan = %d\n",
                job->DCclamp, job->ThreshSpan);
    }
    else
    {
        fprintf(fp,"#Optimization parameters (method: Dynamic Programming): \n");
        fprintf(fp,"#\tMaxBpp = %lf, BppScale = %d, DC clamp = %d, ThreshSpan = %d\n",
                (double) job->opt_method.dp.BppMax,
                job->opt_method.dp.BppScale, job->DCclamp,
                job->ThreshSpan);
    }
    if (job->useCorrection)
    {
        fprintf(fp,"#\tThe predicted Bpp is corrected by adding %lf in the end,\n",job->addToTableBpp);
        fprintf(fp,"#\t\twhich is the error at asking bpp = %lf\n",job->correctionBpp);
    }
    fprintf(fp,"#\n");
}





extern void WriteDescription(char *s, Image *Im, int unum1, int unum2)
{
    sprintf(s,"Image %s color planes %d to %d\0",
            (!strcmp(Im->ImFileName,"-"))?"<stdin>":
            Im->ImFileName,unum1,unum2);
}
