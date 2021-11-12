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



extern void InitOptimJob(OptimJob *job)
{
    int n,i;
    InitImage(&job->TheImage);
    job->NumTables = 1;
    job->BppMax = 1.0;
    job->BppScale = 5000;
    job->MapQ = FALSE;
    job->WeightedCoefs = FALSE;
    job->VerboseLevel = 0;
    job->HistFilePresent = FALSE;
    job->ImFilePresent = FALSE;
    job->DumpStats = FALSE;
    job->OnlyDumpStats = FALSE;
    job->ReadCmdFromFile = FALSE;
    job->Silent = FALSE;


    job->DCclamp = 12;
    for (i=0; i<MAXCOMPONENTS; i++)
    {
        job->BppRangeExists[i] = FALSE;
        for (n=0; n<64; n++)
        {
            job->MinTable[i][n] = 1;
            job->MaxTable[i][n] = QTABENTRYMAX;
        }
    }

    job->DumpPlot = FALSE;
}


extern void GetParams(OptimJob *job, int argc, char *argv[])
{
    int i,j,n;
    int h,v;
    FFLOAT ftemp;
    FILE *fp;
    double b1,b2;

    i=1;

    while (i<argc)
    {

        if (!strncmp(argv[i],"-stats",6))
        {
            job->DumpStats = TRUE;
        }
        else if (!strncmp(argv[i],"-mapq",5))
        {
            job->MapQ = TRUE;
        }
        else if (!strncmp(argv[i],"-silent",7))
        {
            job->Silent = TRUE;
            job->VerboseLevel = 0;
            job->TheImage.Silent = TRUE;
        }
        else if (!strncmp(argv[i],"-onlystats",10))
        {
            job->DumpStats = TRUE;
            job->OnlyDumpStats = TRUE;
        }
        else if (!strncmp(argv[i],"-height",7))
        {
            i++;
            if (i >= argc) BriefUsage();
            if ((job->TheImage.NumRows = atoi(argv[i])) <= 0)
            {
                FatalError("Image height must be > 0!");
            }
        }
        else if (!strncmp(argv[i],"-width",6))
        {
            i++;
            if (i >= argc) BriefUsage();
            if ((job->TheImage.NumCols = atoi(argv[i])) <= 0)
            {
                FatalError("Image width must be > 0!");
            }
        }
        else if (!strncmp(argv[i],"-planes",6))
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
        else if (!strncmp(argv[i],"-subsamp",8))
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
        else if (!strncmp(argv[i],"-insubsamp",10))
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
        else if (!strncmp(argv[i],"-numtables",10))
        {
            i++;
            if (i>=argc) BriefUsage();
            if (((job->NumTables = atoi(argv[i])) <= 0)  ||
                    (job->NumTables > MAXCOMPONENTS))
            {
                FatalError("Bad value in -numtables");
            }
        }
        else if (!strncmp(argv[i],"-clampDC",8))
        {
            i++;
            if (i>=argc) BriefUsage();
            if ((job->DCclamp = atoi(argv[i])) < 0)
                FatalError("Bad value in -clampDC");
        }
        else if (!strncmp(argv[i],"-dontclampDC",12))
        {
            job->DCclamp = QTABENTRYMAX;
        }
        else if (!strncmp(argv[i],"-bppdist",8))
        {
            if ((i+2) >= argc) BriefUsage();
            i++;
            n = atoi(argv[i]);
            if ((n >= MAXCOMPONENTS) || (n < 0))
                FatalError("Bad unit number in -bppdist");
            i++;
            if ((fp = fopen(argv[i],"r")) == NULL)
                FatalError("Could not open file for -bppdist");
            for (j=0; j<64; j++)
            {
                if (fscanf(fp,"%lf:%lf",&b1,&b2) != 2)
                    FatalError("Bad format in file for -bppdist");
                job->BppRange[n][j][0] = b1;
                job->BppRange[n][j][1] = b2;
            }
            fclose(fp);
            job->BppRangeExists[n] = TRUE;
        }
        else if (!strncmp(argv[i],"-mintable",9))
        {
            if ((i+2) >= argc) BriefUsage();
            i++;
            n = atoi(argv[i]);
            if ((n >= MAXCOMPONENTS) || (n < 0))
                FatalError("Bad unit number in -mintable");
            i++;
            if ((fp = fopen(argv[i],"r")) == NULL)
                FatalError("Could not open file for -mintable");
            ReadIntTable(fp,job->MinTable[n]);
            fclose(fp);
        }
        else if (!strncmp(argv[i],"-maxtable",9))
        {
            if ((i+2) >= argc) BriefUsage();
            i++;
            n = atoi(argv[i]);
            if ((n >= MAXCOMPONENTS) || (n < 0))
                FatalError("Bad unit number in -maxtable");
            i++;
            if ((fp = fopen(argv[i],"r")) == NULL)
                FatalError("Could not open file for -maxtable");
            ReadIntTable(fp,job->MaxTable[n]);
            fclose(fp);
        }
        else if (!strncmp(argv[i],"-bppmax",7))
        {
            i++;
            if (i >= argc) BriefUsage();
            job->BppMax = ((FFLOAT) atof(argv[i]));
            if (job->BppMax <= 0.0)
                FatalError("Bad value in -bppmax");
        }
        else if (!strncmp(argv[i],"-bppscale",9))
        {
            i++;
            if (i >= argc) BriefUsage();
            job->BppScale = atoi(argv[i]);
            if (job->BppScale <= 0)
                FatalError("Bad value in -bppscale");
        }
        else if (!strncmp(argv[i],"-weights",8))
        {
            i++;
            if (i >= argc) BriefUsage();
            if ((fp = fopen(argv[i],"r")) == NULL)
                FatalError("Could not open file for -weights");
            ReadRealTable(fp,job->CoefWeights);
            fclose(fp);
            job->WeightedCoefs = TRUE;
            /*** normalize weights ***/
            ftemp = 0.0;
            for (j=0; j<64; j++) ftemp += job->CoefWeights[j];
            if (ftemp <= 0.0) FatalError("Bad weights in -weights");
            for (j=0; j<64; j++) job->CoefWeights[j] /= ftemp;
        }
        else if (!strncmp(argv[i],"-plot",5))
        {
            i++;
            if (i >= argc) BriefUsage();
            strcpy(job->PlotFileName,argv[i]);
            job->DumpPlot = TRUE;
        }
        else if (!strncmp(argv[i],"-help",5))
        {
            Usage();
        }
        else if (!strncmp(argv[i],"-cmdfile",8))
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
        fprintf(stderr,"Cannot give both a histogram and an image as input!\n");
        exit(1);
    }
    if (!job->HistFilePresent && !job->ImFilePresent)
    {
        fprintf(stderr,"Must give a histogram or an image as input!\n");
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
        PeekImage(&job->TheImage);
    }

    /**********************************************************
     Fill derived information in job->...
    ***********************************************************/

    if (job->NumTables > job->TheImage.NumComponents)
    {
        if (!job->Silent)
            fprintf(stderr,"Number of units (%d) reset to number of color planes (%d)\n",job->NumTables,job->TheImage.NumComponents);
        job->NumTables = job->TheImage.NumComponents;
    }

    job->TableRowSize = RoundOff(job->BppMax * ((FFLOAT) job->BppScale)) + 2;


    for (i=0; i<job->TheImage.NumComponents; i++)
    {
        if (job->TheImage.OutSamplingFactor[i][0] <
                job->TheImage.InSamplingFactor[i][0])
        {

            if (!job->Silent)
                fprintf(stderr, "Resetting vertical samp factor for plane %d to %d\n",
                        i,job->TheImage.InSamplingFactor[i][0]);
            job->TheImage.OutSamplingFactor[i][0] =
                job->TheImage.InSamplingFactor[i][0];

        }
        if (job->TheImage.OutSamplingFactor[i][1] <
                job->TheImage.InSamplingFactor[i][1])
        {

            if (!job->Silent)
                fprintf(stderr, "Resetting horizontal samp factor for plane %d to %d\n",
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

    if (job->NumTables == 1)
    {
        job->bppWeight[0] = 1.0;
    }
    else
    {
        for (i=0; i<job->NumTables; i++)
        {
            job->bppWeight[i] = job->NumPixels[i]/job->TotalNumPixels;
        }
    }

    if (job->HistFilePresent)
    {
        if ((job->DumpStats) || (job->OnlyDumpStats))
        {
            if (!job->Silent)
                fprintf(stderr,"No image specified: -stats/-onlystats ignored\n");
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
    fprintf(fp,"#Optimization parameters: \n");
    fprintf(fp,"#\tMaxBpp = %lf, BppScale = %d, DC clamp = %d\n",
            (double) job->BppMax, job->BppScale, job->DCclamp);

    fprintf(fp,"#\n");
}





extern void WriteDescription(char *s, Image *Im, int unum1, int unum2)
{
    sprintf(s,"Image %s color planes %d to %d\0",
            (!strcmp(Im->ImFileName,"-"))?"<stdin>":
            Im->ImFileName,unum1,unum2);
}
