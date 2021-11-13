
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



OptimJob TheJob;

FILE *errfile; /* global stderr */

int main(int argc, char *argv[])
{
    int UnitNum, EndCompNum, n, i;
    FILE *hfile;
    char DesString[STRLENMAX];
    FFLOAT Stemp[64];

    errfile = stderr;

    InitOptimJob(&TheJob);

    /***** Read options, fill derived info in TheJob *********/
    GetParams(&TheJob,argc,argv);

    for (UnitNum=0; UnitNum < TheJob.NumTables; UnitNum++)
    {

        EndCompNum = (UnitNum == (TheJob.NumTables-1)) ?
                     TheJob.TheImage.NumComponents - 1 :
                     UnitNum;
        if (TheJob.VerboseLevel)
        {
            fprintf(errfile,"Unit %d: color planes %d through %d\n",
                    UnitNum,UnitNum,EndCompNum);
            fflush(errfile);
        }

        /************** get stats in TheJob.Histogram ***********/
        if (TheJob.HistFilePresent)
        {
            if (TheJob.VerboseLevel)
            {
                fprintf(errfile,"\tReading statistics from the file %s\n",
                        TheJob.HistFile);
                fflush(errfile);
            }
            if (UnitNum==0)
            {
                if ((hfile = fopen(TheJob.HistFile,"r")) == NULL)
                {
                    FatalError("Could not open Histogram file");
                }
            }
            InitHistogram(TheJob.Histogram);
            ReadHistogram(hfile,TheJob.Histogram,
                          &TheJob.LogSignalSq[UnitNum]);
            TheJob.LogSignalSq[UnitNum] *= TheJob.NumPixels[UnitNum];
            SetSignal(TheJob.Histogram,Stemp);
            if (TheJob.WeightedCoefs[UnitNum])
            {
                for (i=0; i<64; i++)
                    Stemp[i] *= TheJob.CoefWeights[UnitNum][i];
            }
            SigToRemainingSig(Stemp, TheJob.RemainingSignal[UnitNum]);

            if (UnitNum == (TheJob.NumTables-1)) fclose(hfile);

            if (NotEqualReal(TheJob.NumBlocks[UnitNum],
                             GetTotalBlocks(TheJob.Histogram)))
            {
                FatalError("HISTOGRAM file does not match image dimensions");
            }
        }
        else
        {
            InitHistogram(TheJob.Histogram);
            TheJob.LogSignalSq[UnitNum] = 0.0;
            for (i=0; i<64; i++) Stemp[i] = ((FFLOAT) 0.0);
            for (i=UnitNum; i<=EndCompNum; i++)
            {
                if (TheJob.VerboseLevel)
                {
                    fprintf(errfile,"\tReading image color plane #%d\n",i);
                    fflush(errfile);
                }
                if (!ReadImgComp(&TheJob.TheImage, i))
                    FatalError("Could not read image plane");
                if (TheJob.VerboseLevel)
                {
                    fprintf(errfile,"\tComputing statistics for image color plane #%d\n",i);
                    fflush(errfile);
                }
                SetHistogram(TheJob.Histogram,&TheJob.TheImage,i,
                             &TheJob.LogSignalSq[UnitNum], Stemp, TheJob.UseDCDPCM);
                FreeImgComp(&TheJob.TheImage,i);
            }
            if (TheJob.WeightedCoefs[UnitNum])
            {
                for (i=0; i<64; i++)
                    Stemp[i] *= TheJob.CoefWeights[UnitNum][i];
            }
            SigToRemainingSig(Stemp, TheJob.RemainingSignal[UnitNum]);
        }

        for (n=0; n<64; n++)
        {
            while (!TheJob.Histogram[n].PlusCount[TheJob.Histogram[n].PlusSize-1])
                TheJob.Histogram[n].PlusSize--;
            while (!TheJob.Histogram[n].MinusCount[TheJob.Histogram[n].MinusSize-1])
                TheJob.Histogram[n].MinusSize--;
        }

        if (TheJob.DumpStats)
        {
            if (UnitNum == 0)
            {
                if ((hfile = fopen("HISTOGRAM","w")) == NULL)
                {
                    FatalError("Could not open file HISTOGRAM for writing");
                }
            }
            WriteDescription(DesString, &TheJob.TheImage,UnitNum,EndCompNum);
            if (TheJob.VerboseLevel)
            {
                fprintf(errfile,"\tDumping statistics\n");
                fflush(errfile);
            }
            DumpHistogram(hfile,TheJob.Histogram,DesString,
                          TheJob.LogSignalSq[UnitNum]/TheJob.NumPixels[UnitNum]);
            if (UnitNum==(TheJob.NumTables-1)) fclose(hfile);
            if (TheJob.OnlyDumpStats)
            {
                FreeHistogram(TheJob.Histogram);
                continue;
            }
        }



        /*********** fill ERR and BPP **********************/

        if (!TheJob.UseLagrangian)
        {
            PrepareForErrBpp(&TheJob,UnitNum);

            if (TheJob.VerboseLevel)
            {
                fprintf(errfile,"\tFilling errors and entropies\n");
                fflush(errfile);
            }
            SetErrBpp(&TheJob,UnitNum);


            FreeHistogram(TheJob.Histogram);

            if (TheJob.VerboseLevel)
            {
                fprintf(errfile,"\tRunning optimization algorithm\n");
                fflush(errfile);
            }


            Optimize(&TheJob,UnitNum);
        }
        else
        {

            lagrPrepareForErrBpp(&TheJob,UnitNum);

            if (TheJob.VerboseLevel)
            {
                fprintf(errfile,"\tFilling errors and entropies\n");
                fflush(errfile);
            }
            lagrSetErrBpp(&TheJob,UnitNum);


            FreeHistogram(TheJob.Histogram);

            if (TheJob.VerboseLevel)
            {
                fprintf(errfile,"\tSorting R-D tables\n");
                fflush(errfile);
            }

            lagrSortErrBpp(&TheJob,UnitNum);

        }

    }

    if (TheJob.OnlyDumpStats) exit(0);

    if (!TheJob.UseLagrangian)
    {
        /*********************** combine units *****************/
        if (TheJob.NumTables > 1)
        {

            if (TheJob.VerboseLevel)
            {
                fprintf(errfile,"Combining information from all units\n");
                fflush(errfile);
            }

            CombineUnits(&TheJob);
        }

        Epilogue(&TheJob);
        /****************** command interpreter ******************/
        Command(&TheJob);
    }
    else
    {
        lagrEpilogue(&TheJob);
        /****************** command interpreter ******************/
        lagrCommand(&TheJob);
    }


    if (((int) errfile != (int) stderr) && ((int) errfile != (int) stdout)) fclose(errfile);

    return(0);
}
