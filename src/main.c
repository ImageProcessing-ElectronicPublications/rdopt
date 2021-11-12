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

OptimJob TheJob;

int main(int argc, char *argv[])
{
    int UnitNum, EndCompNum, n, i;
    FILE *hfile;
    char DesString[STRLENMAX];

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
            fprintf(stderr,"Unit %d: color planes %d through %d\n",
                    UnitNum,UnitNum,EndCompNum);
        }

        /************** get stats in TheJob.Histogram ***********/
        if (TheJob.HistFilePresent)
        {
            if (TheJob.VerboseLevel)
            {
                fprintf(stderr,"\tReading statistics from the file %s\n",
                        TheJob.HistFile);
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
            for (i=UnitNum; i<=EndCompNum; i++)
            {
                if (TheJob.VerboseLevel)
                {
                    fprintf(stderr,"\tReading image color plane #%d\n",i);
                }
                ReadImgComp(&TheJob.TheImage, i);
                if (TheJob.VerboseLevel)
                {
                    fprintf(stderr,"\tComputing statistics for image color plane #%d\n",i);
                }
                SetHistogram(TheJob.Histogram,&TheJob.TheImage,i,
                             &TheJob.LogSignalSq[UnitNum]);
                FreeImgComp(&TheJob.TheImage,i);
            }
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
                fprintf(stderr,"\tDumping statistics\n");
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


        if (TheJob.BppRangeExists[UnitNum])
        {
            if (TheJob.VerboseLevel)
            {
                fprintf(stderr,"\tTranslating Bpp-ranges to min-max entries\n");
            }
            TranslateBppRange(&TheJob,UnitNum);
        }

        /*********** fill ERR and BPP **********************/
        PrepareForErrBpp(&TheJob,UnitNum);

        if (TheJob.VerboseLevel)
        {
            fprintf(stderr,"\tFilling errors\n");
        }
        SetErr(&TheJob,UnitNum);

        if (TheJob.VerboseLevel)
        {
            fprintf(stderr,"\tFilling entropies\n");
        }
        SetBpp(&TheJob,UnitNum);

        FreeHistogram(TheJob.Histogram);

        if (TheJob.VerboseLevel)
        {
            fprintf(stderr,"\tRunning optimization algorithm\n");
        }


        Optimize(&TheJob,UnitNum);

        free(TheJob.Err[0]);

    }

    if (TheJob.OnlyDumpStats) exit(0);

    /*********************** combine units *****************/
    if (TheJob.NumTables > 1)
    {

        if (TheJob.VerboseLevel)
        {
            fprintf(stderr,"Combining information from all units\n");
        }

        CombineUnits(&TheJob);
    }

    Epilogue(&TheJob);

    /****************** command interpreter ******************/

    Command(&TheJob);


    return(0);
}
