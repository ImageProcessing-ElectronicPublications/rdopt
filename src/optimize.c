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

static void ResetRow(FFLOAT *row, int width)
{
    int i;
    for (i=0; i<width; i++) row[i] = INFINITY;
}

extern void Optimize(OptimJob *job, int unum)
{
    FFLOAT *CurrRow, *PrevRow, *TempRow;
    int LeastLast, MostLast, i, Loc, n, q, k;
    FFLOAT ErrTemp;


    /****** allocate ******/
    if ((CurrRow = (FFLOAT *)
                   calloc(1,job->TableRowSize*sizeof(FFLOAT)))==NULL)
    {
        FatalError("Optimize out of memory");
    }
    if ((PrevRow = (FFLOAT *)
                   calloc(1,job->TableRowSize*sizeof(FFLOAT)))==NULL)
    {
        FatalError("Optimize out of memory");
    }

    ResetRow(CurrRow,job->TableRowSize);
    ResetRow(PrevRow,job->TableRowSize);

    if ((job->QChoice[unum][0] = (Qentry *)
                                 calloc(1,job->TableRowSize*64*sizeof(Qentry)))==NULL)
    {
        FatalError("Optimize out of memory");
    }

    for (n=1; n<64; n++)
    {
        job->QChoice[unum][n] = job->QChoice[unum][n-1] +
                                job->TableRowSize;
    }

    if (job->VerboseLevel > 2)
    {
        fprintf(stderr,"\t\tCalculating row for coeff #%d\n",63);
    }

    LeastLast = job->TableRowSize-1;
    MostLast = 0;
    for (q=job->MaxTable[unum][63]; q>=job->MinTable[unum][63]; q--)
    {
        if ((Loc = job->Bpp[unum][63][q]) >= job->TableRowSize) continue;
        if (job->Err[63][q] < PrevRow[Loc])
        {
            if (Loc < LeastLast) LeastLast = Loc;
            if (Loc > MostLast) MostLast = Loc;
            job->QChoice[unum][63][Loc] = ((Qentry) q);
            PrevRow[Loc] = job->Err[63][q];
        }
    }
    if (job->VerboseLevel > 3)
    {
        fprintf(stderr,"\t\t\tBPP range: %d/%d to %d/%d\n",
                LeastLast,job->BppScale,MostLast,job->BppScale);
    }

    for (n=62; n>=0; n--)
    {

        if (job->VerboseLevel > 2)
        {
            fprintf(stderr,"\t\tCalculating row for coeff #%d\n",n);
        }

        for (i=MostLast; i>=LeastLast; i--)
        {
            if (PrevRow[i] == INFINITY) continue;
            for (q=job->MaxTable[unum][n]; q>=job->MinTable[unum][n]; q--)
            {

                k = job->Bpp[unum][n][q];
                if ((Loc = i+k) >= job->TableRowSize) continue;
                ErrTemp = PrevRow[i] + job->Err[n][q];
                if (ErrTemp < CurrRow[Loc])
                {

                    job->QChoice[unum][n][Loc] = ((Qentry) q);
                    CurrRow[Loc] = ErrTemp;

                }
            }
        }

        while (CurrRow[LeastLast]==INFINITY) LeastLast++;
        MostLast = job->TableRowSize-1;
        while (CurrRow[MostLast]==INFINITY) MostLast--;

        if (job->VerboseLevel > 3)
        {
            fprintf(stderr,"\t\t\tBPP range: %d/%d to %d/%d\n",
                    LeastLast,job->BppScale,MostLast,job->BppScale);
        }


        TempRow = CurrRow;
        CurrRow = PrevRow;
        PrevRow = TempRow;
        ResetRow(CurrRow,job->TableRowSize);
    }

    free(CurrRow);
    job->LeastLast[unum] = LeastLast;
    job->MostLast[unum] = MostLast;

    /*** ensure monotonicity ******/
    ErrTemp = INFINITY;
    for (i=LeastLast; i<=MostLast; i++)
    {
        if (PrevRow[i] >= ErrTemp) PrevRow[i] = INFINITY;
        else ErrTemp = PrevRow[i];
    }

    job->LastRow[unum] = PrevRow;

}


extern void CombineUnits(OptimJob *job)
{
    int i,j,lim, unum, Loc;
    FFLOAT ErrTemp, *PrevRow, *TempRow;

    /********* should not be called if job->NumTables==1 ******/
    /***** allocate ******/
    if ((job->CombinedLastRow = (FFLOAT *)
                                calloc(1,job->TableRowSize*sizeof(FFLOAT)))==NULL)
    {
        FatalError("CombineUnits out of memory");
    }
    ResetRow(job->CombinedLastRow,job->TableRowSize);

    if (job->NumTables > 2)
    {
        if ((PrevRow = (FFLOAT *)
                       calloc(1,job->TableRowSize*sizeof(FFLOAT)))==NULL)
        {
            FatalError("CombineUnits out of memory");
        }
        ResetRow(PrevRow, job->TableRowSize);
    }

    if ((job->ToCombine[1] = (int *)
                             calloc(1,job->TableRowSize*(job->NumTables-1)*sizeof(int)))==NULL)
    {
        FatalError("CombineUnits out of memory");
    }

    for (i=2; i<job->NumTables; i++)
        job->ToCombine[i] = job->ToCombine[i-1] + job->TableRowSize;

    if (job->VerboseLevel)
    {
        fprintf(stderr,"\tIncorporating unit 0\n");
        if (job->VerboseLevel > 1)
        {
            fprintf(stderr,"\t\tBPP range: %d/%d through %d/%d\n",
                    job->LeastLast[0],job->BppScale,job->MostLast[0],job->BppScale);
        }
        fprintf(stderr,"\tIncorporating unit 1\n");
    }


    /** combine rows **/
    for (i=job->LeastLast[0]; i<=job->MostLast[0]; i++)
    {
        lim = job->TableRowSize - i - 1;
        if (job->MostLast[1] < lim) lim = job->MostLast[1];
        for (j=job->LeastLast[1]; j<= lim; j++)
        {
            ErrTemp = job->LastRow[0][i] + job->LastRow[1][j];
            Loc = i+j;
            if (ErrTemp < job->CombinedLastRow[Loc])
            {
                job->CombinedLastRow[Loc] = ErrTemp;
                job->ToCombine[1][Loc] = j;
            }
        }
    }
    job->CombinedLeast=0;
    while (job->CombinedLastRow[job->CombinedLeast] == INFINITY)
        job->CombinedLeast++;
    job->CombinedMost=job->TableRowSize-1;
    while (job->CombinedLastRow[job->CombinedMost] == INFINITY)
        job->CombinedMost--;

    if (job->VerboseLevel > 1)
    {
        fprintf(stderr,"\t\tBPP range: %d/%d through %d/%d\n",
                job->CombinedLeast,job->BppScale,job->CombinedMost,job->BppScale);
    }

    for (unum=2; unum < job->NumTables; unum++)
    {
        if (job->VerboseLevel)
        {
            fprintf(stderr,"\tIncorporating unit %d\n",unum);
        }
        TempRow = job->CombinedLastRow;
        job->CombinedLastRow = PrevRow;
        PrevRow = TempRow;
        ResetRow(job->CombinedLastRow,job->TableRowSize);

        for (i=job->CombinedLeast; i<=job->CombinedMost; i++)
        {
            lim = job->TableRowSize - i - 1;
            if (job->MostLast[unum] < lim) lim = job->MostLast[unum];
            for (j=job->LeastLast[unum]; j<= lim; j++)
            {
                ErrTemp = PrevRow[i] + job->LastRow[unum][j];
                Loc = i+j;
                if (ErrTemp < job->CombinedLastRow[Loc])
                {
                    job->CombinedLastRow[Loc] = ErrTemp;
                    job->ToCombine[unum][Loc] = j;
                }
            }
        }
        while (job->CombinedLastRow[job->CombinedLeast] == INFINITY)
            job->CombinedLeast++;
        job->CombinedMost=job->TableRowSize-1;
        while (job->CombinedLastRow[job->CombinedMost] == INFINITY)
            job->CombinedMost--;

        if (job->VerboseLevel > 1)
        {
            fprintf(stderr,"\t\tBPP range: %d/%d through %d/%d\n",
                    job->CombinedLeast,job->BppScale,job->CombinedMost,job->BppScale);
        }
    }

    if (job->NumTables > 2) free(PrevRow);

    /*** ensure monotonicity ***/
    ErrTemp = INFINITY;
    for (i=job->CombinedLeast; i<=job->CombinedMost; i++)
    {
        if (job->CombinedLastRow[i] >= ErrTemp) job->CombinedLastRow[i] = INFINITY;
        else ErrTemp = job->CombinedLastRow[i];
    }

}

extern void Epilogue(OptimJob *job)
{
    int unum, i, Loc1, Loc2;
    FFLOAT err, psnr, bpp;
    FILE *fp;

    /**** compute LogSignalSq ****/
    job->LogTotalSignalSq = 0.0;
    for (unum=0; unum < job->NumTables; unum++)
    {
        job->LogTotalSignalSq += job->LogSignalSq[unum];
        job->LogSignalSq[unum] = ((FFLOAT) log10((double) job->LogSignalSq[unum]));
    }
    job->LogTotalSignalSq = ((FFLOAT) log10((double) job->LogTotalSignalSq));

    /**** take logarithms of all error rows ***/
    for (unum=0; unum < job->NumTables; unum++)
    {
        for (i=job->LeastLast[unum]; i<= job->MostLast[unum]; i++)
        {
            if (job->LastRow[unum][i] < INFINITY)
                job->LastRow[unum][i] =
                    ((FFLOAT) log10((double) job->LastRow[unum][i]));
        }
    }

    if (job->NumTables > 1)
    {
        for (i=job->CombinedLeast; i<= job->CombinedMost; i++)
        {
            if (job->CombinedLastRow[i] < INFINITY)
                job->CombinedLastRow[i] =
                    ((FFLOAT) log10((double) job->CombinedLastRow[i]));
        }
    }

    if (job->VerboseLevel)
    {
        /** print range **/
        if (job->NumTables > 1)
        {
            Loc1 = job->CombinedLeast;
            while ((Loc1 <= job->CombinedMost) &&
                    (job->CombinedLastRow[Loc1] >= INFINITY)) Loc1++;
            Loc2 = job->CombinedMost;
            while ((Loc2 >= job->CombinedLeast) &&
                    (job->CombinedLastRow[Loc2] >= INFINITY)) Loc2--;

            if ((Loc1 > job->CombinedMost) || (Loc2 < job->CombinedLeast))
            {
                fprintf(stderr,"No feasible quantization tables found in range\n");
            }
            else
            {
                err = job->CombinedLastRow[Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->BppScale);
                fprintf(stderr," Range: bpp (psnr) %lf (%lf) to ",
                        bpp, psnr);
                err = job->CombinedLastRow[Loc2];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc2)/((FFLOAT) job->BppScale);
                fprintf(stderr,"  %lf (%lf)\n", bpp, psnr);
            }
        }
        else
        {
            Loc1 = job->LeastLast[0];
            while ((Loc1 <= job->MostLast[0]) &&
                    (job->LastRow[0][Loc1] >= INFINITY)) Loc1++;
            Loc2 = job->MostLast[0];
            while ((Loc2 >= job->LeastLast[0]) &&
                    (job->LastRow[0][Loc2] >= INFINITY)) Loc2--;

            if ((Loc1 > job->MostLast[0]) || (Loc2 < job->LeastLast[0]))
            {
                fprintf(stderr,"No feasible quantization tables found in range\n");
            }
            else
            {
                err = job->LastRow[0][Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->BppScale);
                fprintf(stderr," Range: bpp (psnr) %lf (%lf) to ",
                        (double) bpp, (double) psnr);
                err = job->LastRow[0][Loc2];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc2)/((FFLOAT) job->BppScale);
                fprintf(stderr,"  %lf (%lf)\n", (double) bpp, (double) psnr);
            }
        }

    }

    if (job->DumpPlot)
    {
        if (!strcmp(job->PlotFileName,"-")) fp = stdout;
        else if ((fp = fopen(job->PlotFileName,"w")) == NULL)
        {
            fprintf(stderr,"Could not open plot file.. ignoring\n");
            return;
        }

        fprintf(fp,"# BPP -- PSNR plot generated by RDOPT\n");
        if (strcmp(job->PlotFileName,"-")) DumpJobChars(job,fp);


        if (job->NumTables > 1)
        {
            Loc2 = (job->CombinedMost - job->CombinedLeast)/PLOT_POINTS;
            Loc1 = job->CombinedLeast;
            for (i=0; i<PLOT_POINTS; i++)
            {
                while ((Loc1 <= job->CombinedMost) &&
                        (job->CombinedLastRow[Loc1] >= INFINITY)) Loc1++;
                if (Loc1 > job->CombinedMost) break;
                err = job->CombinedLastRow[Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->BppScale);
                fprintf(fp,"%lf   %lf\n", (double) bpp, (double) psnr);
                Loc1 += Loc2;
            }
            Loc1 = job->CombinedMost;
            while ((Loc1 >= job->CombinedLeast) &&
                    (job->CombinedLastRow[Loc1] >= INFINITY)) Loc1--;
            if (Loc1 >= job->CombinedLeast)
            {
                err = job->CombinedLastRow[Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->BppScale);
                fprintf(fp,"%lf   %lf\n", (double) bpp, (double) psnr);
            }
        }
        else
        {
            Loc2 = (job->MostLast[0] - job->LeastLast[0])/PLOT_POINTS;
            Loc1 = job->LeastLast[0];
            for (i=0; i<PLOT_POINTS; i++)
            {
                while ((Loc1 <= job->MostLast[0]) &&
                        (job->LastRow[0][Loc1] >= INFINITY)) Loc1++;
                if (Loc1 > job->MostLast[0]) break;
                err = job->LastRow[0][Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->BppScale);
                fprintf(fp,"%lf   %lf\n", (double) bpp, (double) psnr);
                Loc1 += Loc2;
            }
            Loc1 = job->MostLast[0];
            while ((Loc1 >= job->LeastLast[0]) &&
                    (job->LastRow[0][Loc1] >= INFINITY)) Loc1--;
            if (Loc1 >= job->LeastLast[0])
            {
                err = job->LastRow[0][Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->BppScale);
                fprintf(fp,"%lf   %lf\n", (double) bpp, (double) psnr);
            }
        }

        if (strcmp(job->PlotFileName,"-")) fclose(fp);
    }
}

