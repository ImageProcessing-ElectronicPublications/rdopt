
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

static void ResetRow(FFLOAT *row, int width)
{
    int i;
    for (i=0; i<width; i++) row[i] = INFINITY;
}

extern void Optimize(OptimJob *job, int unum)
{
    FFLOAT *CurrRow, *PrevRow, *TempRow;
    int LeastLast, MostLast, i, Loc, Loc2, n, q, k;
    int t2, indx, tspanplus1 = job->ThreshSpan+1;
    FFLOAT ErrTemp;


    /****** allocate ******/
    if ((CurrRow = (FFLOAT *)
                   calloc(1,job->opt_method.dp.TableRowSize*sizeof(FFLOAT)))==NULL)
    {
        FatalError("Optimize out of memory");
    }
    if ((PrevRow = (FFLOAT *)
                   calloc(1,job->opt_method.dp.TableRowSize*sizeof(FFLOAT)))==NULL)
    {
        FatalError("Optimize out of memory");
    }

    ResetRow(CurrRow,job->opt_method.dp.TableRowSize);
    ResetRow(PrevRow,job->opt_method.dp.TableRowSize);

    if ((job->opt_method.dp.QChoice[unum][0] = (Qentry *)
            calloc(1,job->opt_method.dp.TableRowSize*64*sizeof(Qentry)))==NULL)
    {
        FatalError("Optimize out of memory");
    }

    for (n=1; n<64; n++)
    {
        job->opt_method.dp.QChoice[unum][n] = job->opt_method.dp.QChoice[unum][n-1] +
                                              job->opt_method.dp.TableRowSize;
    }

    if (job->ThreshSpan > 0)
    {

        if ((job->opt_method.dp.TChoice[unum][0] = (Qentry *)
                calloc(1,job->opt_method.dp.TableRowSize*64*sizeof(Qentry)))==NULL)
        {
            FatalError("Optimize out of memory");
        }

        for (n=1; n<64; n++)
        {
            job->opt_method.dp.TChoice[unum][n] = job->opt_method.dp.TChoice[unum][n-1] +
                                                  job->opt_method.dp.TableRowSize;
        }

        if (job->opt_method.dp.ErrEncodesRow[unum][63])
        {
            if (job->VerboseLevel > 2)
            {
                fprintf(errfile,"\t\tCalculating row for compacted coeff #%d\n",63);
                fflush(errfile);
            }
            Loc2 = 0;
            for (Loc = job->opt_method.dp.BppOffsetInEncoding[unum][63];
                    Loc <= job->opt_method.dp.BppMaxInEncoding[unum][63]; Loc++,Loc2++)
            {
                PrevRow[Loc] = job->opt_method.dp.Err[unum][63][Loc2];
                job->opt_method.dp.QChoice[unum][63][Loc] =
                    job->opt_method.dp.QforBpp[63][Loc2];
                job->opt_method.dp.TChoice[unum][63][Loc] =
                    job->opt_method.dp.TforBpp[63][Loc2];
            }
            free(job->opt_method.dp.QforBpp[63]);
            free(job->opt_method.dp.TforBpp[63]);
            LeastLast = job->opt_method.dp.BppOffsetInEncoding[unum][63];
            MostLast = job->opt_method.dp.BppMaxInEncoding[unum][63];
        }
        else
        {
            if (job->VerboseLevel > 2)
            {
                fprintf(errfile,"\t\tCalculating row for coeff #%d\n",63);
                fflush(errfile);
            }
            LeastLast = job->opt_method.dp.TableRowSize-1;
            MostLast = 0;
            indx = job->MinTable[unum][63]*tspanplus1;
            for (q=job->MinTable[unum][63]; q<=job->MaxTable[unum][63]; q++)
            {
                for (t2=0; t2<=job->ThreshSpan; t2++,indx++)
                {
                    if ((Loc = job->opt_method.dp.Bpp[unum][63][indx])
                            >= job->opt_method.dp.TableRowSize)
                    {
                        continue;
                    }
                    if (job->opt_method.dp.Err[unum][63][indx] < PrevRow[Loc])
                    {
                        if (Loc < LeastLast) LeastLast = Loc;
                        if (Loc > MostLast) MostLast = Loc;
                        job->opt_method.dp.QChoice[unum][63][Loc] = ((Qentry) q);
                        job->opt_method.dp.TChoice[unum][63][Loc] = ((Qentry) t2);
                        PrevRow[Loc] = job->opt_method.dp.Err[unum][63][indx];
                    }
                }
            }
        }

        if (job->VerboseLevel > 3)
        {
            fprintf(errfile,"\t\t\tBPP range: %d/%d to %d/%d\n",
                    LeastLast,job->opt_method.dp.BppScale,MostLast,job->opt_method.dp.BppScale);
            fflush(errfile);
        }

        for (n=62; n>=0; n--)
        {
            if (job->opt_method.dp.ErrEncodesRow[unum][n])
            {

                if (job->VerboseLevel > 2)
                {
                    fprintf(errfile,"\t\tCalculating row for compacted coeff #%d\n",n);
                    fflush(errfile);
                }


                for (i=MostLast; i>=LeastLast; i--)
                {
                    if (PrevRow[i] == INFINITY) continue;
                    Loc2 = 0;
                    for (k = job->opt_method.dp.BppOffsetInEncoding[unum][n];
                            k <= job->opt_method.dp.BppMaxInEncoding[unum][n]; k++,Loc2++)
                    {
                        if ((Loc = i+k) >= job->opt_method.dp.TableRowSize) continue;
                        ErrTemp = PrevRow[i] + job->opt_method.dp.Err[unum][n][Loc2];
                        if (ErrTemp < CurrRow[Loc])
                        {
                            CurrRow[Loc] = ErrTemp;
                            job->opt_method.dp.QChoice[unum][n][Loc] =
                                job->opt_method.dp.QforBpp[n][Loc2];
                            job->opt_method.dp.TChoice[unum][n][Loc] =
                                job->opt_method.dp.TforBpp[n][Loc2];
                        }
                    }
                }
                free(job->opt_method.dp.QforBpp[n]);
                free(job->opt_method.dp.TforBpp[n]);
            }
            else
            {

                if (job->VerboseLevel > 2)
                {
                    fprintf(errfile,"\t\tCalculating row for coeff #%d\n",n);
                    fflush(errfile);
                }

                for (i=MostLast; i>=LeastLast; i--)
                {
                    if (PrevRow[i] == INFINITY) continue;
                    indx = job->MinTable[unum][n]*tspanplus1;
                    for (q=job->MinTable[unum][n]; q<=job->MaxTable[unum][n]; q++)
                    {
                        for (t2=0; t2<=job->ThreshSpan; t2++,indx++)
                        {

                            k = job->opt_method.dp.Bpp[unum][n][indx];
                            if ((Loc = i+k) >= job->opt_method.dp.TableRowSize) continue;
                            ErrTemp = PrevRow[i] + job->opt_method.dp.Err[unum][n][indx];
                            if (ErrTemp < CurrRow[Loc])
                            {

                                CurrRow[Loc] = ErrTemp;
                                job->opt_method.dp.QChoice[unum][n][Loc] = ((Qentry) q);
                                job->opt_method.dp.TChoice[unum][n][Loc] = ((Qentry) t2);

                            }
                        }
                    }
                }
            }

            while (CurrRow[LeastLast]==INFINITY) LeastLast++;
            MostLast = job->opt_method.dp.TableRowSize-1;
            while (CurrRow[MostLast]==INFINITY) MostLast--;

            if (job->VerboseLevel > 3)
            {
                fprintf(errfile,"\t\t\tBPP range: %d/%d to %d/%d\n",
                        LeastLast,job->opt_method.dp.BppScale,MostLast,job->opt_method.dp.BppScale);
                fflush(errfile);
            }


            TempRow = CurrRow;
            CurrRow = PrevRow;
            PrevRow = TempRow;
            ResetRow(CurrRow,job->opt_method.dp.TableRowSize);
        }
    }
    else
    {
        /* No Thresholding */
        if (job->opt_method.dp.ErrEncodesRow[unum][63])
        {

            if (job->VerboseLevel > 2)
            {
                fprintf(errfile,"\t\tCalculating row for compacted coeff #%d\n",63);
                fflush(errfile);
            }

            Loc2 = 0;
            for (Loc = job->opt_method.dp.BppOffsetInEncoding[unum][63];
                    Loc <= job->opt_method.dp.BppMaxInEncoding[unum][63]; Loc++,Loc2++)
            {
                PrevRow[Loc] = job->opt_method.dp.Err[unum][63][Loc2];
                job->opt_method.dp.QChoice[unum][63][Loc] =
                    job->opt_method.dp.QforBpp[63][Loc2];
            }
            free(job->opt_method.dp.QforBpp[63]);
            LeastLast = job->opt_method.dp.BppOffsetInEncoding[unum][63];
            MostLast = job->opt_method.dp.BppMaxInEncoding[unum][63];
        }
        else
        {

            if (job->VerboseLevel > 2)
            {
                fprintf(errfile,"\t\tCalculating row for coeff #%d\n",63);
                fflush(errfile);
            }

            LeastLast = job->opt_method.dp.TableRowSize-1;
            MostLast = 0;
            for (q=job->MinTable[unum][63]; q<=job->MaxTable[unum][63]; q++)
            {
                if ((Loc = job->opt_method.dp.Bpp[unum][63][q])
                        >= job->opt_method.dp.TableRowSize)
                {
                    continue;
                }
                if (job->opt_method.dp.Err[unum][63][q] < PrevRow[Loc])
                {
                    if (Loc < LeastLast) LeastLast = Loc;
                    if (Loc > MostLast) MostLast = Loc;
                    job->opt_method.dp.QChoice[unum][63][Loc] = ((Qentry) q);
                    PrevRow[Loc] = job->opt_method.dp.Err[unum][63][q];
                }
            }
        }

        if (job->VerboseLevel > 3)
        {
            fprintf(errfile,"\t\t\tBPP range: %d/%d to %d/%d\n",
                    LeastLast,job->opt_method.dp.BppScale,MostLast,job->opt_method.dp.BppScale);
            fflush(errfile);
        }

        for (n=62; n>=0; n--)
        {
            if (job->opt_method.dp.ErrEncodesRow[unum][n])
            {

                if (job->VerboseLevel > 2)
                {
                    fprintf(errfile,"\t\tCalculating row for compacted coeff #%d\n",n);
                    fflush(errfile);
                }

                for (i=MostLast; i>=LeastLast; i--)
                {
                    if (PrevRow[i] == INFINITY) continue;
                    Loc2 = 0;
                    for (k = job->opt_method.dp.BppOffsetInEncoding[unum][n];
                            k <= job->opt_method.dp.BppMaxInEncoding[unum][n]; k++,Loc2++)
                    {
                        if ((Loc = i+k) >= job->opt_method.dp.TableRowSize) continue;
                        ErrTemp = PrevRow[i] + job->opt_method.dp.Err[unum][n][Loc2];
                        if (ErrTemp < CurrRow[Loc])
                        {
                            CurrRow[Loc] = ErrTemp;
                            job->opt_method.dp.QChoice[unum][n][Loc] =
                                job->opt_method.dp.QforBpp[n][Loc2];
                        }
                    }
                }
                free(job->opt_method.dp.QforBpp[n]);
            }
            else
            {

                if (job->VerboseLevel > 2)
                {
                    fprintf(errfile,"\t\tCalculating row for coeff #%d\n",n);
                    fflush(errfile);
                }

                for (i=MostLast; i>=LeastLast; i--)
                {
                    if (PrevRow[i] == INFINITY) continue;
                    for (q=job->MinTable[unum][n]; q<=job->MaxTable[unum][n]; q++)
                    {

                        k = job->opt_method.dp.Bpp[unum][n][q];
                        if ((Loc = i+k) >= job->opt_method.dp.TableRowSize) continue;
                        ErrTemp = PrevRow[i] + job->opt_method.dp.Err[unum][n][q];
                        if (ErrTemp < CurrRow[Loc])
                        {

                            CurrRow[Loc] = ErrTemp;
                            job->opt_method.dp.QChoice[unum][n][Loc] = ((Qentry) q);
                        }
                    }
                }
            }

            while (CurrRow[LeastLast]==INFINITY) LeastLast++;
            MostLast = job->opt_method.dp.TableRowSize-1;
            while (CurrRow[MostLast]==INFINITY) MostLast--;

            if (job->VerboseLevel > 3)
            {
                fprintf(errfile,"\t\t\tBPP range: %d/%d to %d/%d\n",
                        LeastLast,job->opt_method.dp.BppScale,MostLast,job->opt_method.dp.BppScale);
                fflush(errfile);
            }


            TempRow = CurrRow;
            CurrRow = PrevRow;
            PrevRow = TempRow;
            ResetRow(CurrRow,job->opt_method.dp.TableRowSize);
        }
    }

    free(CurrRow);
    job->opt_method.dp.LeastLast[unum] = LeastLast;
    job->opt_method.dp.MostLast[unum] = MostLast;

    /*** ensure monotonicity ******/
    ErrTemp = INFINITY;
    for (i=LeastLast; i<=MostLast; i++)
    {
        if (PrevRow[i] >= ErrTemp) PrevRow[i] = INFINITY;
        else ErrTemp = PrevRow[i];
    }

    job->opt_method.dp.LastRow[unum] = PrevRow;

}


extern void CombineUnits(OptimJob *job)
{
    int i,j,lim, unum, Loc;
    FFLOAT ErrTemp, *PrevRow, *TempRow;

    /********* should not be called if job->NumTables==1 ******/
    /***** allocate ******/
    if ((job->opt_method.dp.CombinedLastRow = (FFLOAT *)
            calloc(1,job->opt_method.dp.TableRowSize*sizeof(FFLOAT)))==NULL)
    {
        FatalError("CombineUnits out of memory");
    }
    ResetRow(job->opt_method.dp.CombinedLastRow,job->opt_method.dp.TableRowSize);

    if (job->NumTables > 2)
    {
        if ((PrevRow = (FFLOAT *)
                       calloc(1,job->opt_method.dp.TableRowSize*sizeof(FFLOAT)))==NULL)
        {
            FatalError("CombineUnits out of memory");
        }
        ResetRow(PrevRow, job->opt_method.dp.TableRowSize);
    }

    if ((job->opt_method.dp.ToCombine[1] = (int *)
                                           calloc(1,job->opt_method.dp.TableRowSize*(job->NumTables-1)*sizeof(int)))==NULL)
    {
        FatalError("CombineUnits out of memory");
    }

    for (i=2; i<job->NumTables; i++)
        job->opt_method.dp.ToCombine[i] = job->opt_method.dp.ToCombine[i-1] + job->opt_method.dp.TableRowSize;

    if (job->VerboseLevel)
    {
        fprintf(errfile,"\tIncorporating unit 0\n");
        fflush(errfile);
        if (job->VerboseLevel > 1)
        {
            fprintf(errfile,"\t\tBPP range: %d/%d through %d/%d\n",
                    job->opt_method.dp.LeastLast[0],job->opt_method.dp.BppScale,job->opt_method.dp.MostLast[0],job->opt_method.dp.BppScale);
            fflush(errfile);
        }
        fprintf(errfile,"\tIncorporating unit 1\n");
        fflush(errfile);
    }


    /** combine rows **/
    for (i=job->opt_method.dp.LeastLast[0]; i<=job->opt_method.dp.MostLast[0]; i++)
    {
        lim = job->opt_method.dp.TableRowSize - i - 1;
        if (job->opt_method.dp.MostLast[1] < lim) lim = job->opt_method.dp.MostLast[1];
        for (j=job->opt_method.dp.LeastLast[1]; j<= lim; j++)
        {
            ErrTemp = job->opt_method.dp.LastRow[0][i] + job->opt_method.dp.LastRow[1][j];
            Loc = i+j;
            if (ErrTemp < job->opt_method.dp.CombinedLastRow[Loc])
            {
                job->opt_method.dp.CombinedLastRow[Loc] = ErrTemp;
                job->opt_method.dp.ToCombine[1][Loc] = j;
            }
        }
    }
    job->opt_method.dp.CombinedLeast=0;
    while (job->opt_method.dp.CombinedLastRow[job->opt_method.dp.CombinedLeast] == INFINITY)
        job->opt_method.dp.CombinedLeast++;
    job->opt_method.dp.CombinedMost=job->opt_method.dp.TableRowSize-1;
    while (job->opt_method.dp.CombinedLastRow[job->opt_method.dp.CombinedMost] == INFINITY)
        job->opt_method.dp.CombinedMost--;

    if (job->VerboseLevel > 1)
    {
        fprintf(errfile,"\t\tBPP range: %d/%d through %d/%d\n",
                job->opt_method.dp.CombinedLeast,job->opt_method.dp.BppScale,job->opt_method.dp.CombinedMost,job->opt_method.dp.BppScale);
        fflush(errfile);
    }

    for (unum=2; unum < job->NumTables; unum++)
    {
        if (job->VerboseLevel)
        {
            fprintf(errfile,"\tIncorporating unit %d\n",unum);
            fflush(errfile);
        }
        TempRow = job->opt_method.dp.CombinedLastRow;
        job->opt_method.dp.CombinedLastRow = PrevRow;
        PrevRow = TempRow;
        ResetRow(job->opt_method.dp.CombinedLastRow,job->opt_method.dp.TableRowSize);

        for (i=job->opt_method.dp.CombinedLeast; i<=job->opt_method.dp.CombinedMost; i++)
        {
            lim = job->opt_method.dp.TableRowSize - i - 1;
            if (job->opt_method.dp.MostLast[unum] < lim) lim = job->opt_method.dp.MostLast[unum];
            for (j=job->opt_method.dp.LeastLast[unum]; j<= lim; j++)
            {
                ErrTemp = PrevRow[i] + job->opt_method.dp.LastRow[unum][j];
                Loc = i+j;
                if (ErrTemp < job->opt_method.dp.CombinedLastRow[Loc])
                {
                    job->opt_method.dp.CombinedLastRow[Loc] = ErrTemp;
                    job->opt_method.dp.ToCombine[unum][Loc] = j;
                }
            }
        }
        while (job->opt_method.dp.CombinedLastRow[job->opt_method.dp.CombinedLeast] == INFINITY)
            job->opt_method.dp.CombinedLeast++;
        job->opt_method.dp.CombinedMost=job->opt_method.dp.TableRowSize-1;
        while (job->opt_method.dp.CombinedLastRow[job->opt_method.dp.CombinedMost] == INFINITY)
            job->opt_method.dp.CombinedMost--;

        if (job->VerboseLevel > 1)
        {
            fprintf(errfile,"\t\tBPP range: %d/%d through %d/%d\n",
                    job->opt_method.dp.CombinedLeast,job->opt_method.dp.BppScale,job->opt_method.dp.CombinedMost,job->opt_method.dp.BppScale);
            fflush(errfile);
        }
    }

    if (job->NumTables > 2) free(PrevRow);

    /*** ensure monotonicity ***/
    ErrTemp = INFINITY;
    for (i=job->opt_method.dp.CombinedLeast; i<=job->opt_method.dp.CombinedMost; i++)
    {
        if (job->opt_method.dp.CombinedLastRow[i] >= ErrTemp) job->opt_method.dp.CombinedLastRow[i] = INFINITY;
        else ErrTemp = job->opt_method.dp.CombinedLastRow[i];
    }

}

extern void Epilogue(OptimJob *job)
{
    int unum, i, Loc1, Loc2, plotUpper;
    FFLOAT err, psnr, bpp, actualBpp;
    FILE *fp;
    int StartingPt[MAXCOMPONENTS];
    int Q[64], T[64], BppDist[64];
    FFLOAT ErrDist[64];
    char qfname[STRLENMAX];

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
        for (i=job->opt_method.dp.LeastLast[unum]; i<= job->opt_method.dp.MostLast[unum]; i++)
        {
            if (job->opt_method.dp.LastRow[unum][i] < INFINITY)
                job->opt_method.dp.LastRow[unum][i] =
                    ((FFLOAT) log10((double)
                                    (job->opt_method.dp.LastRow[unum][i]/job->CompWeights[unum])));
        }
    }

    if (job->NumTables > 1)
    {
        for (i=job->opt_method.dp.CombinedLeast; i<= job->opt_method.dp.CombinedMost; i++)
        {
            if (job->opt_method.dp.CombinedLastRow[i] < INFINITY)
                job->opt_method.dp.CombinedLastRow[i] =
                    ((FFLOAT) log10((double) job->opt_method.dp.CombinedLastRow[i]));
        }
    }

    if (job->useCorrection)
    {
        if (job->TheImage.ImFileName[0] == '\0')
        {
            fprintf(errfile,"Ignoring -correct: image read off stdin\n");
            job->useCorrection = FALSE;
            goto afterCorrection;
        }
        bpp = job->correctionBpp*((FFLOAT) job->opt_method.dp.BppScale);
        Loc1 = RoundOff(bpp);
        if (job->NumTables==1)
        {
            Loc2 = ScanRowForBpp(job->opt_method.dp.LastRow[0], Loc1,
                                 job->opt_method.dp.LeastLast[0],
                                 job->opt_method.dp.MostLast[0]);
            if (Loc2 < 0)
            {
                /* should not happen */
                fprintf(errfile,"Error applying correction.. ignoring -correct\n");
                job->useCorrection = FALSE;
                goto afterCorrection;
            }
        }
        else
        {
            Loc2 = ScanRowForBpp(job->opt_method.dp.CombinedLastRow, Loc1,
                                 job->opt_method.dp.CombinedLeast,
                                 job->opt_method.dp.CombinedMost);
            if (Loc2 < 0)
            {
                /* should not happen */
                fprintf(errfile,"Error applying correction.. ignoring -correct\n");
                job->useCorrection = FALSE;
                goto afterCorrection;
            }
            else
            {
                GetStartingPoints(job->opt_method.dp.ToCombine, StartingPt,
                                  Loc2, job->NumTables);
            }
        }
        /* continue applying correction.. */
        sprintf(qfname,"/tmp/rdopt.%s.correctionQ\0",
                FileNameTail(job->TheImage.ImFileName));
        if ((fp = fopen(qfname,"w")) == NULL)
        {
            fprintf(errfile,"Could not open %s\n",qfname);
            job->useCorrection = FALSE;
            goto afterCorrection;
        }

        if (job->NumTables==1)
        {
            if (job->ThreshSpan > 0)
            {
                RecoverQandT(Loc2, 0, job, Q, T, BppDist, ErrDist);
                WriteIntTable(fp, Q, "");
                WriteThreshTable(fp, T, "#T ");
            }
            else
            {
                RecoverQ(Loc2, 0, job, Q, BppDist, ErrDist);
                WriteIntTable(fp, Q, "");
            }
        }
        else
        {
            for (unum=0; unum < job->NumTables; unum++)
            {
                if (job->ThreshSpan > 0)
                {
                    RecoverQandT(StartingPt[unum], unum, job, Q, T, BppDist, ErrDist);
                    WriteIntTable(fp, Q, "");
                    WriteThreshTable(fp, T, "#T ");
                }
                else
                {
                    RecoverQ(StartingPt[unum], unum, job, Q, BppDist, ErrDist);
                    WriteIntTable(fp, Q, "");
                }
            }
        }
        fclose(fp);

        if (GetActualBpp(job, qfname, &actualBpp))
        {
            job->correctionBpp = ((FFLOAT) Loc2)/((FFLOAT) job->opt_method.dp.BppScale);
            job->addToTableBpp = actualBpp - job->correctionBpp;
        }
        else
        {
            job->useCorrection = FALSE;
        }
        remove(qfname);
    }


afterCorrection:
    if (!job->Silent)
    {
        /** print range **/
        if (job->NumTables > 1)
        {
            Loc1 = job->opt_method.dp.CombinedLeast;
            while ((Loc1 <= job->opt_method.dp.CombinedMost) &&
                    (job->opt_method.dp.CombinedLastRow[Loc1] >= INFINITY)) Loc1++;
            Loc2 = job->opt_method.dp.CombinedMost;
            while ((Loc2 >= job->opt_method.dp.CombinedLeast) &&
                    (job->opt_method.dp.CombinedLastRow[Loc2] >= INFINITY)) Loc2--;

            if ((Loc1 > job->opt_method.dp.CombinedMost) || (Loc2 < job->opt_method.dp.CombinedLeast))
            {
                fprintf(errfile,"No feasible quantization tables found in range\n");
                fflush(errfile);
            }
            else
            {
                err = job->opt_method.dp.CombinedLastRow[Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->opt_method.dp.BppScale);
                if (job->useCorrection) bpp += job->addToTableBpp;
                fprintf(errfile," Range: bpp (psnr) %lf (%lf) to ",
                        bpp, psnr);
                err = job->opt_method.dp.CombinedLastRow[Loc2];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc2)/((FFLOAT) job->opt_method.dp.BppScale);
                if (job->useCorrection) bpp += job->addToTableBpp;
                fprintf(errfile,"  %lf (%lf)\n", bpp, psnr);
                fflush(errfile);
            }
        }
        else
        {
            Loc1 = job->opt_method.dp.LeastLast[0];
            while ((Loc1 <= job->opt_method.dp.MostLast[0]) &&
                    (job->opt_method.dp.LastRow[0][Loc1] >= INFINITY)) Loc1++;
            Loc2 = job->opt_method.dp.MostLast[0];
            while ((Loc2 >= job->opt_method.dp.LeastLast[0]) &&
                    (job->opt_method.dp.LastRow[0][Loc2] >= INFINITY)) Loc2--;

            if ((Loc1 > job->opt_method.dp.MostLast[0]) || (Loc2 < job->opt_method.dp.LeastLast[0]))
            {
                fprintf(errfile,"No feasible quantization tables found in range\n");
                fflush(errfile);
            }
            else
            {
                err = job->opt_method.dp.LastRow[0][Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->opt_method.dp.BppScale);
                if (job->useCorrection) bpp += job->addToTableBpp;
                fprintf(errfile," Range: bpp (psnr) %lf (%lf) to ",
                        (double) bpp, (double) psnr);
                err = job->opt_method.dp.LastRow[0][Loc2];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc2)/((FFLOAT) job->opt_method.dp.BppScale);
                if (job->useCorrection) bpp += job->addToTableBpp;
                fprintf(errfile,"  %lf (%lf)\n", (double) bpp, (double) psnr);
                fflush(errfile);
            }
        }

    }

    if (job->DumpPlot)
    {
        plotUpper = ((int) ((FFLOAT)
                            ((FFLOAT) job->opt_method.dp.BppScale)*job->PlotBppMax+0.5));
        if (!strcmp(job->PlotFileName,"-")) fp = stdout;
        else if ((fp = fopen(job->PlotFileName,"w")) == NULL)
        {
            fprintf(errfile,"Could not open plot file.. ignoring\n");
            return;
            fflush(errfile);
        }

        fprintf(fp,"# BPP -- PSNR plot generated by RDOPT\n");
        if (strcmp(job->PlotFileName,"-")) DumpJobChars(job,fp);


        if (job->NumTables > 1)
        {
            if (job->opt_method.dp.CombinedMost < plotUpper)
                plotUpper = job->opt_method.dp.CombinedMost;
            Loc2 = (plotUpper - job->opt_method.dp.CombinedLeast)/job->PlotPoints;
            Loc1 = job->opt_method.dp.CombinedLeast;
            for (i=0; i<job->PlotPoints; i++)
            {
                while ((Loc1 <= job->opt_method.dp.CombinedMost) &&
                        (job->opt_method.dp.CombinedLastRow[Loc1] >= INFINITY)) Loc1++;
                if (Loc1 > job->opt_method.dp.CombinedMost) break;
                err = job->opt_method.dp.CombinedLastRow[Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->opt_method.dp.BppScale);
                if (job->useCorrection) bpp += job->addToTableBpp;
                fprintf(fp,"%lf   %lf\n", (double) bpp, (double) psnr);
                Loc1 += Loc2;
            }
            Loc1 = plotUpper;
            while ((Loc1 >= job->opt_method.dp.CombinedLeast) &&
                    (job->opt_method.dp.CombinedLastRow[Loc1] >= INFINITY)) Loc1--;
            if (Loc1 >= job->opt_method.dp.CombinedLeast)
            {
                err = job->opt_method.dp.CombinedLastRow[Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->opt_method.dp.BppScale);
                if (job->useCorrection) bpp += job->addToTableBpp;
                fprintf(fp,"%lf   %lf\n", (double) bpp, (double) psnr);
            }
        }
        else
        {
            if (job->opt_method.dp.MostLast[0] < plotUpper)
                plotUpper = job->opt_method.dp.MostLast[0];
            Loc2 = (plotUpper - job->opt_method.dp.LeastLast[0])/job->PlotPoints;
            Loc1 = job->opt_method.dp.LeastLast[0];
            for (i=0; i<job->PlotPoints; i++)
            {
                while ((Loc1 <= job->opt_method.dp.MostLast[0]) &&
                        (job->opt_method.dp.LastRow[0][Loc1] >= INFINITY)) Loc1++;
                if (Loc1 > job->opt_method.dp.MostLast[0]) break;
                err = job->opt_method.dp.LastRow[0][Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->opt_method.dp.BppScale);
                if (job->useCorrection) bpp += job->addToTableBpp;
                fprintf(fp,"%lf   %lf\n", (double) bpp, (double) psnr);
                Loc1 += Loc2;
            }
            Loc1 = plotUpper;
            while ((Loc1 >= job->opt_method.dp.LeastLast[0]) &&
                    (job->opt_method.dp.LastRow[0][Loc1] >= INFINITY)) Loc1--;
            if (Loc1 >= job->opt_method.dp.LeastLast[0])
            {
                err = job->opt_method.dp.LastRow[0][Loc1];
                psnr = (job->LogTotalPeakSq - err)*10.0;
                bpp = ((FFLOAT) Loc1)/((FFLOAT) job->opt_method.dp.BppScale);
                if (job->useCorrection) bpp += job->addToTableBpp;
                fprintf(fp,"%lf   %lf\n", (double) bpp, (double) psnr);
            }
        }

        if (strcmp(job->PlotFileName,"-")) fclose(fp);
    }
}

