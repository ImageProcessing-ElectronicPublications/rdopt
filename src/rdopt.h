
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
#include <math.h>

#include "version.h"
#include "precision.h"

#define CQTAB "RDOPT.Q" VERSION
#define CQHIST "RDOPT.H" VERSION

#include "Image.h"

extern FILE *errfile;

/* quantization table entries in the range 1..QTABENTRYMAX */
#if (QTABBITS==8)
#define QTABENTRYMAX 255
typedef unsigned char Qentry;
#elif (QTABBITS==10)
#define QTABENTRYMAX 1023
typedef unsigned short Qentry;
#elif (QTABBITS==12)
#define QTABENTRYMAX 4095
typedef unsigned short Qentry;
#elif (QTABBITS==16)
#define QTABENTRYMAX 65535
typedef unsigned short Qentry;
#else
/** Force syntax error **/
QTABBITS must be 8,12,or16
#endif


#ifndef STRLENMAX
#define STRLENMAX 200
#endif

#ifndef boolean
#define boolean int
#endif
#ifndef TRUE
#define TRUE 1
#define FALSE 0
#endif

#define FFLOAT double

#define INFINITY ((FFLOAT) 1e300) /* machine-dependent */

#define ZERO_CORRECTION ((FFLOAT) ((MAXSAMPLE+1)<<2))

#define PLOT_POINTS 20 /* default # of points to plot with DumpPlot set */

#define NUM_HOOKS 20 /* default # of lamda intervals in
			lagrangian opt */
#include "opers.h"

/*********************************
type definitions
**********************************/

/******* histograms **********/
struct HistStruct
{
    /*
      PlusCount[v] will store the number of times the DCT
      coefficient (for which this struct is storing the
      counts) takes a postive value that gets Discretized
      to the integer v. (For definition of "Discretization"
      see opers.h)
      MinusCount[v] will store the number of times the DCT
      coefficient takes a negative value that gets Discretized
      to the integer -v.
      The array bounds will be 0..PlusSize and 0..MinusSize,
      respectively.
    */

    int PlusSize;
    int MinusSize;
    int *PlusCount;
    int *MinusCount;
};
typedef struct HistStruct Hist;



/************** OptimizationJob ***********/
struct OptimJobStruct
{
    Image TheImage;
    int NumTables; /* Number of quantization tables to be found.
		    If NumTables < TheImage.NumComponents, then
		    the last table will be used for all the
		    color components numbered (NumTables - 1)
		    thru (TheImage.NumComponents - 1). We call
		    each set of components being quantized together
		    a "Unit"
                 */

    boolean MapQ;
    boolean Silent;
    int VerboseLevel;
    boolean ReadCmdFromFile;
    char CmdFileName[STRLENMAX];

    int ThreshSpan;
    /* for qentry=q, thresholds upto
     q/2 + ThreshSpan/2 will be tried */

    boolean UseLagrangian;

    boolean useCorrection;
    FFLOAT correctionBpp;
    FFLOAT addToTableBpp;

    union
    {
        struct
        {
            /*** components to be used for dynamic programming  ***/

            FFLOAT BppMax;
            int BppScale;
            int TableRowSize;
            FFLOAT *LastRow[MAXCOMPONENTS]; /* last row of Dyn Prog table
  				     for each unit */
            int LeastLast[MAXCOMPONENTS]; /* lower and */
            int MostLast[MAXCOMPONENTS];  /* upper bounds on the extent of LastRow */
            FFLOAT *CombinedLastRow;
            int CombinedLeast;
            int CombinedMost;
            int *ToCombine[MAXCOMPONENTS]; /* ToCombine[i][bits] will give
  				    the starting point in LastRow[i] */

            Qentry *QChoice[MAXCOMPONENTS][64];
            Qentry *TChoice[MAXCOMPONENTS][64];

            FFLOAT *Err[MAXCOMPONENTS][64];
            int *Bpp[MAXCOMPONENTS][64]; /* BPP for each unit needs
  				     to be retained */
            /* explain the complexity of encoded use of Err and Bpp some day */
            boolean ErrEncodesRow[MAXCOMPONENTS][64];
            Qentry *QforBpp[64];
            Qentry *TforBpp[64];

            int BppOffsetInEncoding[MAXCOMPONENTS][64];
            int BppMaxInEncoding[MAXCOMPONENTS][64];


        }  dp;

        struct
        {
            /***** components for Lagrangian method  ***/
            FFLOAT ErrScale; /* each error will be divided by this */
            FFLOAT *Err[MAXCOMPONENTS][64];
            FFLOAT *Bpp[MAXCOMPONENTS][64];
            int    *SortedIndex[MAXCOMPONENTS][64];
            int     IndexEntries[MAXCOMPONENTS][64];
            FFLOAT *Lambda[MAXCOMPONENTS][64];
            FFLOAT LambdaMax;
            FFLOAT BppMax;
            FFLOAT BppMin;
            FFLOAT ErrMax;
            FFLOAT ErrMin;

            FFLOAT *BppHook;
            FFLOAT *LambdaHook;
            FFLOAT *ErrHook;
            int NumHooks;

        }  lagr;


    } opt_method;

    boolean WeightedComps;
    FFLOAT CompWeights[MAXCOMPONENTS];

    boolean WeightedCoefs[MAXCOMPONENTS];
    FFLOAT CoefWeights[MAXCOMPONENTS][64];

    FFLOAT bppWeight[MAXCOMPONENTS]; /* bpp = \sum_{i=0..NumTables-1}
					       bpp-for-ith-unit *
					       bppWeight[i] */
    FFLOAT RemainingSignal[MAXCOMPONENTS][64];
    FFLOAT LogSignalSq[MAXCOMPONENTS];
    FFLOAT LogTotalSignalSq;
    FFLOAT LogPeakSq[MAXCOMPONENTS];
    FFLOAT LogTotalPeakSq;
    FFLOAT NumPixels[MAXCOMPONENTS];
    FFLOAT NumBlocks[MAXCOMPONENTS];
    FFLOAT TotalNumPixels;
    FFLOAT NumPixelsForBpp;

    boolean ImFilePresent;
    boolean HistFilePresent;
    char HistFile[STRLENMAX];

    Hist Histogram[64];
    boolean DumpStats;
    boolean OnlyDumpStats;

    int MinTable[MAXCOMPONENTS][64];
    int MaxTable[MAXCOMPONENTS][64];
    int DCclamp;

    boolean UseDCDPCM;

    boolean DumpPlot;
    int PlotPoints;
    FFLOAT PlotBppMax;
    char PlotFileName[STRLENMAX];

    int bppplane; /* -1 implies all planes */


};

typedef struct OptimJobStruct OptimJob;


#define HEADER_BYTES 200
/** approximate number of header bytes needed for each unit **/



/********************** prototypes ************************/



extern void SetHistogram(Hist *H, Image *Im, int cnum, FFLOAT *mssq,
                         FFLOAT *SigTab, boolean UseDCDPCM);

extern void HistIncrCount(Hist *H, int coefnum, int val,
                          boolean waspositive);

extern void InitHistogram(Hist *H);
extern void ReadHistogram(FILE *hfile, Hist *H, FFLOAT *mssq);

extern void DumpHistogram(FILE *hfile, Hist *H, char *Description, FFLOAT mssq);

extern void FreeHistogram(Hist *H);

extern FFLOAT GetTotalBlocks(Hist *H);


extern void WriteDescription(char *s, Image *Im, int unum1, int unum2);

extern void InitOptimJob(OptimJob *job);

extern void GetParams(OptimJob *job, int argc, char *argv[]);

extern void Usage(void);

extern void BriefUsage(void);

extern boolean GetActualBpp(OptimJob *job, char *qfile, FFLOAT *ans);

extern boolean Compress(OptimJob *job, char *qfile, char *cfile, int *nbytes);

extern int newline(char s[], FILE *fp);
extern int newlinefd(char s[], int fd);

extern void FatalError(char *s);

extern void ReadIntTable(FILE *fp, int *tab);

extern void ReadRealTable(FILE *fp, FFLOAT *tab);

extern void WriteIntTable(FILE *fp, int *tab,  char *prefix);

extern void WriteThreshTable(FILE *fp, int *tab,  char *prefix);

extern boolean NotEqualReal(FFLOAT x, FFLOAT y);

extern void PrepareForErrBpp(OptimJob *job, int unum);

extern void lagrPrepareForErrBpp(OptimJob *job, int unum);

extern void SetErrBpp(OptimJob *job, int unum);

extern void lagrSetErrBpp(OptimJob *job, int unum);

extern void lagrSortErrBpp(OptimJob *job, int unum);

extern void Optimize(OptimJob *job, int unum);

extern void CombineUnits(OptimJob *job);

extern void Epilogue(OptimJob *job);

extern void lagrEpilogue(OptimJob *job);

extern void DumpJobChars(OptimJob *job, FILE *fp);

extern void Command(OptimJob *job);

extern void lagrCommand(OptimJob *job);

extern boolean lagrFindRate(OptimJob *job, FFLOAT *bpp, FFLOAT *err, FFLOAT delta, FFLOAT lambdamin, FFLOAT lambdamax, FFLOAT *lambda, int (*anstab)[64]);

extern boolean lagrFindErr(OptimJob *job, FFLOAT *bpp, FFLOAT *err, FFLOAT delta, FFLOAT lambdamin, FFLOAT lambdamax, FFLOAT *lambda, int (*anstab)[64]);

extern void WriteErrBppDist(FILE *fp, int bppscale, int *bpptab,
                            FFLOAT *errtab, char *prefix);

extern void WriteErrFBppDist(FILE *fp, FFLOAT *bpptab,
                             FFLOAT *errtab, char *prefix);


extern void SigToRemainingSig(FFLOAT *Sig, FFLOAT *RemainingSig);

extern void SetSignal(Hist *H, FFLOAT *Sig);

extern int ScanRowForErr(FFLOAT *row, FFLOAT err, int least, int most);

extern int ScanRowForBpp(FFLOAT *row, int bpp, int least, int most);

extern void GetStartingPoints(int *ToCombine[], int *StartingPt,
                              int Loc, int NumUnits);

extern void RecoverQ(int Loc, int unum, OptimJob *job, int *Q, int *Bpp,
                     FFLOAT *Err);

extern void RecoverQandT(int Loc, int unum, OptimJob *job, int *Q,
                         int *T, int *Bpp, FFLOAT *Err);

extern void lagrRecoverQandT(OptimJob *job, int (*QandT)[64],
                             int (*Q)[64], int (*T)[64], FFLOAT (*Bpp)[64], FFLOAT (*Err)[64]);

extern char * FileNameTail(char *s);


