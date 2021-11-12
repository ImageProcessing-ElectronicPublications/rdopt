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
#include <math.h>

#include "version.h"
#include "precision.h"

#define CQTAB "RDOPT.Q" VERSION
#define CQHIST "RDOPT.H" VERSION

#include "Image.h"

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

#define PLOT_POINTS 10 /* # of points to plot with DumpPlot set */

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
    FFLOAT *Err[64]; /* ERR can be thrown away after use */
    int *Bpp[MAXCOMPONENTS][64]; /* BPP for each unit needs
				     to be retained */

    boolean WeightedCoefs;
    FFLOAT CoefWeights[64];

    FFLOAT bppWeight[MAXCOMPONENTS]; /* bpp = \sum_{i=0..NumTables-1}
					       bpp-for-ith-unit *
					       bppWeight[i] */
    FFLOAT LogSignalSq[MAXCOMPONENTS];
    FFLOAT LogTotalSignalSq;
    FFLOAT LogPeakSq[MAXCOMPONENTS];
    FFLOAT LogTotalPeakSq;
    FFLOAT NumPixels[MAXCOMPONENTS];
    FFLOAT NumBlocks[MAXCOMPONENTS];
    FFLOAT TotalNumPixels;

    boolean ImFilePresent;
    boolean HistFilePresent;
    char HistFile[STRLENMAX];

    Hist Histogram[64];
    boolean DumpStats;
    boolean OnlyDumpStats;

    boolean BppRangeExists[MAXCOMPONENTS];
    FFLOAT BppRange[MAXCOMPONENTS][64][2];
    int MinTable[MAXCOMPONENTS][64];
    int MaxTable[MAXCOMPONENTS][64];
    int DCclamp;

    boolean DumpPlot;
    char PlotFileName[STRLENMAX];

};

typedef struct OptimJobStruct OptimJob;


#define HEADER_BYTES 200
/** approximate number of header bytes needed for each unit **/



/********************** prototypes ************************/



extern void SetHistogram(Hist *H, Image *Im, int cnum, FFLOAT *mssq);

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

extern int newline(char s[], FILE *fp);
extern int newlinefd(char s[], int fd);

extern void FatalError(char *s);

extern void ReadIntTable(FILE *fp, int *tab);

extern void ReadRealTable(FILE *fp, FFLOAT *tab);

extern void WriteIntTable(FILE *fp, int *tab);

extern boolean NotEqualReal(FFLOAT x, FFLOAT y);

extern void PrepareForErrBpp(OptimJob *job, int unum);

extern void SetErr(OptimJob *job, int unum);

extern void SetBpp(OptimJob *job, int unum);

extern void Optimize(OptimJob *job, int unum);

extern void CombineUnits(OptimJob *job);

extern void Epilogue(OptimJob *job);

extern void DumpJobChars(OptimJob *job, FILE *fp);

extern void Command(OptimJob *job);

extern FFLOAT GetBpp(OptimJob *job, int unum, int n, int q);

extern void TranslateBppRange(OptimJob *job, int unum);

