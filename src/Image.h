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

#ifndef IMAGE_H_INCLUDED
#define IMAGE_H_INCLUDED
#include <stdio.h>

#include "precision.h"

/* pixel values in the range 0..MAXSAMPLE */
#ifndef MAXSAMPLE

#if (SAMPLEBITS==8)

#define MAXSAMPLE 255
typedef unsigned char Pixel;

#elif (SAMPLEBITS==10)

#define MAXSAMPLE 1023
typedef unsigned short Pixel;

#elif (SAMPLEBITS==12)

#define MAXSAMPLE 4095
typedef unsigned short Pixel;

#elif (SAMPLEBITS==16)

#define MAXSAMPLE 65535
typedef unsigned short Pixel;

#else
/** Force syntax error **/
SAMPLEBITS must be 8,10,12,or16
#endif

#endif /* MAXSAMPLE */

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




#define MAXCOMPONENTS 5

/****** image kinds **********/


#define IM_UNKNOWN 0 /* not an actual kind */
#define IM_RAW  1
#define IM_PGM_ASCII 2
#define IM_PGM_RAW 3
#define IM_PPM_ASCII 4
#define IM_PPM_RAW  5
#define IM_PGM 6     /* not an actual kind */
#define IM_PPM 7     /* not an actual kind */
#define IM_PNM 8     /* not an actual kind */

#define ImKindString(k) \
  ((k==IM_RAW)? "RAW" : \
  ((k==IM_PGM_ASCII)? "PGM (ascii)" : \
  ((k==IM_PGM_RAW)? "PGM (raw)" : \
  ((k==IM_PPM_ASCII)? "PPM (ascii)" : \
  ((k==IM_PPM_RAW)? "PPM (raw)" : \
  ((k==IM_PGM)? "PGM" : \
  ((k==IM_PPM)? "PPM" : \
  ((k==IM_PNM)? "PNM" : \
		"Unknown"))))))))

/** ImKind inclusions:
  apart from <kind> < <kind>,

  PGM_ASCII, PGM_RAW < PGM, PNM
  PPM_ASCII, PPM_RAW < PPM, PNM
**/

#define ImKindSubSetOf(a,b) \
     ((a==b) ? 1 : \
     ((((a==IM_PGM_ASCII) || (a==IM_PGM_RAW)) && \
       ((b==IM_PGM) || (b==IM_PNM))) ? 1 : \
     ((((a==IM_PPM_ASCII) || (a==IM_PPM_RAW)) && \
       ((b==IM_PPM) || (b==IM_PNM))) ? 1 : \
	     0)))



typedef enum ColorConvTypeEnum {NONE,RGBtoYCC,RGBto2YCC} ColorConvType;

/******* images ****************/
struct ImageStruct
{
    char ImFileName[STRLENMAX]; /* Default: '\0', implying stdin */
    int Silent; /* default 0 */
    int NumRows; /* default 512 */
    int NumCols; /* default 512 */
    int NumComponents; /* default 1 */
    int InSamplingFactor[MAXCOMPONENTS][2]; /* default: all 1's */
    /* InSamplingFactor[i][0] is the sampling factor
    in the vertical direction for the ith component,
    for the image being read.
    That is, the ith component has
    NumRows/InSamplingFactor[i][0] rows only.
    Similarly, InSamplingFactor[i][1] is the sampling
    factor in the horizontal direction.
    */
    int OutSamplingFactor[MAXCOMPONENTS][2]; /* default: all 1's */
    /* OutSamplingFactors dictate subsampling to be done
    on the image */
    int ImKind; /* may be set by user, in which case
		 an error will be reported if the actual
		 kind is not the same (or not a subkind) */

    ColorConvType ColorConvKind; /* default: NONE */
    boolean IsErrImage; /* Default: FALSE */

    /********** PRIVATE DATA*****************/
    boolean ImExists[MAXCOMPONENTS];
    int ImFileFd;
    unsigned char FirstTwoChars[2];
    Pixel *Im[MAXCOMPONENTS];
    /* each 2-D component will be stored in row-major order in the
       appropriate Im[][] array
    */
};


typedef struct ImageStruct Image;


/********************** prototypes ************************/


extern void PrintImgChars(Image *Im);

extern void InitImage(Image *Im);
extern void PeekImage(Image *Im);

extern void ReadImgComp(Image *Im, int cnum);

extern void FreeImgComp(Image *Im, int cnum);


static int GetImKind(Image *Im);

static void PNMGetParams(Image *Im, int *maxv);

static void PNMReadImage(Image *Im);

#endif /* IMAGE_H_INCLUDED */


