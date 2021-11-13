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
#define STRLENMAX 400
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
#define IM_NOT_RAW 9 /* not an actual kind */

#define ImKindString(k) \
  ((k==IM_RAW)? "RAW" : \
  ((k==IM_PGM_ASCII)? "PGM (ascii)" : \
  ((k==IM_PGM_RAW)? "PGM (raw)" : \
  ((k==IM_PPM_ASCII)? "PPM (ascii)" : \
  ((k==IM_PPM_RAW)? "PPM (raw)" : \
  ((k==IM_PGM)? "PGM" : \
  ((k==IM_PPM)? "PPM" : \
  ((k==IM_PNM)? "PNM" : \
  ((k==IM_NOT_RAW)? "Not-RAW" : \
		"Unknown")))))))))

/** ImKind inclusions:
  apart from <kind> < <kind>,

  PGM_ASCII, PGM_RAW < PGM, PNM
  PPM_ASCII, PPM_RAW < PPM, PNM
   everything < NOT_RAW, except RAW
**/

#define ImKindSubSetOf(a,b) \
     ((a==b) ? 1 : \
     ((((a==IM_PGM_ASCII) || (a==IM_PGM_RAW)) && \
       ((b==IM_PGM) || (b==IM_PNM))) ? 1 : \
     ((((a==IM_PPM_ASCII) || (a==IM_PPM_RAW)) && \
       ((b==IM_PPM) || (b==IM_PNM))) ? 1 : \
     (((b==IM_NOT_RAW) && (a!=IM_RAW) && (a!=IM_UNKNOWN)) ? 1 : \
	     0))))



typedef enum ColorConvTypeEnum {NONE,RGBtoYCC,RGBto2YCC,YCCtoRGB,YCC2toRGB} ColorConvType;
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

    /** info for saving image **/
    int SaveImKind; /* Default: IM_RAW options: IM_PNM or IM_PGM or IM_PPM */
    ColorConvType SaveColorConv; /* default: NONE options YCCtoRGB, YCC2toRGB */

    FILE * errfile; /* fatal error messages printed here.. stderr by default */

    int ImFileFd;
    int UserGaveFd;
    unsigned long ImBytes; /* needs to be set if and only if all the
		       following are true:
		       1. Image is coming on stdin
		       2. Image needs a filter to be forked
			  (not PNM and not RAW, currently)
		       3. Other end of pipe is *not* going
			  to be closed immediately after
			  sending the image data (for example,
			  if the other end is going to read
			  some info back on the same fd) */
    /********** PRIVATE DATA*****************/
    boolean ImExists[MAXCOMPONENTS];
    char ImFilterUsed[STRLENMAX];
    unsigned char FirstFiveChars[5];
    Pixel *Im[MAXCOMPONENTS];
    /* each 2-D component will be stored in row-major order in the
       appropriate Im[][] array
    */
};


typedef struct ImageStruct Image;


static int GetImKind(Image *Im);

static void PNMGetParams(Image *Im, int *maxv);

static int PNMReadImage(Image *Im);


/********************** prototypes ************************/


extern void PrintImgChars(Image *Im);

extern void InitImage(Image *Im);
extern int PeekImage(Image *Im);

extern int ReadImgComp(Image *Im, int cnum);

extern void FreeImgComp(Image *Im, int cnum);

extern int SaveImg(Image *Im, char * fname);

extern int SaveImgComp(Image *Im, int cnum, char * fname); /* will not do conversions,
						    sampling */


extern int DoSubSampling(Image *Im, int cnum);

#endif /* IMAGE_H_INCLUDED */


