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
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>



#include "Image.h"

static void FatalError(char *s)
{
    fprintf(stderr,"%s\n",s);
    exit(1);
}

static int MyRead(int fd, unsigned char *buff, int num)
{
    unsigned char *temp;
    int tocome = num;
    int thistime;

    temp = buff;

    while (tocome > 0)
    {
        thistime = read(fd,temp,tocome);
        if (thistime <= 0) return(num-tocome);
        tocome -= thistime;
        temp = temp + thistime;
    }

    return(num);
}


/***********************************************
Dependencies inside the image structure:
   ColConvKind = RGBtoYCC implies
     NumComponents = 3
     InSampFac[0,1,2][0,1] = 1
   ColConvKind = RGBto2YCC implies
     NumComponents = 3
     InSampFac[0,1,2][0,1] = 1
     OutSampFac[0][0,1] = 1
     OutSampFac[1,2][0,1] = 2
   ImKind = IM_PGM_* implies
     ColConvKind = NONE
     NumComponents = 1
   ImKind = IM_PPM_* implies
     NumComponents = 3
***********************************************/

extern void InitImage(Image *Im)
{
    int i;
    Im->NumRows = 512;
    Im->NumCols = 512;
    Im->NumComponents = 1;
    for (i=0; i<MAXCOMPONENTS; i++)
    {
        Im->ImExists[i] = FALSE;
        Im->InSamplingFactor[i][0] = 1;
        Im->InSamplingFactor[i][1] = 1;
        Im->OutSamplingFactor[i][0] = 1;
        Im->OutSamplingFactor[i][1] = 1;
    }
    Im->ImKind = IM_UNKNOWN;
    Im->ImFileFd = 0;
    Im->ImFileName[0] = '\0';

    Im->ColorConvKind = NONE;
    Im->IsErrImage = FALSE;
    Im->Silent = FALSE;
}

extern void PrintImgChars(Image *Im)
{
    int i;
    if (!Im->Silent)
    {
        fprintf(stderr,"%s: %s image file: ",
                ((Im->ImFileName[0]=='\0')?"<stdin>":Im->ImFileName),
                ImKindString(Im->ImKind));
        fprintf(stderr,"%d x %d / ",Im->NumCols,Im->NumRows);
        for (i=0; i<Im->NumComponents-1; i++)
        {
            fprintf(stderr, "(%d x %d):",Im->InSamplingFactor[i][1],
                    Im->InSamplingFactor[i][0]);
        }
        fprintf(stderr, "(%d x %d)\n",Im->InSamplingFactor[Im->NumComponents-1][1],
                Im->InSamplingFactor[Im->NumComponents-1][0]);
    }
}


extern void PeekImage(Image *Im)
{
    int maxval,i;
    int tempKind;


    if (Im->ImFileName[0] != '\0')
    {
        if ((Im->ImFileFd = open(Im->ImFileName, O_RDONLY, 0)) < 0)
        {
            FatalError("Could not open image file");
        }
    }
    /* check the first few characters. if recognizable
       format, then use the info to fill values
       in Im->.. appropriately.
       also, use Im->ColorConvKind
    */

    tempKind = GetImKind(Im);

    if (Im->ImKind != IM_UNKNOWN)
    {
        if (!ImKindSubSetOf(tempKind,Im->ImKind))
        {
            fprintf(stderr,"Image kind (%s) not a subset of (%s)",
                    ImKindString(tempKind),ImKindString(Im->ImKind));
            FatalError("");
        }
    }
    Im->ImKind = tempKind;

    if ((SAMPLEBITS > 8) && (Im->ImKind != IM_RAW))
    {
        fprintf(stderr,"Image kind (%s) does not permit %d-bit samples",
                ImKindString(Im->ImKind),SAMPLEBITS);
        FatalError("");
    }

    if ((Im->ImKind == IM_PGM_ASCII) || (Im->ImKind == IM_PGM_RAW))
    {
        Im->NumComponents = 1;
        Im->ColorConvKind = NONE;
        PNMGetParams(Im, &maxval);
        if (maxval > MAXSAMPLE)
            FatalError("Maxval in PGM file too high");
    }
    else if ((Im->ImKind == IM_PPM_ASCII) || (Im->ImKind == IM_PPM_RAW))
    {
        Im->NumComponents = 3;
        PNMGetParams(Im, &maxval);
        if (maxval > MAXSAMPLE)
            FatalError("Maxval in PGM file too high");
    }
    else if (Im->ImKind == IM_RAW)
    {
        /*** no action needed ***/
    }
    else
    {
        FatalError("Unknown image file format");
    }
    if (Im->ColorConvKind == RGBtoYCC)
    {
        if (Im->NumComponents != 3)
            FatalError("Number of components must be 3 for RGBtoYCC");
        for (i=0; i<3; i++)
        {
            if ((Im->InSamplingFactor[i][0] != 1) ||
                    (Im->InSamplingFactor[i][1] != 1))
                FatalError("Sampling factors must be all 1 for RGBtoYCC");
            Im->OutSamplingFactor[i][0] = 1;
            Im->OutSamplingFactor[i][1] = 1;
        }
    }
    if (Im->ColorConvKind == RGBto2YCC)
    {
        if (Im->NumComponents != 3)
            FatalError("Number of components must be 3 for RGBto2YCC");
        for (i=1; i<3; i++)
        {
            if ((Im->InSamplingFactor[i][0] != 1) ||
                    (Im->InSamplingFactor[i][1] != 1))
                FatalError("Input sampling factors must be 1 for RGBto2YCC");
            Im->OutSamplingFactor[i][0] = 2;
            Im->OutSamplingFactor[i][1] = 2;
        }
        if ((Im->InSamplingFactor[0][0] != 1) ||
                (Im->InSamplingFactor[0][1] != 1))
            FatalError("Input sampling factors must be 1 for RGBto2YCC");
    }
    PrintImgChars(Im);
}

static void DoRGBtoYCC(Image *Im)
{
    int y,cb,cr;
    int imsize, i;

    imsize = Im->NumRows * Im->NumCols;

    for (i=0; i<imsize; i++)
    {
        y  = ((int) ((double)
                     (((double) Im->Im[0][i]) * ((double) 0.299)) +
                     (((double) Im->Im[1][i]) * ((double) 0.587)) +
                     (((double) Im->Im[2][i]) * ((double) 0.114)) +
                     ((double) 0.5)));
        cb = ((int) ((double)
                     (((double) Im->Im[0][i]) * ((double) -0.16784)) +
                     (((double) Im->Im[1][i]) * ((double) -0.33126)) +
                     (((double) Im->Im[2][i]) * ((double) 0.5))      +
                     ((double) 0.5))) + MAXSAMPLE/2;
        cr = ((int) ((double)
                     (((double) Im->Im[0][i]) * ((double) 0.5)) +
                     (((double) Im->Im[1][i]) * ((double) -0.41869)) +
                     (((double) Im->Im[2][i]) * ((double) -0.08131))      +
                     ((double) 0.5))) + MAXSAMPLE/2;
        if (y > MAXSAMPLE) y = MAXSAMPLE;
        if (cb > MAXSAMPLE) cb = MAXSAMPLE;
        if (cr > MAXSAMPLE) cr = MAXSAMPLE;
        /*
        if (y < 0) y = 0;
        */
        if (cb < 0) cb = 0;
        if (cr < 0) cr = 0;
        Im->Im[0][i] = (Pixel) y;
        Im->Im[1][i] = (Pixel) cb;
        Im->Im[2][i] = (Pixel) cr;
    }
}


static void DoSubSampling(Image *Im, int cnum)
{
    int inrows, incols, inrowskip, incolskip, outrows, outcols;
    int insize, outsize, outbytes;
    int outi,outj,inj,inptrincr;
    Pixel *buffarea, *outptr, *inptr, *inrowptr;
    int boxsize, boxsizeby2;
    int val,i,j;

    outrows = Im->NumRows/Im->OutSamplingFactor[cnum][0];
    outcols = Im->NumCols/Im->OutSamplingFactor[cnum][1];

    outsize = outrows*outcols;

    outbytes = outsize*sizeof(Pixel);
    if ((buffarea = (Pixel *) calloc(1,outbytes))==NULL)
    {
        FatalError("DoSubSampling out of memory");
    }

    inrows = Im->NumRows/Im->InSamplingFactor[cnum][0];
    incols = Im->NumCols/Im->InSamplingFactor[cnum][1];

    inrowskip = inrows/outrows;
    incolskip = incols/outcols;
    boxsize = inrowskip*incolskip;
    boxsizeby2 = boxsize/2;
    inptrincr = inrowskip*incols;

    outptr = buffarea;
    inptr = Im->Im[cnum];
    for (outi=0; outi < outrows; outi++)
    {
        for (outj=0,inj=0; outj < outcols; outj++,inj+=incolskip)
        {
            val = boxsizeby2;
            inrowptr = inptr;
            for (i=0; i< inrowskip; i++)
            {
                for (j=0; j< incolskip; j++)
                {
                    val += (int) inrowptr[inj+j];
                }
                inrowptr += incols;
            }
            val /= boxsize;
            *outptr++ = ((Pixel) val);
        }
        inptr += inptrincr;
    }

    free(Im->Im[cnum]);
    Im->Im[cnum] = buffarea;
}

static void ReadRawComponent(Image *Im, int cnum)
{
    int i,imsize, imbytes;
    unsigned char *buffarea;
    unsigned short temp;
    unsigned int itemp;


    imsize = Im->NumRows * Im->NumCols/
             (Im->InSamplingFactor[cnum][0] * Im->InSamplingFactor[cnum][1]);
    imbytes = imsize*sizeof(Pixel);
    if ((Im->Im[cnum] = (Pixel *) calloc(1,imbytes))==NULL)
    {
        FatalError("ReadImgComp out of memory");
    }

    buffarea = (unsigned char *) Im->Im[0];

    if (cnum==0)
    {
        /* FirstTwoChars were already read in */
        buffarea[0] = Im->FirstTwoChars[0];
        buffarea[1] = Im->FirstTwoChars[1];
        buffarea = buffarea + 2;
        imbytes -= 2;
    }

    if (MyRead(Im->ImFileFd,buffarea,imbytes) != imbytes)
    {
        FatalError("Raw image file seems to be too small");
    }

    if ((cnum==(Im->NumComponents - 1)) && (Im->ImFileFd != 0))
    {
        close(Im->ImFileFd);
    }

    Im->ImExists[cnum] = TRUE;

#ifdef SHORTFLIPPED
#if  (SAMPLEBITS > 8)
    /** type Pixel is unsigned short **/
    for (i=0; i<imsize; i++)
    {
        itemp = (unsigned int) Im->Im[cnum][i];
        temp = (unsigned short) ((unsigned int) (itemp & 0xFF) << 8);
        temp += (unsigned short) ((unsigned int) itemp >> 8);
        Im->Im[cnum][i] = temp;
    }
#endif
#endif


}



extern void ReadImgComp(Image *Im, int cnum)
{
    int imsize, imbytes,i;
    int LastRead;
    unsigned short temp;
    unsigned int itemp;
    unsigned char *buffarea;
    int numloops, remaining;

    if (Im->ImExists[cnum])
    {
        return;
    }


    switch(Im->ImKind)
    {
    case IM_RAW:
        ReadRawComponent(Im,cnum);
        LastRead = cnum;
        if ((Im->ColorConvKind == RGBtoYCC) ||
                (Im->ColorConvKind == RGBto2YCC))
        {
            /* cnum better be 0 */
            if (cnum != 0) FatalError("Logical error in ReadImg!!");
            ReadRawComponent(Im,1);
            ReadRawComponent(Im,2);
            LastRead = 2;
        }
        break;
    case IM_PGM_RAW:
    case IM_PGM_ASCII:
    case IM_PPM_RAW:
    case IM_PPM_ASCII:
        if (cnum == 0)
        {
            imsize = Im->NumRows * Im->NumCols;
            imbytes = imsize*sizeof(Pixel);
            for (i=0; i<Im->NumComponents; i++)
            {
                if ((Im->Im[i] = (Pixel *) calloc(1,imbytes))==NULL)
                {
                    FatalError("ReadImgComp out of memory");
                }
                Im->ImExists[i] = TRUE;
            }
            PNMReadImage(Im);
            if (Im->ImFileFd != 0) close(Im->ImFileFd);
            LastRead = Im->NumComponents-1;
        }
        else
        {
            FatalError("Logical error in ReadImg!!");
        }
        break;
    default:
        FatalError("Logical error in ReadImg!");
    }
    if ((Im->ColorConvKind == RGBtoYCC) ||
            (Im->ColorConvKind == RGBto2YCC))
    {
        DoRGBtoYCC(Im);
    }
    for (i=cnum; i<=LastRead; i++)
    {
        if ((Im->InSamplingFactor[i][0] != Im->OutSamplingFactor[i][0]) ||
                (Im->InSamplingFactor[i][0] != Im->OutSamplingFactor[i][0]))
        {
            DoSubSampling(Im,i);
        }
    }

}

extern void FreeImgComp(Image *Im, int cnum)
{
    Im->ImExists[cnum] = FALSE;
    free(Im->Im[cnum]);
}


extern int GetImKind(Image *Im)
{
    int nread;

    nread = MyRead(Im->ImFileFd, Im->FirstTwoChars, 2);
    if (nread != 2) return(IM_UNKNOWN);
    if (Im->FirstTwoChars[0] != 'P') return(IM_RAW);
    if (Im->FirstTwoChars[1] == '2') return(IM_PGM_ASCII);
    if (Im->FirstTwoChars[1] == '5') return(IM_PGM_RAW);
    if (Im->FirstTwoChars[1] == '3') return(IM_PPM_ASCII);
    if (Im->FirstTwoChars[1] == '6') return(IM_PPM_RAW);
    return(IM_RAW);
}

static int getaline(char s[], int fd)
{
    int i;
    unsigned char c;
    i = 0;
    while ( MyRead(fd,&c,1) == 1 )
    {
        if (c == '\n')
        {
            s[i] = '\0';
            return (i+1);
        }
        s[i] = c;
        i++;
        if (i==(STRLENMAX-1))
        {
            s[i] = '\0';
            return (i);
        }
    }
    return (i);
}

static void PNMGetParams(Image *Im, int *maxv)
{
    char line[STRLENMAX];
    int temp, numread = 0;
    int imfile;

    imfile = Im->ImFileFd;

    /*** file must have been confirmed earlier as a PNM file */
    while (numread < 3)
    {

        getaline(line,imfile);
        while (line[0] == '#') getaline(line,imfile);
        if ((line[0] == 'P') && ('1' < line[1]) &&
                ('7' > line[1]) && (!numread))
        {
            line[0] = ' ';
            line[1] = ' ';
        }
        if (!numread) temp = sscanf(line,"%d%d%d",
                                        &Im->NumCols,&Im->NumRows,maxv);
        else if (numread==1) temp = sscanf(line,"%d%d",
                                               &Im->NumRows,maxv);
        else /* (numread==2) */ temp = sscanf(line,"%d",maxv);
        if (temp < 0) temp = 0;
        numread += temp;
    }
}

static void PNMReadImage(Image *Im)
{
    int i,j,n,val;
    unsigned char c;
    unsigned char buffer[1011]; /* buffer size multpl of 3 */
    int numiters, remaining, totalbytes;
    int kind,  imfile, ht, wd;
    unsigned char *plane1, *plane2, *plane3, *buffarea;

#if (SAMPLEBITS==8)

    kind = Im->ImKind;
    imfile = Im->ImFileFd;
    ht = Im->NumRows;
    wd = Im->NumCols;
    plane1 = Im->Im[0];
    plane2 = Im->Im[1];
    plane3 = Im->Im[2];

    totalbytes = ht*wd;
    if ((kind == IM_PPM_ASCII) || (kind == IM_PPM_RAW))
        totalbytes *= 3;

    /* assumes all the prelude has been read */
    if (kind == IM_PGM_ASCII)
    {
        for (i=0; i<totalbytes; i++)
        {
            c = ' ';
            while ((c==' ') || (c=='\n') || (c=='\t'))
            {
                if (MyRead(imfile,&c,1) != 1) FatalError("Premature EOF in image file");
                while (c== '#')
                {
                    if (MyRead(imfile,&c,1) != 1) FatalError("Premature EOF in image file");
                    while (c != '\n')
                    {
                        if (MyRead(imfile,&c,1) != 1) FatalError("Premature EOF in image file");
                    }
                    if (MyRead(imfile,&c,1) != 1) FatalError("Premature EOF in image file");
                }
            }

            val = 0;
            while ((c >= '0') && (c <= '9'))
            {
                val *= 10;
                val += (c - '0');
                if ((MyRead(imfile,&c,1) != 1) && (i != (totalbytes-1)))
                    FatalError("Premature EOF in image file");

            }
            plane1[i] = (unsigned char) val;
        }
        return;
    }

    if (kind == IM_PGM_RAW)
    {
        if (MyRead(imfile,plane1,totalbytes) != totalbytes)
        {
            FatalError("PGM (raw) file seems too small");
        }
        return;
    }

    if (kind == IM_PPM_ASCII)
    {
        for (i=0; i<totalbytes; i++)
        {
            c = ' ';
            while ((c==' ') || (c=='\n') || (c=='\t'))
            {
                if (MyRead(imfile,&c,1) != 1) FatalError("Premature EOF in image file");
                while (c== '#')
                {
                    if (MyRead(imfile,&c,1) != 1) FatalError("Premature EOF in image file");
                    while (c != '\n')
                    {
                        if (MyRead(imfile,&c,1) != 1) FatalError("Premature EOF in image file");
                    }
                    if (MyRead(imfile,&c,1) != 1) FatalError("Premature EOF in image file");
                }
            }

            val = 0;
            while ((c >= '0') && (c <= '9'))
            {
                val *= 10;
                val += (c - '0');
                if ((MyRead(imfile,&c,1) != 1) && (i != (totalbytes-1)))
                    FatalError("Premature EOF in image file");

            }
            if ((i%3)==0) plane1[i/3] = (unsigned char) val;
            else if ((i%3)==1) plane2[i/3] = (unsigned char) val;
            else /* ((i%3)==2) */ plane3[i/3] = (unsigned char) val;
        }
        return;
    }

    if (kind == IM_PPM_RAW)
    {
        numiters = totalbytes/1011;
        remaining = totalbytes - (numiters*1011);
        n = 0;
        for (i=0; i<numiters; i++)
        {
            if (MyRead(imfile,buffer,1011) < 1011)
                FatalError("Premature EOF in image file");
            j=0;
            while (j<1011)
            {
                plane1[n] = buffer[j++];
                plane2[n] = buffer[j++];
                plane3[n] = buffer[j++];
                n++;
            }
        }
        if (MyRead(imfile,buffer,remaining) < remaining)
            FatalError("Premature EOF in image file");
        j=0;
        while (j<remaining)
        {
            plane1[n] = buffer[j++];
            plane2[n] = buffer[j++];
            plane3[n] = buffer[j++];
            n++;
        }
        return;
    }

    FatalError("Unknown PNM file format");

#endif /* SAMPLEBITS==8 */


}




