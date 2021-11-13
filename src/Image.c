#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fcntl.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <limits.h>
#include <unistd.h>
#include <signal.h>

#ifndef SA_NOCLDWAIT
#define SA_NOCLDWAIT 0
#endif


#include "Image.h"


static int InnerFatalError(FILE *fp, char *s)
{
    fprintf(fp,"%s\n",s);
    return(0);
}

#define FatalError(s) return(InnerFatalError(Im->errfile,(s)))

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

static int MyWrite(int fd, unsigned char *buff, int num)
{
    unsigned char *temp;
    int tocome = num;
    int thistime;

    temp = buff;

    while (tocome > 0)
    {
        thistime = write(fd,temp,tocome);
        if (thistime <= 0) return(num-tocome);
        tocome -= thistime;
        temp = temp + thistime;
    }

    return(num);
}


#define FILTER_BUF_SIZE 4096

static int ForkFilter(char *filter, char *InFile,
                      int InFd, int *OutFd,
                      int NumFirstBytes, char *FirstBytes, int errfd,
                      unsigned long ImBytes)
{
    int pipetofilter[2], pipefromfilter[2];
    int pid, subpid,  len;
    char buf[FILTER_BUF_SIZE];
    static int ForkFiltDidSigAct = 0;
    unsigned long bytes_read;

#ifdef _DEC
    struct sigaction act = {0,0,SA_NOCLDWAIT | SA_NOCLDSTOP};
#else
    struct sigaction act = {0,0,0,SA_NOCLDWAIT | SA_NOCLDSTOP};
#endif

    char *FilterArgs[3];

    if (!ForkFiltDidSigAct)
    {
        /* do not want to create zombies */
        sigaction(SIGCHLD, &act, 0);
        ForkFiltDidSigAct = 1;
    }

    if (InFile)
    {


        /* no need for subfilter */
        if (pipe(pipefromfilter) == -1) return(0);

        if ((pid=fork()) == -1)
        {
            close(pipefromfilter[0]);
            close(pipefromfilter[1]);
            return(0);
        }

        if (pid==0)
        {
            /* child */
            close(pipefromfilter[0]);
            if (dup2(pipefromfilter[1],1) == -1)
            {
                exit(1);
            }
            FilterArgs[0] = filter;
            FilterArgs[1] = InFile;
            FilterArgs[2] = (char *) 0;
            dup2(errfd,2); /* close stderr */
            execvp(FilterArgs[0], FilterArgs);
        }
        else
        {
            /* parent */
            close(pipefromfilter[1]);
            *OutFd = pipefromfilter[0];
            return(1);
        }

    }
    else
    {

        /* If we get here, then InFile was NULL. */
        if (pipe(pipefromfilter) == -1)
        {
            return(0);
        }

        if (pipe(pipetofilter) == -1)
        {
            close(pipefromfilter[0]);
            close(pipefromfilter[1]);
            return(0);
        }


        if ((pid=fork()) == -1)
        {
            close(pipefromfilter[0]);
            close(pipefromfilter[1]);
            close(pipetofilter[0]);
            close(pipetofilter[1]);
            return(0);
        }

        if (pid==0)
        {
            /* child */
            close(pipefromfilter[0]);
            close(pipetofilter[1]);
            if (dup2(pipefromfilter[1],1) == -1)
            {
                exit(1);
            }
            if (dup2(pipetofilter[0],0) == -1)
            {
                exit(1);
            }
            FilterArgs[0] = filter;
            FilterArgs[1] = (char *) 0;

            dup2(errfd,2); /* redirect stderr */
            execvp(FilterArgs[0], FilterArgs);
        }
        else
        {
            /* parent */
            close(pipefromfilter[1]);
            close(pipetofilter[0]);

            if ((subpid=fork()) == -1)
            {
                close(pipefromfilter[0]);
                close(pipetofilter[1]);
                return(0);
            }

            if (subpid==0)
            {
                /* child */
                close(pipefromfilter[0]);
                bytes_read = 0;
                if (NumFirstBytes > 0)
                {
                    MyWrite(pipetofilter[1], FirstBytes,
                            NumFirstBytes);
                    bytes_read += NumFirstBytes;
                }
                while (((ImBytes==0)||
                        (bytes_read < ImBytes)) &&
                        (len=read(InFd, buf, FILTER_BUF_SIZE)) > 0)
                {
                    MyWrite(pipetofilter[1], buf, len);
                    bytes_read += len;
                }
                close(pipetofilter[1]);
                close(InFd);
                exit(0);

            }
            else
            {
                /* parent */
                close(pipetofilter[1]);
                *OutFd = pipefromfilter[0];
                return(1);
            }
        }
    }
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

   SaveImKind = IM_PGM implies
     ColConvKind = NONE
     SaveColConv = NONE
     NumComponents = 1
   SaveColorConv = YCCtoRGB or YCC2toRGB implies
     NumComponents = 3

***********************************************/

extern void InitImage(Image *Im)
{
    int i;
    Im->errfile = stderr;
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
    Im->UserGaveFd = 0;
    Im->ImFileName[0] = '\0';
    Im->SaveImKind = IM_RAW;
    Im->ImBytes = 0;
    strcpy(Im->ImFilterUsed,"none");


    Im->ColorConvKind = NONE;
    Im->SaveColorConv = NONE;

    Im->IsErrImage = FALSE;
    Im->Silent = FALSE;
}

extern void PrintImgChars(Image *Im)
{
    int i;
    fprintf(Im->errfile,"%s: %s image file: ",
            (Im->UserGaveFd ? "<User-Given-Fd>" : (
                 (Im->ImFileFd==0)?"<stdin>": (
                     (Im->ImFileName[0]=='\0')?"<pipe>":Im->ImFileName))),
            ImKindString(Im->ImKind));
    fprintf(Im->errfile,"%d x %d / ",Im->NumCols,Im->NumRows);
    for (i=0; i<Im->NumComponents-1; i++)
    {
        fprintf(Im->errfile, "(%d x %d):",Im->InSamplingFactor[i][1],
                Im->InSamplingFactor[i][0]);
    }
    fprintf(Im->errfile, "(%d x %d) ",Im->InSamplingFactor[Im->NumComponents-1][1],
            Im->InSamplingFactor[Im->NumComponents-1][0]);
    fprintf(Im->errfile, "Filter: %s\n",Im->ImFilterUsed);
    fflush(Im->errfile);
}



extern int PeekImage(Image *Im)
{
    int maxval,i;
    int tempKind;


    if (!Im->UserGaveFd)
    {
        if (Im->ImFileName[0] != '\0')
        {
            if ((Im->ImFileFd = open(Im->ImFileName, O_RDONLY, 0)) < 0)
            {
                FatalError("Could not open image file");
            }
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
            Im->ImKind = IM_UNKNOWN;
            if (!Im->UserGaveFd && Im->ImFileFd) close(Im->ImFileFd);
            return(0);
        }
    }

    Im->ImKind = tempKind;

    if ((SAMPLEBITS > 8) && (Im->ImKind != IM_RAW))
    {
        Im->ImKind = IM_UNKNOWN;
        if (!Im->UserGaveFd && Im->ImFileFd) close(Im->ImFileFd);
        return(0);
    }

    if ((Im->ImKind == IM_PGM_ASCII) || (Im->ImKind == IM_PGM_RAW))
    {
        Im->NumComponents = 1;
        Im->ColorConvKind = NONE;
        PNMGetParams(Im, &maxval);
        if (maxval > MAXSAMPLE)
        {
            Im->ImKind = IM_UNKNOWN;
            if (!Im->UserGaveFd && Im->ImFileFd) close(Im->ImFileFd);
            return(0);
        }
    }
    else if ((Im->ImKind == IM_PPM_ASCII) || (Im->ImKind == IM_PPM_RAW))
    {
        Im->NumComponents = 3;
        PNMGetParams(Im, &maxval);
        if (maxval > MAXSAMPLE)
        {
            Im->ImKind = IM_UNKNOWN;
            if (!Im->UserGaveFd && Im->ImFileFd) close(Im->ImFileFd);
            return(0);
        }
    }
    else if (Im->ImKind == IM_RAW)
    {
        /*** no action needed ***/
    }
    else
    {
        Im->ImKind = IM_UNKNOWN;
        if (!Im->UserGaveFd && Im->ImFileFd) close(Im->ImFileFd);
        return(0);
    }

    if (Im->ColorConvKind == RGBtoYCC)
    {
        if (Im->NumComponents != 3)
        {
            Im->ImKind = IM_UNKNOWN;
            if (!Im->UserGaveFd && Im->ImFileFd) close(Im->ImFileFd);
            return(0);
        }
        for (i=0; i<3; i++)
        {
            if ((Im->InSamplingFactor[i][0] != 1) ||
                    (Im->InSamplingFactor[i][1] != 1))
            {
                Im->ImKind = IM_UNKNOWN;
                if (!Im->UserGaveFd && Im->ImFileFd) close(Im->ImFileFd);
                return(0);
            }
            Im->OutSamplingFactor[i][0] = 1;
            Im->OutSamplingFactor[i][1] = 1;
        }
    }
    if (Im->ColorConvKind == RGBto2YCC)
    {
        if (Im->NumComponents != 3)
        {
            Im->ImKind = IM_UNKNOWN;
            if (!Im->UserGaveFd && Im->ImFileFd) close(Im->ImFileFd);
            return(0);
        }
        for (i=1; i<3; i++)
        {
            if ((Im->InSamplingFactor[i][0] != 1) ||
                    (Im->InSamplingFactor[i][1] != 1))
            {
                Im->ImKind = IM_UNKNOWN;
                if (!Im->UserGaveFd && Im->ImFileFd) close(Im->ImFileFd);
                return(0);
            }
            Im->OutSamplingFactor[i][0] = 2;
            Im->OutSamplingFactor[i][1] = 2;
        }
        if ((Im->InSamplingFactor[0][0] != 1) ||
                (Im->InSamplingFactor[0][1] != 1))
        {
            Im->ImKind = IM_UNKNOWN;
            if (!Im->UserGaveFd && Im->ImFileFd) close(Im->ImFileFd);
            return(0);
        }
    }
    if (!Im->Silent) PrintImgChars(Im);
    return(1);
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

        if (cb < 0) cb = 0;
        if (cr < 0) cr = 0;
        Im->Im[0][i] = (Pixel) y;
        Im->Im[1][i] = (Pixel) cb;
        Im->Im[2][i] = (Pixel) cr;
    }
}


extern int DoSubSampling(Image *Im, int cnum)
{
    int inrows, incols, rowskip, colskip, outrows, outcols;
    int insize, outsize, outbytes;
    int outi,outj,ini,inj,inptrincr,outptrincr;
    Pixel *buffarea, *outptr, *inptr, *inrowptr, *outrowptr;
    int boxsize, boxsizeby2;
    int val,i,j;

    outrows = Im->NumRows/Im->OutSamplingFactor[cnum][0];
    outcols = Im->NumCols/Im->OutSamplingFactor[cnum][1];

    inrows = Im->NumRows/Im->InSamplingFactor[cnum][0];
    incols = Im->NumCols/Im->InSamplingFactor[cnum][1];

    if ((outrows == inrows) && (outcols == incols)) return(1);

    outsize = outrows*outcols;

    outbytes = outsize*sizeof(Pixel);
    if ((buffarea = (Pixel *) calloc(1,outbytes))==NULL)
    {
        FatalError("DoSubSampling out of memory");
    }

    if ((inrows >= outrows) && (incols >= outcols))
    {
        rowskip = inrows/outrows;
        colskip = incols/outcols;
        boxsize = rowskip*colskip;
        boxsizeby2 = boxsize/2;
        inptrincr = rowskip*incols;

        outptr = buffarea;
        inptr = Im->Im[cnum];
        for (outi=0; outi < outrows; outi++)
        {
            for (outj=0,inj=0; outj < outcols; outj++,inj+=colskip)
            {
                val = boxsizeby2;
                inrowptr = inptr;
                for (i=0; i< rowskip; i++)
                {
                    for (j=0; j< colskip; j++)
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
    }
    else if ((inrows <= outrows) && (incols <= outcols))
    {
        /* expand */
        rowskip = outrows/inrows;
        colskip = outcols/incols;
        outptr = buffarea;
        inptr = Im->Im[cnum];
        outptrincr = rowskip*outcols;
        for (ini=0; ini<inrows; ini++)
        {
            for (inj=0,outj=0; inj<incols; inj++,outj+=colskip)
            {

                outrowptr = outptr;
                for (i=0; i< rowskip; i++)
                {
                    for (j=0; j< colskip; j++)
                    {
                        outrowptr[outj + j] = inptr[inj];
                    }
                    outrowptr += outcols;
                }
            }
            outptr += outptrincr;
            inptr += incols;
        }
    }

    free(Im->Im[cnum]);
    Im->Im[cnum] = buffarea;
    return(1);
}

static int ReadRawComponent(Image *Im, int cnum)
{
    int i,imsize, imbytes;
    unsigned char *buffarea;
    unsigned short temp;
    unsigned long itemp;


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
        /* FirstFiveChars were already read in */
        for (i=0; i<5; i++) buffarea[i] = Im->FirstFiveChars[i];
        buffarea = buffarea + 5;
        imbytes -= 5;
    }

    if (MyRead(Im->ImFileFd,buffarea,imbytes) != imbytes)
    {
        FatalError("Raw image file seems to be too small");
    }

    if ((!Im->UserGaveFd) &&
            (cnum==(Im->NumComponents - 1)) && (Im->ImFileFd != 0))
    {
        close(Im->ImFileFd);
    }

    Im->ImExists[cnum] = TRUE;

#ifdef SHORTFLIPPED
#if  (SAMPLEBITS > 8)
    /** type Pixel is unsigned short **/
    for (i=0; i<imsize; i++)
    {
        itemp = (unsigned long) Im->Im[cnum][i];
        temp = (unsigned short) ((unsigned long) (itemp & 0xFF) << 8);
        temp += (unsigned short) ((unsigned long) itemp >> 8);
        Im->Im[cnum][i] = temp;
    }
#endif
#endif

    return(1);
}



extern int ReadImgComp(Image *Im, int cnum)
{
    int imsize, imbytes,i;
    int LastRead;
    unsigned short temp;
    unsigned long itemp;
    unsigned char *buffarea;
    int numloops, remaining;

    if (Im->ImExists[cnum])
    {
        return(1);
    }


    switch(Im->ImKind)
    {
    case IM_RAW:
        if (!ReadRawComponent(Im,cnum)) return(0);
        LastRead = cnum;
        if ((Im->ColorConvKind == RGBtoYCC) ||
                (Im->ColorConvKind == RGBto2YCC))
        {
            /* cnum better be 0 */
            if (cnum != 0) FatalError("Logical error in ReadImg!!");
            if (!ReadRawComponent(Im,1)) return(0);
            if (!ReadRawComponent(Im,2)) return(0);
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
            if (!PNMReadImage(Im)) return(0);
            if ((!Im->UserGaveFd) && (Im->ImFileFd != 0))
                close(Im->ImFileFd);
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
                (Im->InSamplingFactor[i][1] != Im->OutSamplingFactor[i][1]))
        {
            if (!DoSubSampling(Im,i)) return(0);
        }
    }

    return(1);
}

extern void FreeImgComp(Image *Im, int cnum)
{
    if (Im->ImExists[cnum])
    {
        Im->ImExists[cnum] = FALSE;
        free(Im->Im[cnum]);
    }
}

static int ImFilter(Image *Im, char *filter)
{
    int newFd;


    if (Im->ImFileFd && !Im->UserGaveFd)
    {
        /* image was read from a file */
        if (Im->ImFileName[0] == '\0') return(0);
        if (ForkFilter(filter, Im->ImFileName, 0, &newFd, 0, (char *) 0,
                       fileno(Im->errfile), 0))
        {
            close(Im->ImFileFd);
            Im->ImFileFd = newFd;
            strcpy(Im->ImFilterUsed,filter);
            return(1);
        }
        else
        {
            return(0);
        }
    }
    else
    {
        /* must read from existing fd */
        if (ForkFilter(filter, (char *) 0, Im->ImFileFd, &newFd, 5,
                       Im->FirstFiveChars, fileno(Im->errfile), Im->ImBytes))
        {
            Im->ImFileFd = newFd;
            strcpy(Im->ImFilterUsed,filter);
            return(1);
        }
        else
        {
            return(0);
        }

    }

}



extern int GetImKind(Image *Im)
{
    int nread;

    nread = MyRead(Im->ImFileFd, Im->FirstFiveChars, 5);
    if (nread != 5) return(IM_UNKNOWN);
    if (Im->ImKind==IM_RAW) return(IM_RAW);
    if (Im->FirstFiveChars[0] == 'P')
    {
        if (Im->FirstFiveChars[1] == '2') return(IM_PGM_ASCII);
        if (Im->FirstFiveChars[1] == '5') return(IM_PGM_RAW);
        if (Im->FirstFiveChars[1] == '3') return(IM_PPM_ASCII);
        if (Im->FirstFiveChars[1] == '6') return(IM_PPM_RAW);
    }
#ifdef HAVE_GIFTOPNM
    if ((Im->FirstFiveChars[0] == 'G') && (Im->FirstFiveChars[1] == 'I') &&
            (Im->FirstFiveChars[2] == 'F'))
    {
        if (!ImFilter(Im,"giftopnm")) return(IM_UNKNOWN);
        return(GetImKind(Im));
    }
#endif
#ifdef HAVE_TIFFTOPNM
    if (((Im->FirstFiveChars[0] == '\115') &&
            (Im->FirstFiveChars[1] == '\115')) ||
            ((Im->FirstFiveChars[0] == '\111') &&
             (Im->FirstFiveChars[1] == '\111')))
    {
        if (!ImFilter(Im,"tifftopnm")) return(IM_UNKNOWN);
        return(GetImKind(Im));
    }
#endif
#ifdef HAVE_DJPEG
    if ((Im->FirstFiveChars[0] == 0xFF) &&
            (Im->FirstFiveChars[1] == 0xD8) &&
            (Im->FirstFiveChars[2] == 0xFF) &&
            ((Im->FirstFiveChars[3] == 0xE0) ||
             (Im->FirstFiveChars[3] == 0xEE)))
    {
        if (!ImFilter(Im,"djpeg")) return(IM_UNKNOWN);
        return(GetImKind(Im));
    }
#endif
#ifdef HAVE_RASTTOPNM
    if ((Im->FirstFiveChars[0] == 0x59) &&
            (Im->FirstFiveChars[1] == 0xA6) &&
            (Im->FirstFiveChars[2] == 0x6A) &&
            (Im->FirstFiveChars[3] == 0x95))
    {
        if (!ImFilter(Im,"rasttopnm")) return(IM_UNKNOWN);
        return(GetImKind(Im));
    }
#endif
#ifdef HAVE_FITSTOPNM
    if ((Im->FirstFiveChars[0] == 'S') &&
            (Im->FirstFiveChars[1] == 'I') &&
            (Im->FirstFiveChars[2] == 'M') &&
            (Im->FirstFiveChars[3] == 'P') &&
            (Im->FirstFiveChars[4] == 'L'))
    {
        if (!ImFilter(Im,"fitstopnm")) return(IM_UNKNOWN);
        return(GetImKind(Im));
    }
#endif
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
    char *tline;
    int temp, numread = 0;
    int imfile,i,first;

    imfile = Im->ImFileFd;
    first = 1;

    /*** file must have been confirmed earlier as a PNM file */
    while (numread < 3)
    {

        if (first)
        {
            first = 0;
            for (i=0; i<5; i++) line[i] = Im->FirstFiveChars[i];
            tline = &line[0] + 5;
        }
        else
        {
            tline = &line[0];
        }

        getaline(tline,imfile);

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

static int PNMReadImage(Image *Im)
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
        return(1);
    }

    if (kind == IM_PGM_RAW)
    {
        if (MyRead(imfile,plane1,totalbytes) != totalbytes)
        {
            FatalError("PGM (raw) file seems too small");
        }
        return(1);
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
        return(1);
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
        return(1);
    }


#endif /* SAMPLEBITS==8 */


    FatalError("Unknown PNM file format");

}




extern int SaveImgComp(Image *Im, int cnum, char * fname)
{
    int fd, ans;
    unsigned long sz, i, itemp;
    unsigned short temp;

    if (!fname || (fname[0]=='\0') || (!strcmp(fname,"-"))) fd = 1;
    else if ((fd = open(fname, O_CREAT | O_WRONLY | O_TRUNC,
                        S_IRUSR | S_IWUSR)) < 0) return(0);

    sz = Im->NumRows/Im->OutSamplingFactor[cnum][0] *
         Im->NumCols/Im->OutSamplingFactor[cnum][1] * sizeof(Pixel);


#ifdef SHORTFLIPPED
#if  (SAMPLEBITS > 8)
    /** type Pixel is unsigned short so flip before saving **/
    for (i=0; i<sz; i++)
    {
        itemp = (unsigned long) Im->Im[cnum][i];
        temp = (unsigned short) ((unsigned long) (itemp & 0xFF) << 8);
        temp += (unsigned short) ((unsigned long) itemp >> 8);
        Im->Im[cnum][i] = temp;
    }
#endif
#endif

    ans = 1;
    if (MyWrite(fd, (unsigned char *) Im->Im[cnum], sz) < sz) ans = 0;
    if (fd != 1) close(fd);

#ifdef SHORTFLIPPED
#if  (SAMPLEBITS > 8)
    /** type Pixel is unsigned short so flip before saving **/
    for (i=0; i<sz; i++)
    {
        itemp = (unsigned long) Im->Im[cnum][i];
        temp = (unsigned short) ((unsigned long) (itemp & 0xFF) << 8);
        temp += (unsigned short) ((unsigned long) itemp >> 8);
        Im->Im[cnum][i] = temp;
    }
#endif
#endif

    return(ans);
}


static unsigned long SubSampleIndex(int numc, unsigned long idx)
{
    int r, c, rby2, cby2;
    int ncby2;

    r = idx/numc;
    c = idx - (r*numc);
    rby2 = r/2;
    cby2 = c/2;
    ncby2 = numc/2;

    return((((unsigned long) ncby2) * rby2) + cby2);
}

static unsigned char YCCtoR(float y, float cr, float cb)
{
    float fc;
    int ic;

    fc = y;
    fc += (cr - 128.0) * 1.402;
    ic = ((int) ((float) fc + 0.5));
    if (ic < 0) ic = 0;
    if (ic > 255) ic = 255;
    return((unsigned char) ic);
}

static unsigned char YCCtoG(float y, float cr, float cb)
{
    float fc;
    int ic;

    fc =  y;
    fc -= (cb - 128.0) * 0.34414;
    fc -= (cr - 128.0) * 0.71414;

    ic = ((int) ((float) fc + 0.5));
    if (ic < 0) ic = 0;
    if (ic > 255) ic = 255;
    return((unsigned char) ic);
}

static unsigned char YCCtoB(float y, float cr, float cb)
{
    float fc;
    int ic;

    fc = y;
    fc += (cb - 128.0) * 1.772;
    ic = ((int) ((float) fc + 0.5));
    if (ic < 0) ic = 0;
    if (ic > 255) ic = 255;
    return((unsigned char) ic);
}


extern int SaveImg(Image *Im, char * fname)
{
    int fd, ans, cnum;
    unsigned long sz, i, itemp, n, ni, nj;
    unsigned short temp;
    unsigned char buff[3000];
    int numiters, remaining, iter, remby3;
    float y, cr, cb;

    if (!fname || (fname[0]=='\0') || (!strcmp(fname,"-"))) fd = 1;
    else if ((fd = open(fname, O_CREAT | O_WRONLY | O_TRUNC,
                        S_IRUSR | S_IWUSR)) < 0) return(0);


    for (cnum=0; cnum<Im->NumComponents; cnum++)
    {
        if (!Im->ImExists[cnum]) FatalError("SaveImg: some plane nonexistent");
    }

    if (SAMPLEBITS != 8)
    {
        /* must save in raw format */
        ans = 1;
        for (cnum = 0; cnum <Im->NumComponents; cnum++)
        {

            sz = Im->NumRows/Im->OutSamplingFactor[cnum][0] *
                 Im->NumCols/Im->OutSamplingFactor[cnum][1] * sizeof(Pixel);


#ifdef SHORTFLIPPED
#if  (SAMPLEBITS > 8)
            /** type Pixel is unsigned short so flip before saving **/
            for (i=0; i<sz; i++)
            {
                itemp = (unsigned long) Im->Im[cnum][i];
                temp = (unsigned short) ((unsigned long) (itemp & 0xFF) << 8);
                temp += (unsigned short) ((unsigned long) itemp >> 8);
                Im->Im[cnum][i] = temp;
            }
#endif
#endif

            if (MyWrite(fd, (unsigned char *) Im->Im[cnum], sz) < sz) ans = 0;


#ifdef SHORTFLIPPED
#if  (SAMPLEBITS > 8)
            /** type Pixel is unsigned short so flip before saving **/
            for (i=0; i<sz; i++)
            {
                itemp = (unsigned long) Im->Im[cnum][i];
                temp = (unsigned short) ((unsigned long) (itemp & 0xFF) << 8);
                temp += (unsigned short) ((unsigned long) itemp >> 8);
                Im->Im[cnum][i] = temp;
            }
#endif
#endif
        }

        if (fd != 1) close(fd);

        return(ans);
    }

    else if ((Im->SaveImKind == IM_RAW) &&
             ((Im->SaveColorConv == NONE) || (Im->NumComponents != 3)))
    {
        /* SAMPLEBITS == 8 */
        /* must save in raw format */
        ans = 1;
        for (cnum = 0; cnum <Im->NumComponents; cnum++)
        {

            sz = Im->NumRows/Im->OutSamplingFactor[cnum][0] *
                 Im->NumCols/Im->OutSamplingFactor[cnum][1] * sizeof(Pixel);


            if (MyWrite(fd, (unsigned char *) Im->Im[cnum], sz) < sz) ans = 0;


        }

        if (fd != 1) close(fd);

        return(ans);
    }

    else if (Im->SaveImKind == IM_RAW)
    {
        /*  SAMPLEBITS == 8 */
        /*	((Im->SaveColorConv != NONE) && (Im->NumComponents == 3)) */
        /* save as raw after YCCtoRGB or YCC2toRGB */

        ans = 1;

        sz = Im->NumRows * Im->NumCols;

        numiters = sz/3000;
        remaining = sz - (numiters * 3000);

        if ((Im->OutSamplingFactor[0][0] == 1) &&
                (Im->OutSamplingFactor[0][1] == 1) &&
                (Im->OutSamplingFactor[1][0] == 1) &&
                (Im->OutSamplingFactor[1][1] == 1) &&
                (Im->OutSamplingFactor[2][0] == 1) &&
                (Im->OutSamplingFactor[2][1] == 1) &&
                (Im->SaveColorConv == YCCtoRGB))
        {

            /* do the red plane */
            ni = 0;
            for (iter = 0; iter < numiters; iter++)
            {
                for (i=0; i<3000; i++)
                {
                    y  = (float) Im->Im[0][ni];
                    cb = (float) Im->Im[1][ni];
                    cr = (float) Im->Im[2][ni++];
                    buff[i] = YCCtoR(y,cr,cb);
                }
                if (MyWrite(fd, buff, 3000 ) < 3000)
                {
                    ans = 0;
                    break;
                }
            }
            for (i=0; i<remaining; i++)
            {
                y  = (float) Im->Im[0][ni];
                cb = (float) Im->Im[1][ni];
                cr = (float) Im->Im[2][ni++];
                buff[i] = YCCtoR(y,cr,cb);
            }
            if (MyWrite(fd, buff, remaining ) < remaining) ans = 0;

            /* do the green plane */
            ni = 0;
            for (iter = 0; iter < numiters; iter++)
            {
                for (i=0; i<3000; i++)
                {
                    y  = (float) Im->Im[0][ni];
                    cb = (float) Im->Im[1][ni];
                    cr = (float) Im->Im[2][ni++];
                    buff[i] = YCCtoG(y,cr,cb);
                }
                if (MyWrite(fd, buff, 3000 ) < 3000)
                {
                    ans = 0;
                    break;
                }
            }
            for (i=0; i<remaining; i++)
            {
                y  = (float) Im->Im[0][ni];
                cb = (float) Im->Im[1][ni];
                cr = (float) Im->Im[2][ni++];
                buff[i] = YCCtoG(y,cr,cb);
            }
            if (MyWrite(fd, buff, remaining ) < remaining) ans = 0;

            /* do the blue plane */
            ni = 0;
            for (iter = 0; iter < numiters; iter++)
            {
                for (i=0; i<3000; i++)
                {
                    y  = (float) Im->Im[0][ni];
                    cb = (float) Im->Im[1][ni];
                    cr = (float) Im->Im[2][ni++];
                    buff[i] = YCCtoB(y,cr,cb);
                }
                if (MyWrite(fd, buff, 3000 ) < 3000)
                {
                    ans = 0;
                    break;
                }
            }
            for (i=0; i<remaining; i++)
            {
                y  = (float) Im->Im[0][ni];
                cb = (float) Im->Im[1][ni];
                cr = (float) Im->Im[2][ni++];
                buff[i] = YCCtoB(y,cr,cb);
            }
            if (MyWrite(fd, buff, remaining ) < remaining) ans = 0;

        }
        else if ((Im->OutSamplingFactor[0][0] == 1) &&
                 (Im->OutSamplingFactor[0][1] == 1) &&
                 (Im->OutSamplingFactor[1][0] == 2) &&
                 (Im->OutSamplingFactor[1][1] == 2) &&
                 (Im->OutSamplingFactor[2][0] == 2) &&
                 (Im->OutSamplingFactor[2][1] == 2) &&
                 ((Im->SaveColorConv == YCCtoRGB) || (Im->SaveColorConv == YCCtoRGB)))
        {

            /* do the red plane */
            ni = 0;
            for (iter = 0; iter < numiters; iter++)
            {
                for (i=0; i<3000; i++)
                {
                    nj = SubSampleIndex(Im->NumCols, ni); /* inefficient.. sigh */
                    y  = (float) Im->Im[0][ni++];
                    cb = (float) Im->Im[1][nj];
                    cr = (float) Im->Im[2][nj];
                    buff[i] = YCCtoR(y,cr,cb);
                }
                if (MyWrite(fd, buff, 3000 ) < 3000)
                {
                    ans = 0;
                    break;
                }
            }
            for (i=0; i<remaining; i++)
            {
                nj = SubSampleIndex(Im->NumCols, ni); /* inefficient.. sigh */
                y  = (float) Im->Im[0][ni++];
                cb = (float) Im->Im[1][nj];
                cr = (float) Im->Im[2][nj];
                buff[i] = YCCtoR(y,cr,cb);
            }
            if (MyWrite(fd, buff, remaining ) < remaining) ans = 0;

            /* do the green plane */
            ni = 0;
            for (iter = 0; iter < numiters; iter++)
            {
                for (i=0; i<3000; i++)
                {
                    nj = SubSampleIndex(Im->NumCols, ni); /* inefficient.. sigh */
                    y  = (float) Im->Im[0][ni++];
                    cb = (float) Im->Im[1][nj];
                    cr = (float) Im->Im[2][nj];
                    buff[i] = YCCtoG(y,cr,cb);
                }
                if (MyWrite(fd, buff, 3000 ) < 3000)
                {
                    ans = 0;
                    break;
                }
            }
            for (i=0; i<remaining; i++)
            {
                nj = SubSampleIndex(Im->NumCols, ni); /* inefficient.. sigh */
                y  = (float) Im->Im[0][ni++];
                cb = (float) Im->Im[1][nj];
                cr = (float) Im->Im[2][nj];
                buff[i] = YCCtoG(y,cr,cb);
            }
            if (MyWrite(fd, buff, remaining ) < remaining) ans = 0;

            /* do the blue plane */
            ni = 0;
            for (iter = 0; iter < numiters; iter++)
            {
                for (i=0; i<3000; i++)
                {
                    nj = SubSampleIndex(Im->NumCols, ni); /* inefficient.. sigh */
                    y  = (float) Im->Im[0][ni++];
                    cb = (float) Im->Im[1][nj];
                    cr = (float) Im->Im[2][nj];
                    buff[i] = YCCtoB(y,cr,cb);
                }
                if (MyWrite(fd, buff, 3000 ) < 3000)
                {
                    ans = 0;
                    break;
                }
            }
            for (i=0; i<remaining; i++)
            {
                nj = SubSampleIndex(Im->NumCols, ni); /* inefficient.. sigh */
                y  = (float) Im->Im[0][ni++];
                cb = (float) Im->Im[1][nj];
                cr = (float) Im->Im[2][nj];
                buff[i] = YCCtoB(y,cr,cb);
            }
            if (MyWrite(fd, buff, remaining ) < remaining) ans = 0;

        }
        else
        {
            ans = 0;
        }

        if (fd != 1) close(fd);
        return(ans);
    }


    /* 8-bit samples to be saved as PGM or PPM */

    else if (Im->NumComponents == 1)
    {
        /* save as PGM */
        ans = 1;

        sprintf(buff,"P5\n%d %d\n255\n\0",
                Im->NumCols/Im->OutSamplingFactor[0][1],
                Im->NumRows/Im->OutSamplingFactor[0][0]);

        if (MyWrite(fd, buff, (i=strlen(buff))) < i) ans = 0;

        sz = Im->NumRows/Im->OutSamplingFactor[0][0] *
             Im->NumCols/Im->OutSamplingFactor[0][1];

        if (MyWrite(fd, (unsigned char *) Im->Im[0], sz) < sz) ans = 0;

        if (fd != 1) close(fd);

        return(ans);
    }

    else if (Im->NumComponents == 3)
    {

        /* save as PPM */

        /* possibilities allowed:
           1> Im->OutSamplingFactor[0,1,2][0,1] == 1 & Im->SaveColorConv == NONE
           2> Im->OutSamplingFactor[0,1,2][0,1] == 1 & Im->SaveColorConv == YCCtoRGB
           3> Im->OutSamplingFactor[0][0,1] == 1 &
              Im->OutSamplingFactor[1,2][0,1] == 2 &
              Im->SaveColorConv == YCC2toRGB
        */

        ans = 1;

        sprintf(buff,"P6\n%d %d\n255\n\0", Im->NumCols, Im->NumRows);

        if (MyWrite(fd, buff, (i=strlen(buff))) < i) ans = 0;

        sz = Im->NumRows * Im->NumCols * 3;

        numiters = sz/3000;
        remaining = sz - (numiters * 3000);
        remby3 = remaining/3;

        if ((Im->OutSamplingFactor[0][0] == 1) &&
                (Im->OutSamplingFactor[0][1] == 1) &&
                (Im->OutSamplingFactor[1][0] == 1) &&
                (Im->OutSamplingFactor[1][1] == 1) &&
                (Im->OutSamplingFactor[2][0] == 1) &&
                (Im->OutSamplingFactor[2][1] == 1) &&
                (Im->SaveColorConv == NONE))
        {

            ni = 0;

            for (iter = 0; iter < numiters; iter++)
            {
                n = 0;
                for (i=0; i<1000; i++)
                {
                    buff[n++] = Im->Im[0][ni];
                    buff[n++] = Im->Im[1][ni];
                    buff[n++] = Im->Im[2][ni++];
                }
                if (MyWrite(fd, buff, 3000 ) < 3000)
                {
                    ans = 0;
                    break;
                }
            }
            n = 0;
            for (i=0; i<remby3; i++)
            {
                buff[n++] = Im->Im[0][ni];
                buff[n++] = Im->Im[1][ni];
                buff[n++] = Im->Im[2][ni++];
            }
            if (MyWrite(fd, buff, remaining ) < remaining) ans = 0;

        }
        else if ((Im->OutSamplingFactor[0][0] == 1) &&
                 (Im->OutSamplingFactor[0][1] == 1) &&
                 (Im->OutSamplingFactor[1][0] == 1) &&
                 (Im->OutSamplingFactor[1][1] == 1) &&
                 (Im->OutSamplingFactor[2][0] == 1) &&
                 (Im->OutSamplingFactor[2][1] == 1) &&
                 (Im->SaveColorConv == YCCtoRGB))
        {

            ni = 0;

            for (iter = 0; iter < numiters; iter++)
            {
                n = 0;
                for (i=0; i<1000; i++)
                {
                    y  = (float) Im->Im[0][ni];
                    cb = (float) Im->Im[1][ni];
                    cr = (float) Im->Im[2][ni++];
                    buff[n++] = YCCtoR(y,cr,cb);
                    buff[n++] = YCCtoG(y,cr,cb);
                    buff[n++] = YCCtoB(y,cr,cb);
                }
                if (MyWrite(fd, buff, 3000 ) < 3000)
                {
                    ans = 0;
                    break;
                }
            }
            n = 0;
            for (i=0; i<remby3; i++)
            {
                y  = (float) Im->Im[0][ni];
                cb = (float) Im->Im[1][ni];
                cr = (float) Im->Im[2][ni++];
                buff[n++] = YCCtoR(y,cr,cb);
                buff[n++] = YCCtoG(y,cr,cb);
                buff[n++] = YCCtoB(y,cr,cb);
            }
            if (MyWrite(fd, buff, remaining ) < remaining) ans = 0;

        }
        else if ((Im->OutSamplingFactor[0][0] == 1) &&
                 (Im->OutSamplingFactor[0][1] == 1) &&
                 (Im->OutSamplingFactor[1][0] == 2) &&
                 (Im->OutSamplingFactor[1][1] == 2) &&
                 (Im->OutSamplingFactor[2][0] == 2) &&
                 (Im->OutSamplingFactor[2][1] == 2) &&
                 ((Im->SaveColorConv == YCC2toRGB) || (Im->SaveColorConv == YCCtoRGB)))
        {

            ni = 0;

            for (iter = 0; iter < numiters; iter++)
            {
                n = 0;
                for (i=0; i<1000; i++)
                {
                    nj = SubSampleIndex(Im->NumCols, ni); /* inefficient.. sigh */
                    y  = (float) Im->Im[0][ni++];
                    cb = (float) Im->Im[1][nj];
                    cr = (float) Im->Im[2][nj];
                    buff[n++] = YCCtoR(y,cr,cb);
                    buff[n++] = YCCtoG(y,cr,cb);
                    buff[n++] = YCCtoB(y,cr,cb);
                }
                if (MyWrite(fd, buff, 3000 ) < 3000)
                {
                    ans = 0;
                    break;
                }
            }
            n = 0;
            for (i=0; i<remby3; i++)
            {
                nj = SubSampleIndex(Im->NumCols, ni);
                y  = (float) Im->Im[0][ni++];
                cb = (float) Im->Im[1][nj];
                cr = (float) Im->Im[2][nj];
                buff[n++] = YCCtoR(y,cr,cb);
                buff[n++] = YCCtoG(y,cr,cb);
                buff[n++] = YCCtoB(y,cr,cb);
            }
            if (MyWrite(fd, buff, remaining ) < remaining) ans = 0;


        }
        else
        {
            fprintf(Im->errfile,"SaveImg as PPM: Bad sampling factors/ color conversion");
            ans = 0;
        }

        if (fd != 1) close(fd);
        return(ans);
    }


    FatalError("SaveImg: Format unknown/not handled");



}
