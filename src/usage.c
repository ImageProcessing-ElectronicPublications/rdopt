
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

extern void Usage(void)
{

    printf("rdopt [flags] {-im imfile [image_flags] | -hist histfile}\n");
    printf("where,\n");
    printf("flags can be:\n\n");
    printf("\t -v: verbose output: the more -v's, the more verbose the output\n");
    printf("\t -silent: Normal stderr output will be suppressed\n");
    printf("\t -cmdfile cfile: after optimization, read commands from cfile instead of stdin\n");
    printf("\t -height H: nominal height of the image is H (default 512)\n");
    printf("\t -width W: nominal width of the image is W (default 512)\n");
    printf("\t -planes N: there are N color planes in the image (default 1)\n");
    printf("\t            these color planes are numbered 0..(N-1)\n");
    printf("\t -numtables t: t tables will be used for quantization.\n");
    printf("\t               1 table each for color planes numbered 0..(t-2)\n");
    printf("\t               Color planes numbered (t-1) thru (N-1) will be\n");
    printf("\t               quantized together. (default 1)\n");
    printf("\t -clampDC d: don't consider qtable entries greater than d\n");
    printf("\t             for the DC coefficient (default 255 or QTABENTRYMAX for >8-bit prec)\n");
    printf("\t -dontclampDC: consider all qtable entries for DC too\n");
    printf("\t -dcdpcm:  use an approximation of dpcm coding done for\n");
    printf("\t            the DC coef to estimate its entropy.\n");
    printf("\t -bppmax b: consider bits per pixel values upto b (default 1.0)\n");
    printf("\t -bppscale B: discretize bpp by the integer B (default 5000)\n");
    printf("\t              the slowness of the program grows linearly with b*B\n");
    printf("\t -thresh T: try thresholds upto T/2 greater than q/2. (default 0)\n");
    printf("\t -weights n cwfile: for unit n, the errors in various DCT coefficients\n");
    printf("\t                  are to be weighted according to the weights\n");
    printf("\t                  listed in RMWSS order in cwfile. these weights\n");
    printf("\t                  will be normalized so that they add to 64\n");
    printf("\t                  Useful for using perceptually weighted distortion measures\n");
    printf("\t -pweights w0[,w1,..]: for aggregate error,  unit n\n");
    printf("\t               will have weight wn. Missing weights will be set\n");
    printf("\t               to the last specified weight. weights will be\n");
    printf("\t               normalized to add up to numtables\n");
    printf("\t -subsamp n h w: color plane number n is to be subsampled to height H/h and \n");
    printf("\t                 width W/w. defaults are h=w=1 for each\n");
    printf("\t                 color plane.\n");
    printf("\t -insubsamp n h w: color plane number n in input file has height H/h and \n");
    printf("\t                 width W/w. defaults are h=w=1 for each\n");
    printf("\t                 color plane. H(W) must be divisible by h(w)\n");
    printf("\t -mintable n fname: qtable # n has to have entries no lesser\n");
    printf("\t                    than those in the table listed in\n");
    printf("\t                    row-major white-space-separated (RMWSS)\n");
    printf("\t                    order in the file fname\n");
    printf("\t -maxtable n fname: qtable # n has to have entries no more\n");
    printf("\t                    than those in the table listed in\n");
    printf("\t                    RMWSS order in the file fname\n");
    printf("\t -correct c: find actual bpp at predicted-bpp-c, apply correction thereafter\n");
    printf("\t -plot pfile: dump a plot of bpp-psnr pairs in the\n");
    printf("\t              file pfile\n");
    printf("\t -points n: use n points in -plot (default 20)\n");
    printf("\t -pbppmax b: plot bpp range 0-b (default min of bppmax,1.5)\n");
    printf("\t -bppplane n: useful for images with > 1 planes. Use plane\n");
    printf("\t             number n to calculate # of pixels for reporting\n");
    printf("\t             bits per pixel. Default is to use the sum of\n");
    printf("\t           #pixels in each plane\n");
    printf("\t -errfile fname: send stderr to file fname (can be - for stdout)\n");
    printf("\n\t ** Obscure flags that you shouldn't really need:\n\n");
    printf("\t -stats: histograms will be dumped in the file HISTOGRAM\n");
    printf("\t -onlystats: only histogram-dumping, no optimization\n");
    printf("\t -mapq: qentries tried will increase logarithmically\n");
    printf("\t        This is useful for speeding up rdopt with 12+ bit samples\n\n");
    printf("** Either the -im flag or the -hist flag must be present. In case\n");
    printf("of -hist, the file histfile must contain a histogram generated\n");
    printf("by rdopt -stats or -onlystats. You shouldn't need to use this option, really.\n");
    printf("If the image file name is given as -, it means stdin\n");
    printf("The image file can be raw row-major bytes (or 16-bits for\n");
    printf("10- 12- and 16-bit images) with color planes stored one after\n");
    printf("the other. In this case, the -height, -width, -planes, and -insubsamp\n");
    printf("options must have been correctly set whenever they differ from\n");
    printf("the default values. OR, the image can be in pgm or ppm format\n");
    printf("in which csae none of these options need to be set as their values\n");
    printf("will be read from the image file\n");

    printf("\n** Image_flags can be:\n");
    printf("\t -rgbtoycc: the input image is in RGB color format. these are\n");
    printf("\t              to be converted to YCrCb before use\n");
    printf("\t -rgbto2ycc: convert RGB to YCbCr, subsample 2:1 in horiz and\n");
    printf("\t             vert directions for Cb and Cr planes\n");

    exit(1);
}

extern void BriefUsage(void)
{

    printf("rdopt [flags] {-im imfile [image_flags] | -hist histfile}\n");
    printf("  flags can be:\n");
    printf("	 [-v] [-silent] [-cmdfile cfile] [-height H] [-width W] [-correct c]\n");
    printf("	 [-planes N] [-numtables t] [-clampDC d] [-dontclampDC] [-pweights w0,..]\n");
    printf("	 [-bppmax b] [-bppscale B] [-weights n cwfile] [-subsamp n h w]\n");
    printf("	 [-insubsamp n h w] [-mintable n fname] [-maxtable n fname]\n");
    printf("	 [-plot pfile] [-points n] [-pbppmax b] [-dcdpcm]\n");
    printf("       [-bppplane n] [-errfile fname] [-thresh T] [-help]\n");
    printf("  image_flags can be:\n");
    printf("	 [-rgbtoycc] [-rgbto2ycc]\n");

    exit(1);
}
