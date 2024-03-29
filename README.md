![GitHub release (latest by date)](https://img.shields.io/github/v/release/ImageProcessing-ElectronicPublications/rdopt)
![GitHub Release Date](https://img.shields.io/github/release-date/ImageProcessing-ElectronicPublications/rdopt)
![GitHub repo size](https://img.shields.io/github/repo-size/ImageProcessing-ElectronicPublications/rdopt)
![GitHub all releases](https://img.shields.io/github/downloads/ImageProcessing-ElectronicPublications/rdopt/total)
![GitHub](https://img.shields.io/github/license/ImageProcessing-ElectronicPublications/rdopt)

# rdopt

### Description

Rdopt is a program to generate customized DCT quantization tables for use in JPEG compression. The input is an image file (raw row-major bytes or pgm or ppm). The program runs an optimization algorithm and then prompts the user to enter compression or SNR/PSNR specifications whereupon it outputs a Quantization Table that when used in JPEG compression with customized Huffman Tables (eg., Independent JPEG group's `cjpeg` program with the options `cjpeg -optimize -qtables <Q-FILE>`) compresses the image to within those specifications optimally.

**NOTE**: If you avoid reading the CHANGES file, the main thing you need to know is:
```
  use flag "-method lagrangian" for faster speed
  (except when.. heck, read the CHANGES file, or the
  "MORE INFO" part of this file)
```
### Compiling

Copy the appropriate `makefile.{DEC,HP,RS6000,PENTIUM,SUN}` file onto makefile. Edit makefile: especially, set BINPATH to wherever you want to put the executable. Do `make` and pray!

### Usage

```sh
rdopt [flags] {-im imfile [image_flags] | -hist histfile}
```
where, flags can be:
```
     -v: verbose output: the more -v's, the more verbose the output
     -silent: Normal stderr output will be suppressed
     -cmdfile cfile: after optimization, read commands from cfile instead of stdin
     -height H: nominal height of the image is H (default 512)
     -width W: nominal width of the image is W (default 512)
     -planes N: there are N color planes in the image (default 1)
                these color planes are numbered 0..(N-1)
     -numtables t: t tables will be used for quantization.
                   1 table each for color planes numbered 0..(t-2)
                   Color planes numbered (t-1) thru (N-1) will be
                   quantized together (default 1).
     -clampDC d: don't consider qtable entries greater than d
                 for the DC coefficient (default 255 or QTABENTRYMAX for >8-bit prec)
     -dontclampDC: consider all qtable entries for DC too
     -dcdpcm:  use an approximation of dpcm coding done for
                the DC coef to estimate its entropy.
     -bppmax b: consider bits per pixel values upto b (default 1.0)
     -bppscale B: discretize bpp by the integer B (default 5000)
                  the slowness of the program grows linearly with b*B
     -thresh T: try thresholds upto T/2 greater than q/2. (default 0)
     -weights n cwfile: for unit n, the errors in various DCT coefficients
                      are to be weighted according to the weights
                      listed in RMWSS order in cwfile. these weights
                      will be normalized so that they add to 64
                      Useful for using perceptually weighted distortion measures
     -pweights w0[,w1,..]: for aggregate error,  unit n
                   will have weight wn. Missing weights will be set
                   to the last specified weight. weights will be
                   normalized to add up to numtables
     -subsamp n h w: color plane number n is to be subsampled to height H/h and
                     width W/w. defaults are h=w=1 for each
                     color plane.
     -insubsamp n h w: color plane number n in input file has height H/h and
                     width W/w. defaults are h=w=1 for each
                     color plane. H(W) must be divisible by h(w)
     -mintable n fname: qtable # n has to have entries no lesser
                        than those in the table listed in
                        row-major white-space-separated (RMWSS)
                        order in the file fname
     -maxtable n fname: qtable # n has to have entries no more
                        than those in the table listed in
                        RMWSS order in the file fname
     -correct c: find actual bpp at predicted-bpp-c, apply correction thereafter
     -plot pfile: dump a plot of bpp-psnr pairs in the
                  file pfile
     -points n: use n points in -plot (default 20)
     -pbppmax b: plot bpp range 0-b (default min of bppmax,1.5)
     -bppplane n: useful for images with > 1 planes. Use plane
                 number n to calculate # of pixels for reporting
                 bits per pixel. Default is to use the sum of
               #pixels in each plane
     -errfile fname: send stderr to file fname (can be - for stdout)

     ** Obscure flags that you shouldn't really need:

     -stats: histograms will be dumped in the file HISTOGRAM
     -onlystats: only histogram-dumping, no optimization
     -mapq: qentries tried will increase logarithmically
            This is useful for speeding up rdopt with 12+ bit samples
```

** Either the `-im` flag or the `-hist` flag must be present. In case of `-hist`, the file histfile must contain a histogram generated by `rdopt -stats` or `-onlystats`. You shouldn't need to use this option, really.
If the image file name is given as `-`, it means stdin.
The image file can be raw row-major bytes (or 16-bits for 10- 12- and 16-bit images) with color planes stored one after the other. In this case, the `-height`, `-width`, `-planes`, and `-insubsamp` options must have been correctly set whenever they differ from the default values. OR, the image can be in pgm or ppm format in which csae none of these options need to be set as their values will be read from the image file

** Image_flags can be:
```
     -rgbtoycc: the input image is in RGB color format. these are
                  to be converted to YCrCb before use
     -rgbto2ycc: convert RGB to YCbCr, subsample 2:1 in horiz and
                 vert directions for Cb and Cr planes
```

### Command Interface

After the optimization algorithm is done running, `rdopt` prompts the user with
```
Command>
```
unless the `-cmdfile` option was used in which case commands are read from the file named in that option. In response to most commands, quantization tables `Q` (and thresholding tables, if `-thresh` was specified) will be generated.

Valid commands are:
```
     [compress] size <target_size>
       find Q to get compressed size <= target_size bytes

     [compress] bpp <target_bpp>
       find Q to get compressed size <= target_bpp bits per pixel

     [compress] psnr <target_psnr>
       find Q to get psnr >= target_psnr dB

     [compress] snr <target_snr>
       find Q to get snr >= target_snr dB

     [compress] rmse <target_rmse>
       find Q to get rmse <= target_rmse

     qfile fname
       output subsequent next qtables generated by  commands to
       the file fname instead of the default file name (explained below)

     cfile fname
       output subsequent compressed files to the file fname instead
       of the default

     stats
       If this command is issued, all subsequent qtables will
       also contain coefficient-wise breakup of rate and distortion.
       This is useful for splitting the jpeg file into scans to
       use the "progressive" mode. This directory also contains
       a program (QTtoSF.c) which will convert this break-up info
       into a "scans file" that can be used by IJG's cjpeg. See below
       (in the section on supporting software) for details.
       Printing these stats can be disabled by specifying:

     nostats
       Do not print r-d breakups in qtable files (default).

     correct <use_bpp>
       You need cjpeg for this. The qtables found will be used to
       actually compress the image to get the actual bits per pixel
       (real_bpp). Subsequently, all reported bpp's will be corrected
       by the amount real_bpp - use_bpp. This can be disabled with:

     nocorrect
       Do not apply bpp corrections (default).

  A file RDOPT.Q.<command>.<target> will be generated in response
  to the commands "size", "bpp", "psnr", "snr", and "rmse"
  (unless a qfile command has been issued earlier), if the "compress"
  prefix is not given. This file will contain the quantization
  table(s) and the threshoding table(s) (if -thresh was used).
  If the "compress" prefix is specified, then  IJG's cjpeg will be
  used to compress the image according to the quantization table(s)
  found. The compressed file will have the name <imagename>.jpg,
  unless a cfile command has been issued earlier.

  In case the lagrangian method is used, three additional commands
  can be issued at the command interface:

    deltabpp <value>
      While searching for a target bpp, stop as soon as you
      have come within plusminus <value>. Default is 0.001, you
      can set this to zero.

    deltamse <value>
      While searching for a target distortion, stop as soon as you
      have come within plusminus <value> of the mse.
      Default is 0.01, you can set this to zero.

    deltalambda <value>
      While using the bisection method to search for the
      Lagrangian parameter lambda that will give you rate/distortion
      closest to the target, stop when the lambda interval
      becomes smaller than <value>. Default is 1e-8. Avoid tinkering
      with this. Especially, do not set this to zero!
```

## MORE INFO:

When should I use the Lagrangian method for optimization?

  Almost always. This would result in substantial improvement
  in speed. However, for thresholding tmax greater than about
  100 (never really needed), use dynamic programming.

What is unit_num?

  Unit number is a number between 0 and (t-1) where t is the
  number of units in the image. Each unit is a set of color planes
  quantized with the same qtable.

  Examples:
```
  rdopt -bpplane 0 -meth lagr -numtables 2 -im lena.ppm -rgbto2ycc
```
    Here, there are three color planes (Y, Cb, Cr) and
    two units. Unit # 0 consists of the Y plane and
    Unit # 1 consists of the Cb and Cr planes.
```
  rdopt -im lena.pgm
```
    Here, there is just one plane, and one unit.
```
  rdopt -meth lagr -planes 5 -numtables 3 -height 120 -width 120 -im some.raw.image
```
    Here, there are 5 color planes and 3 units. Unit # i has
    only the color plane # i, for i=0,1. Unit # 2 consists of
    color planes #2, #3, and #4.

### Accuracy of results

The `psnr/snr/rmse` values predicted by rdopt are *very* accurate. The bpp predicted are *usually* within 0.02 bits per pixel of the actual rate achieved by any good jpeg compressor using customized huffman tables.

### Using the generated Q tables

The Q tables can be used directly by the cjpeg program of the Independent JPEG group. For example,
```sh
rdopt -im lena.pgm

Command> qfile tab1
Command> bpp 1.0
Command> quit
```
can be followed with
```sh
cjpeg -dct float -optimize -qtables tab1 -grayscale lena.pgm > lena.jpg
```
We recommend the use of `-dct` float and `-optimize` options in the cjpeg program, for accuracy of the psnr/bpp predicted by `rdopt`.

### Using thresholding

For each unit, zeroing thresholds for DCT coefficients can be determined too. For example,
```sh
rdopt -method lagrangian -thresh 20 -im lena.pgm

Command> qfile tab2
Command> bpp 1.0
Command> quit
```
This will produce the file tab2 which will look like:
```
#RDOPT.Qv2.0
#Image read from PGM (raw) file ../../../archive/lena2/IM
#Width x Height:512 x 512 / (1 x 1)
#Color-conversion applied: none
#Number of qtables used: 1
#Optimization parameters (method: Lagrangian):
#        DC clamp = 255, ThreshSpan = 20
#
#Bits per pixel for each unit and the whole image
#    reported as num-bits/262144.000000, the denominator
#    being the sum of num pixels in all planes
#
#Quantization table for color planes 0 thru 0
9 10 9 9 9 10 11 11
9 10 11 9 11 10 11 10
9 10 11 11 10 11 12 10
10 10 10 11 10 10 11 11
11 10 11 10 10 12 11 10
11 10 12 11 11 11 10 6
11 10 10 12 10 11 11 12
12 6 11 11 11 11 4 13
# Table of thresholds
#T   4.5   5.5   5.5   5.5   6.0   6.5   7.5   8.0
#T   4.5   6.0   6.5   6.0   7.0   6.5   7.5   7.5
#T   6.5   6.0   7.0   7.0   7.5   7.5   8.0   8.0
#T   6.5   6.5   7.0   7.5   7.0   7.5   8.0   8.5
#T   7.5   7.0   7.0   7.0   7.5   8.0   8.0   8.5
#T   7.5   7.5   8.0   7.5   8.5   8.5   8.5   9.5
#T   8.5   8.5   8.0   9.0   8.5   9.0   9.0   9.5
#T   9.0   9.5   9.5   9.0   9.0   9.0  10.0  11.5
#Bpp = 1.000247 bits per pixel
#Size = 32776 bytes
#SNR = 33.584820 dB
#PSNR = 39.262947 dB
#RMSE = 2.775830
#END
```
The threshold table is given in lines beginning with
```
#T .......
```
To use the threshold table (and it'd better be used if you want to achieve the specified rate/distortion!), you need to modify your favorite JPEG compressor so that it can read these threshold tables, and zero off every DCT coefficient that is less than the corresponding threshold. IJG's cjpeg can be easily modified to do so: see below and read the file cjpeg_thresh. Note that there will be one threshold table for every unit.

## SUPPORTING SOFTWARE

1. Independent JPEG Group's `cjpeg`

Use ver6 or higher, preferably. `rdopt` will try to fork off `cjpeg` if `-correct` is used, or if a `compress` or `correct` command is issued at the command interface. `cjpeg` does not support thresholding.
A simple modified version of `cjpeg` (to support thresholding) can be built by following the instructions in the file `cjpeg_thresh`. The modified cjpeg can take the command line option `-thresh <tfile>` tfile must contain thresholding table(s), on lines beginning with `#T `. `rdopt` puts quantization and thresholding tables in the same file, but lines containing  thresh tables begin with `#T `.

2. `QTtoSF`

This program converts coefficient-wise rate/distortion statistics present in the qtable file rdopt produces (need to issue a `stats` command at rdopt command interface for this) into progressive scan breakups such that each AC scan will have approximately the same size. The output is a `scans file` that can be used by `cjpeg` in its `-scans <file>` command line option.

Usage is:
```sh
QTtoSF <numACscans> qtabfile [scanfile]
```

### 10- 12- and 16- bit images

`rdopt` can be compiled for these images too, by changing `SAMPLEBITS` and `QTABBITS` in `precision.h` These images can only be in the raw form, stored as 2 bytes per pixel.

For 12- and 16- bit images it's advisable to use `-mapq` option for greater speed.

I haven't really tested rdopt much with 10, 12 or 16 bit images. Comments/bug-reports etc. will be appreciated.


### Comments, Suggestions, Bugs

Please send email to ratnakar@cs.wisc.edu

Viresh Ratnakar
University of Wisconsin-Madison
Computer Sciences Dept
1210 W Dayton St
Madison, WI 53706
