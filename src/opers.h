
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

/* This file is included by rdopt.h */

/* Rounding-off for positive numbers */
#define RoundOff(r) \
    ((int) ((FFLOAT) (r) + 0.5))

/* Truncation for positive numbers */
#define Truncate(r) \
    ((int) ((FFLOAT) (r)))

/* Discretization of a real number, for storing the histogram */
#define Discretize(r) \
   (((r) >= 0.0) ? Truncate((r)*2.0) : \
   ((int) 0 - Truncate((r)*(-2.0))))

/* Undscretization: the positive flag is needed because 0
   can be either undiscretized to 0.25 or to -0.25 */
#define UnDiscretize(v,positive) \
   ((FFLOAT) (((FFLOAT) (v))/2.0) + (((positive)) ? 0.25 : -0.25))

#define UnDiscretizePlus(v) \
   ((FFLOAT) (((FFLOAT) (v))/2.0) + (0.25))

#define UnDiscretizeMinus(v) \
   ((FFLOAT) (((FFLOAT) (v))/2.0) - (0.25))

#define Quantize(r,q) \
  (((r) >= 0.0) ? RoundOff((r)/((FFLOAT) (q))) : \
  ((int) 0 - RoundOff((r)/((FFLOAT) (0-(q))))))

#define QuantizeDis(v,q) \
  (((v) >= 0) ? RoundOff(((FFLOAT) (v))/((FFLOAT) ((q)*2))) : \
  ((int) 0 - RoundOff(((FFLOAT) (v))/((FFLOAT) (0-((q)*2))))))


#define IntUnQuantize(v,q) ((v)*(q))

#define RealUnQuantize(v,q) ((FFLOAT) ((v)*(q)))

#define QuantToLowDis(q,qval) ((((qval)<<1)-1)*(q))
#define QuantToHighDis(q,qval) (((((qval)<<1)+1)*(q)) -1)


