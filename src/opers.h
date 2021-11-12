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

#define Quantize(r,q) \
  (((r) >= 0.0) ? RoundOff((r)/((FFLOAT) (q))) : \
  ((int) 0 - RoundOff((r)/((FFLOAT) (0-(q))))))

#define QuantizeDis(v,q) \
  (((v) >= 0) ? RoundOff(((FFLOAT) (v))/((FFLOAT) ((q)*2))) : \
  ((int) 0 - RoundOff(((FFLOAT) (v))/((FFLOAT) (0-((q)*2))))))


#define IntUnQuantize(v,q) ((v)*(q))

#define RealUnQuantize(v,q) ((FFLOAT) ((v)*(q)))


