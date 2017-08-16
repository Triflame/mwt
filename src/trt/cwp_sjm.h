#ifndef CWP_SJM
#define CWP_SJM

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>

#include <fcntl.h>      /* non-ANSI */
#include <unistd.h>     /* non-ANSI */
#include <sys/types.h>  /* non-ANSI */

#ifdef CADDR_T_NOT_DEFINED
typedef char *          caddr_t;
#endif

/* TYPEDEFS */
#ifdef CRAY
typedef struct _complexStruct { /* complex number */
	float r,i;
}  cwp_complex;
typedef struct _dcomplexStruct { /* double-precision complex number */
	double r,i;
}  cwp_dcomplex;
#define complex cwp_complex
#define dcomplex cwp_dcomplex
#define cadd cwp_cadd
#define csub cwp_csub
#define cmul cwp_cmul
#define cdiv cwp_cdiv
#define rcabs cwp_rcabs
#define cmplx cwp_cmplx
#define conjg cwp_conjg
#define cneg cwp_cneg
#define cinv cwp_cinv
#define csqrt cwp_csqrt
#define cexp cwp_cexp
#define crmul cwp_crmul
#define cipow cwp_cipow
#define crpow cwp_crpow
#define rcpow cwp_rcpow
#define ccpow cwp_ccpow
#define ccos cwp_ccos
#define csin cwp_csin
#define ccosh cwp_ccosh
#define csinh cwp_csinh
#define cexp1 cwp_cexp1
#define clog cwp_clog

#else
 
#ifndef __cplusplus /* if not C++, define the C struct complex */
#ifndef complex
typedef struct _complexStruct { /* complex number */
	float r,i;
} complex;
#endif/* complex */

#ifndef dcomplex
typedef struct _dcomplexStruct { /* double-precision complex number */
	double r,i;
} dcomplex;
#endif/* dcomplex */

#else /* if C++, define the C++ class complex */
#include "Complex.h"
#endif /* C++ */

#endif

/* DEFINES */
/* uncomment the next block if you are installing */
/* under ultrix, but not using the GCC compiler */

/*
#ifdef ultrix
#define const
#endef
*/

/* CWP ENDIAN */
#ifdef CWP_BIG_ENDIAN
#define CWPENDIAN 1
#endif
#ifdef CWP_LITTLE_ENDIAN
#define CWPENDIAN 0
#endif


#ifndef NULL
#define NULL	((void *)0)
#endif
#ifndef EXIT_FAILURE
#define EXIT_FAILURE (1)
#endif
#ifndef EXIT_SUCCESS
#define EXIT_SUCCESS (0)
#endif
#ifndef SEEK_SET
#define SEEK_SET (0)
#endif
#ifndef SEEK_CUR
#define SEEK_CUR (1)
#endif
#ifndef SEEK_END
#define SEEK_END (2)
#endif
#ifndef PI
#define PI (3.141592653589793)
#endif
#ifndef GOLDEN_RATIO 
#define GOLDEN_RATIO (1.618034)   /* the golden ratio */
#endif
#ifndef TRUE
#define TRUE (1)
#endif
#ifndef FALSE
#define FALSE (0)
#endif
#ifndef YES
#define YES (1)
#endif
#ifndef NO
#define NO (0)
#endif
#ifndef SGN
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#endif
#ifndef POW2
#define POW2(x) ((x)*(x))
#endif
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#ifndef MOD
#define MOD(x,y) ((int)fmod((x),(y)))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define CLOSETO(x, y) ((ABS((x) - (y)) <= FLT_EPSILON*ABS(y))?cwp_true:cwp_false)
#define ISODD(n) ((n) & 01)
#define ISIZE sizeof(int)
#define FSIZE sizeof(float)
#define DSIZE sizeof(double)
#define	STREQ(s,t) (strcmp(s,t) == 0)
#define	STRLT(s,t) (strcmp(s,t) < 0)
#define	STRGT(s,t) (strcmp(s,t) > 0)
#define	DIM(a) (sizeof(a)/sizeof(a[0]))
#ifndef EQ
#define EQ(x,y) ((ABS((x) - (y)) <= DBL_MIN) ? TRUE : FALSE) 
#endif
#ifndef NE
#define NE(x,y) ((ABS((x) - (y)) > DBL_MIN) ? TRUE : FALSE) 
#endif
#ifndef LE
#define LE(x,y) (((x) <= (y)+DBL_MIN) ? TRUE : FALSE) 
#endif
#ifndef GE
#define GE(x,y) (((x) >= (y)-DBL_MIN) ? TRUE : FALSE) 
#endif
#ifndef BT
#define BT(x,y,z) ((((x) < (y)) && ((y) < (z))) ? TRUE : FALSE) 
#endif
#ifndef BTE
#define BTE(x,y,z) ((LE((x),(y)) && LE((y),(z))) ? TRUE : FALSE) 
#endif
#ifndef RAD
#define RAD(x) ((x)/180.0*PI) 
#endif


/* FUNCTION PROTOTYPES */

#ifdef __cplusplus /* if C++, specify external linkage to C functions */
extern "C" {
#endif

#ifdef __cplusplus /* if C++, end external linkage specification */
}
#endif


#endif
