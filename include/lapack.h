#ifndef _LAPACK_
#define _LAPACK_
extern double ddot_(int* N, double *X, int* incx, double *Y, int *incy);
extern void dgemm_( char *transa, char *transb, int *m, int *n, int *k,
                    double * alpha, double * A, int * lda, double * B,
                    int * ldb, double * beta, double * C, int * ldc);
extern double dnrm2_(int *N, double *X, int *incx);
extern void dcopy_(int *N, double *X, int *incx, double *Y, int *incy );
extern void dlarnv_(int *idist, int *iseed, int *n, double *X);
extern void daxpy_(int* n, double *da, double *dx, int *incx, double *dy, int *incy);
extern void dscal_(int *n, double *da, double *dx, int *incx);
#endif
