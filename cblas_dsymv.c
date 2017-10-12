#include"matrix.h"
#include"lapack.h"

/* YVec = AMat * XVec 
   AMat represents a double symmetric matrix with N*(N+1)/2 elements */
/*---------------------------------------------------------------------*/
void CBLAS_DSYMV(double *YVec, double *AMat, double *XVec, int N)
/*---------------------------------------------------------------------*/
{
    double Alpha, Beta;
    char UPLO[]="L";
  
    int i, j, k;
    int ind;
    int incx, incy;

    Alpha = 1.0;
    Beta  = 0.0;
    incx  = 1;
    incy  = 1;

    double WA[N*N];
    
    k=0;
    for (i=0; i<N; i++) {
        for (j=0; j<=i; j++) {
            /* fortran (i,j) -> ind */
            ind = j * N + i;
            WA[ind] = AMat[k];
            k++;
        }
    }
    /*Yuncai Mei: WA is a lower triangular matrix saved in column major,
     *Upper part of WA is zero.*/
    
    dsymv_(UPLO, &N, &Alpha, WA, &N, XVec, &incx, &Beta, YVec, &incy);
}
     

