#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<stdarg.h>
#include"vector.h"
#include"lapack.h"
#include"matrix.h"

Dvector *CreateDvector(int dim)
{
    Dvector *pVect;
    pVect = calloc(1, sizeof(Dvector));
    pVect->dim = dim;
    pVect->data = calloc(1, sizeof(double)*dim);
    return pVect;
}

void DeleteDvector(Dvector **pX)
{
    if(pX != NULL)
    {
        free((*pX)->data);
        free(*pX);
        *pX = NULL;
    }
}

void DeleteDvectorList(int count, Dvector **pX, ...)
{
    DeleteDvector(pX);
    int i;
    va_list list;
    va_start(list, pX);
    for(i=0; i < count-1; i++)
        DeleteDvector( va_arg(list, Dvector **) );
    va_end(list);
}

void InitDvector(Dvector *pX, double *pArray)
{
    if(pX == NULL) 
    {
        printf("Initial Dvector failed. Allocate Dvector first!\n");
        return;
    }
    int inc = 1;
    dcopy_(&(pX->dim), pArray, &inc, pX->data, &inc);
}

void ShowDvector(Dvector *pVect)
{
    if(pVect == NULL) 
    {
        printf("Show Dvector failed. Allocate Dvector first!\n");
        return;
    }
    int dim, i;
    dim = pVect->dim;
    for(i=0; i<dim; i++)
    {
        printf("%5f ", pVect->data[i]);
    }
    printf("\n");
}


void DvectorArithmetic(Dvector *pX, double alpha, Dvector *pY,
    double beta, Dvector *pZ)
/* X := alpha*Y + beta*Z;*/
{
   if(pX == NULL || pY == NULL || pZ == NULL)
    {
        printf("DvectorArithmetic: Vector(s) is NULL!\n");
        return;
    }
    if(pY->dim != pZ->dim)
    {
        printf("DvectorArithmetic: Vectors have different size!\n");
        return;
    }
    double temp[pZ->dim];
    int dim_temp = pZ->dim;
    int inc_temp = 1;
    dcopy_(&dim_temp, pZ->data, &inc_temp, temp, &inc_temp);
    dscal_(&dim_temp, &beta, temp, &inc_temp);
    daxpy_(&dim_temp, &alpha, pY->data, &inc_temp, temp, &inc_temp);
    dcopy_(&dim_temp, temp, &inc_temp, pX->data, &inc_temp);
}

void DvectorScalar(Dvector *pX, double alpha, Dvector *pY)
/* X := alpha * Y */
{
    if(pX == NULL || pY == NULL)
    {
        printf("DvectorScalar: Vector(s) is NULL!\n");
        return;
    }
    if(pX->dim != pY->dim)
    {
        printf("DvectorScalar: OutVectors have different size with input vector!\n");
        return;
    }
    double temp[pY->dim];
    int dim_temp = pY->dim;
    int inc_temp = 1;
    dcopy_(&dim_temp, pY->data, &inc_temp, temp, &inc_temp);
    dscal_(&dim_temp, &alpha, temp, &inc_temp);
    dcopy_(&dim_temp, temp, &inc_temp, pX->data, &inc_temp);
}

double DvectorIP(Dvector *pX, Dvector *pY)
    /*Inner product*/
{
    if(pX == NULL || pY == NULL)
    {
        printf("DvectorIP: Vector(s) is NULL!\n");
        return 0;
    }
    if(pX->dim != pY->dim)
    {
        printf("Dvector: Vectors have different size!\n");
        return 0;
    }
    int i=1;
    return ddot_(&(pX->dim), pX->data, &i, pY->data, &i);
}

double Dvector_2Norm(Dvector *pX)
    /*return sqrt(X' * X)*/
{
    if(pX == NULL)
    {
        printf("Dvector_2Norm: Vector is NULL!\n");
        return 0;
    }
    int incx = 1;
    return dnrm2_(&(pX->dim), pX->data, &incx);
}

void CopyDvector(Dvector *pX, Dvector *pY)
    /*X := Y*/
{
    if(pX == NULL || pY == NULL)
    {
        printf("CopyDvector: Vector(s) is NULL! Calloc vectors first!");
        return;
    }
    if(pX->dim != pY->dim)
    {
        printf("CopyDvector: warning, the size of two vector is different!\n"); 
    }
    int inc = 1;
    dcopy_(&(pY->dim), pY->data, &inc, pX->data, &inc);
}

void RandomDvector(Dvector *pX)
    /*generate random double vector
     *all entries have normal distrubution between (-1,1)*/
{
    int seed[4];
    time_t t;
    srand((unsigned) time(&t));
    int i;
    for(i=0; i< 4; i++){
        if(i!=3)
            seed[i] = rand() % 4095;
        else
            seed[i] = (rand() % 4095)/2+1;
    }
    int idist=2;
    dlarnv_(&idist, seed, &(pX->dim), pX->data);
}


void Lapack_Dgemv(Dvector *pY, Dmatrix *pA, Dvector *pX)
    /* Y = A*X */
{
    double beta = 0;
    double alpha = 1;
    int inc = 1;
    dgemv_("N", &(pA->nDimRow), &(pA->nDimCol), &alpha,
            pA->data, &(pA->nDimRow), pX->data, &inc, &beta,
            pY->data, &inc);
} 
