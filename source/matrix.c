#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<stdarg.h>
#include<time.h>
#include"lapack.h"
#include"matrix.h"
#include"vector.h"

Dmatrix *CreateDmatrix(int nDimRow, int nDimCol)
{
    Dmatrix *_pMatrix;
    _pMatrix = calloc(1, sizeof(Dmatrix));
    _pMatrix->nDimRow = nDimRow;
    _pMatrix->nDimCol = nDimCol;
    _pMatrix->data = calloc(1, sizeof(double)*nDimRow*nDimCol);
    return _pMatrix;
}

void DeleteDmatrix(Dmatrix **pMatrix)
{   
    if(pMatrix != NULL)
    {
        free((*pMatrix)->data);
        free(*pMatrix);
        *pMatrix = NULL;
    }
}

void DeleteDmatrixList(int count, Dmatrix **pMatrix, ...)
{
    DeleteDmatrix(pMatrix);

    va_list list;
    int i;
    va_start(list, pMatrix);
    for(i=0; i < count-1; i++)
        DeleteDmatrix( va_arg(list, Dmatrix**) );
    va_end(list);
}

void RandomDmatrix(Dmatrix *pMatrix)
    /*generate random double matrix
     *all entries have normal distrubution between (-1,1)*/
{
    int seed[4];
    time_t t;
    srand((unsigned) time(&t));
    int i;
    int size;
    size = pMatrix->nDimCol * pMatrix->nDimRow;
    for(i=0; i< 4; i++){
        if(i!=3)
            seed[i] = rand() % 4095;
        else
            seed[i] = (rand() % 4095)/2+1;
    }
    int idist=2;
    dlarnv_(&idist, seed, &(size), pMatrix->data);
}

void InitDmatrix(Dmatrix *pMatrix, double *pArray)
{
    if(pMatrix == NULL){
        printf("InitDmatrix: Initial Matrix failed, Matrix is NULL!\n");
        return;
    }
    int dim = pMatrix->nDimRow*pMatrix->nDimCol;
    int inc = 1;
    dcopy_(&dim, pArray, &inc, pMatrix->data, &inc);
}

void ShowDmatrix(Dmatrix *pMatrix)
{
    if(pMatrix == NULL)
    {
        printf("ShowDmatrix: Dmatrix is NULL!\n");
        return;
    }
    int i, j;
    for(i=0; i<pMatrix->nDimRow; i++)
    {
        for(j=0; j<pMatrix->nDimCol; j++) 
        {
            printf("%.3f ", pMatrix->data[j*(pMatrix->nDimRow) + i]);
        }
        printf("\n");
    }
}

void InitDmatrixRow(Dmatrix *pMatrix, double *pArray, int Row)
{
    /*Initial certain row of a matrix from an array*/ 
    int inc = 1;
    dcopy_(&(pMatrix->nDimCol), pArray, &inc,
            &(pMatrix->data[Row]), &(pMatrix->nDimRow));
}


void CopyDmatrix(Dmatrix *pA, Dmatrix *pB)
    /*A= B*/
{
    if(pA == NULL || pB == NULL)
    {
        printf("CopyDmatrix: Copy matrix failed! \n");
        return;
    }

    pA->nDimCol = pB->nDimCol;
    pA->nDimRow = pB->nDimRow;

    int inc = 1;
    int size = (pB->nDimCol) * (pB->nDimRow);
    dcopy_(&size, pB->data, &inc, pA->data, &inc);
}


void InitDmatrixCol(Dmatrix *pMatrix, double *pArray, int Col)
{
    /*Initial certain column of a matrix from an array*/ 
    int inc = 1;
    dcopy_(&(pMatrix->nDimRow), pArray, &inc,
            &(pMatrix->data[Col*pMatrix->nDimRow]), &(inc));
}

void ExtractDmatrixRow(Dvector *pVector, Dmatrix *pMatrix, int Row)
{
    if(pMatrix == NULL)
    {
        printf("ExtractDmatrixRow: Extracting Row vector failed, matrix is NULL!\n");
        return;
    }
    int inc = 1;
    dcopy_(&(pMatrix->nDimCol), &(pMatrix->data[Row]), &(pMatrix->nDimRow),
            pVector->data, &inc);
}

void ExtractDmatrixCol(Dvector *pVector, Dmatrix *pMatrix, int Col)
{
    if(pMatrix == NULL)
    {
        printf("ExtractDmatrixCol: Extracting Col vector failed, matrix is NULL!\n");
        return;
    }
    int inc = 1;
    dcopy_(&(pMatrix->nDimRow), &(pMatrix->data[Col * pMatrix->nDimRow]), &(inc),
            pVector->data, &inc);
}

void DmatrixMulti_NN(Dmatrix *pA, Dmatrix *pB, Dmatrix *pC)
    /* A := B * C;*/
{
    if(pA == NULL || pB == NULL || pC == NULL)
    {
        printf("DmatrixMulti_NN: matrix is null, terminated matriix multi!\n"); 
        return;
    }
    if(pA->nDimCol != pC->nDimCol || pA->nDimRow != pB->nDimRow 
            || pB->nDimCol != pC->nDimRow)
    {
        printf("DmatrixMulti_NN: size not matched, terminated matrix multi!\n");
        return;
    }
    double alpha=1, beta=0;
        dgemm_("N", "N", &(pB->nDimRow), &(pC->nDimCol), &(pB->nDimCol),
                &alpha, pB->data, &(pB->nDimRow), pC->data, &(pC->nDimRow),
                    &(beta), pA->data, &(pA->nDimRow));
        return;
}

void DmatrixMulti_TN(Dmatrix *pA, Dmatrix *pB, Dmatrix *pC)
    /* A := (B)^T * C;*/
{
    if(pA == NULL || pB == NULL || pC == NULL)
    {
        printf("DmatrixMulti_TN: matrix is null, terminated matriix multi!\n"); 
        return;
    }
    if(pA->nDimCol != pC->nDimCol || pA->nDimRow != pB->nDimCol
            || pB->nDimRow != pC->nDimRow)
    {
        printf("DmatrixMulti_TN: size not matched, terminated matrix multi!\n"); 
        return;
    }
    double alpha=1, beta=0;
        dgemm_("T", "N", &(pB->nDimCol), &(pC->nDimCol), &(pB->nDimRow),
                &alpha, pB->data, &(pB->nDimRow), pC->data, &(pC->nDimRow),
                    &(beta), pA->data, &(pA->nDimRow));
        return;
}

void DmatrixMulti_NT(Dmatrix *pA, Dmatrix *pB, Dmatrix *pC)
    /* A := B * (C)^T;*/
{
    if(pA == NULL || pB == NULL || pC == NULL)
    {
        printf("DmatrixMulti_TN: matrix is null, terminated matriix multi!\n"); 
        return;
    }
    if(pA->nDimCol != pC->nDimRow || pA->nDimRow != pB->nDimRow
            || pB->nDimCol != pC->nDimCol)
    {
        printf("DmatrixMulti_TN: size not matched, terminated matrix multi!\n"); 
        return;
    }
    double alpha=1, beta=0;
        dgemm_("N", "T", &(pB->nDimRow), &(pC->nDimRow), &(pB->nDimCol),
                &alpha, pB->data, &(pB->nDimRow), pC->data, &(pC->nDimRow),
                    &(beta), pA->data, &(pA->nDimRow));
        return;
}

void DmatrixMulti_TT(Dmatrix *pA, Dmatrix *pB, Dmatrix *pC)
    /* A := (B)^T * (C)^T;*/
{
    if(pA == NULL || pB == NULL || pC == NULL)
    {
        printf("DmatrixMulti_TN: matrix is null, terminated matriix multi!\n"); 
        return;
    }
    if(pA->nDimCol != pC->nDimRow || pA->nDimRow != pB->nDimCol
            || pB->nDimRow != pC->nDimCol)
    {
        printf("DmatrixMulti_TN: size not matched, terminated matrix multi!\n"); 
        return;
    }
    double alpha=1, beta=0;
        dgemm_("T", "T", &(pB->nDimCol), &(pC->nDimRow), &(pB->nDimRow),
                &alpha, pB->data, &(pB->nDimRow), pC->data, &(pC->nDimRow),
                    &(beta), pA->data, &(pA->nDimRow));
        return;
}

void TransposeDmatrix(Dmatrix *pB, Dmatrix *pA)
    /*B := A(T)*/
{
    if(pA == NULL || pB == NULL){
        printf("TransposeDmatrix: Matrix is NULL! Transpose input matrix fail!\n");
        return;
    }
    int i;
    Dvector *temp;
    Dmatrix *pMatrix_tem;

    pB->nDimCol=pA->nDimRow;
    pB->nDimRow=pA->nDimCol;

    temp = CreateDvector(pA->nDimCol);
    if(pA == pB){
        pMatrix_tem = CreateDmatrix(pA->nDimRow, pA->nDimCol);
        CopyDmatrix(pMatrix_tem, pA);
    }
    else
        pMatrix_tem = pA;

    for(i=0; i<pA->nDimRow; i++)
    {
        ExtractDmatrixRow(temp, pMatrix_tem, i);
        InitDmatrixCol(pB, temp->data, i);
    }
    DeleteDvector(&temp);
    if(pA == pB)
        DeleteDmatrix(&pMatrix_tem);
}

void DmatrixExpansionByCol(Dmatrix *pA, Dmatrix *pB)
/* A := A(left) + B(right); matrix expansion along column*/
{
    if(pA == NULL || pB== NULL)
    {
        printf("MatrixExpansionByCol: Matrix expansion by column failed! Matrix(s) is NULL!\n");
        return;
    }
    if(pA->nDimRow != pB->nDimRow)
    {
        printf("MatrixExpansionByCol: Matrixs do not share same size in row! Expansion failed!\n"); 
        return;
    }
    int dim = pB->nDimRow*pB->nDimCol;
    int inc = 1;
    pA->data = (double *)realloc(pA->data, (dim+pA->nDimCol*pA->nDimRow)*sizeof(double));
    dcopy_(&dim, pB->data, &inc, &(pA->data[pA->nDimCol*pA->nDimRow]), &inc);
    pA->nDimCol = pA->nDimCol+pB->nDimCol;
}

void DmatrixExpansionByRow(Dmatrix *pA, Dmatrix *pB)
/* A := A(up) + B(down); matrix expansion along row; */
{
    if(pA == NULL || pB== NULL)
    {
        printf("MatrixExpansionByRow: Matrix expansion by row failed! Matrix(s) is NULL!\n");
        return;
    }
    if(pA->nDimCol != pB->nDimCol)
    {
        printf("MatrixExpansionByRow: Matrixs do not share same size in Col! Expansion failed!\n"); 
        return;
    }
    int dim = (pA->nDimRow+pB->nDimRow)*(pA->nDimCol);
    int j;
    int inc = 1;
    Dmatrix *ptem;
    ptem = CreateDmatrix(pA->nDimRow+pB->nDimRow, pA->nDimCol);
    for(j=0; j<pA->nDimCol; j++)
    {
        dcopy_(&(pA->nDimRow), &(pA->data[j*pA->nDimRow]), &inc, &(ptem->data[j*ptem->nDimRow]), &inc);
        dcopy_(&(pB->nDimRow), &(pB->data[j*pB->nDimRow]), &inc, &(ptem->data[j*ptem->nDimRow+pA->nDimRow]), &inc);
    }
    pA->data = (double *)realloc(pA->data, dim * sizeof(double));
    dcopy_(&dim, ptem->data, &inc, pA->data, &inc);
    pA->nDimRow = pA->nDimRow+pB->nDimRow;
    DeleteDmatrix(&ptem);
}


void SolveUpTriMatrixEq(Dmatrix *pR, Dvector *pX, Dvector *pb)
    /* RX = b
     * Solve the linear matrix equation for X.
     * R is a full rank upper triangular matrix;
     * X, b is a double vector;
     * */
{
    int i, j;
    int dim;

    dim = pR->nDimCol;

    for(i=dim-1; i >= 0; i--)
    {
        for(j=dim-1; j> i; j--)
        {
            pb->data[i] -= pR->data[j*dim+i] * pX->data[j];
        }
        pX->data[i] = pb->data[i] / pR->data[i*dim+i];
   }
}
