/* the matrix is saved in 1-d array by column-major*/
#ifndef _MATRIX
#define _MATRIX

#include"vector.h"

typedef struct
{
    int nDimRow;
    int nDimCol;
    double *data;
}Dmatrix;

Dmatrix *CreateDmatrix(int DimRow, int DimCol);
void DeleteDmatrix(Dmatrix **pMatrix);
void DeleteDmatrixList(int count, Dmatrix **pMatrix, ...);

void RandomDmatrix(Dmatrix *pMatrix);
void InitDmatrix(Dmatrix *pMatrix, double *pArray);
void ShowDmatrix(Dmatrix *pMatrix);
void InitDmatrixRow(Dmatrix *pMatrix, double *pArray, int Row);
void InitDmatrixCol(Dmatrix *pMatrix, double *pArray, int Row);

void ExtractDmatrixRow(Dvector *pVector, Dmatrix *pMatrix, int Row);
void ExtractDmatrixCol(Dvector *pVector, Dmatrix *pMatrix, int Col);
void DmatrixMulti_NN(Dmatrix *pMatrixOut, Dmatrix *pMatrix1, Dmatrix *pMatrix2);
void DmatrixMulti_TN(Dmatrix *pMatrixOut, Dmatrix *pMatrix1, Dmatrix *pMatrix2);
void DmatrixMulti_NT(Dmatrix *pMatrixOut, Dmatrix *pMatrix1, Dmatrix *pMatrix2);
void DmatrixMulti_TT(Dmatrix *pMatrixOut, Dmatrix *pMatrix1, Dmatrix *pMatrix2);
void TransposeDmatrix(Dmatrix *pMatrixOut, Dmatrix *pMatrix);
void DmatrixExpansionByCol(Dmatrix *pMatrix1, Dmatrix *pMatrix2);
void DmatrixExpansionByRow(Dmatrix *pMatrix1, Dmatrix *pMatrix2);

void CBLAS_DSYMV(double *YVec, double *AMat, double *XVec, int N);
void CopyDmatrix(Dmatrix *pA, Dmatrix *pB);
void Lapack_Dtrtri(Dmatrix *pA, Dmatrix *pB);
void Lapack_Dgemv(Dvector *pY, Dmatrix *pA, Dvector *pX);

void SolveUpTriMatrixEq(Dmatrix *pR, Dvector *pX, Dvector *pb);
#endif
