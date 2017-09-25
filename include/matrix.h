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

/*
DoubleMatrix *CreateDoubleMatrix(int nDimRow, int nDimCol);
void InitDoubleMatrix(DoubleMatrix *pMatrix, double *pArray);
void InitDoubleMatrixRow(DoubleMatrix *pMatrix, double *pArray, int Row);
void InitDoubleMatrixCol(DoubleMatrix *pMatrix, double *pArray, int Col);
void DeleteDoubleMatrix(DoubleMatrix *pMatrix);
void ShowDoubleMatrix(DoubleMatrix *pMatrix);

DoubleVector *ExtractMatrixRow(DoubleMatrix *pMatrix, int Row);
DoubleVector *ExtractMatrixCol(DoubleMatrix *pMatrix, int Col);
DoubleMatrix *MatrixTranspose(DoubleMatrix *pMatrix);
DoubleMatrix *MatrixExpansionByCol(DoubleMatrix *pMatrix1, DoubleMatrix *pMatrix2);
DoubleMatrix *MatrixExpansionByRow(DoubleMatrix *pMatrix1, DoubleMatrix *pMatrix2);

DoubleMatrix *MatrixMultiplication(DoubleMatrix *pMatrix_1, DoubleMatrix *pMatrix_2);
*/
#endif
