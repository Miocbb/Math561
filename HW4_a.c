#include<stdio.h>
#include<stdlib.h>
#include"matrix.h"
#include"vector.h"
#include"QR.h"


void ReadCSV( Dmatrix *pMatrix, char * file);
void ReadCSV_Col( Dvector *pVector, char * file, int Col);

int main()
{
    /* AX ~ b 
     * least square problem
     */
    Dmatrix *pA, *pQ_A, *pR_A;
    Dmatrix *pR_Inv_A, *pQ_T_A;
    Dvector *pb, *pX, *pQb;

    pA = CreateDmatrix(4898,11);
    pb = CreateDvector(4898);
    pX = CreateDvector(11);
    pQ_A = CreateDmatrix(4898,11);
    pR_A = CreateDmatrix(11,11);
    pR_Inv_A = CreateDmatrix(11,11);
    pQ_T_A = CreateDmatrix(11, 4898);
    pQb = CreateDvector(11);
    ReadCSV(pA, "datasets/HW4.csv");
    ReadCSV_Col(pb, "datasets/HW4.csv", 12);

    /* do modified GS_QR */
    MGS_QR(pA, pQ_A, pR_A);

    //Lapack_Dtrtri(pR_Inv_A, pR_A); /*get the inverse of R*/
    TransposeDmatrix(pQ_T_A, pQ_A); /*get transpose of Q*/
    Lapack_Dgemv(pQb, pQ_T_A, pb); /* Qb = Q*T * b */
    //Lapack_Dgemv(pX, pR_Inv_A, pQb); /* X = R^-1 * Qb */

    
    SolveUpTriMatrixEq(pR_A, pX, pQb);
    printf("coefficient vector Beta:\n");
    SolveUpTriMatrixEq(pR_A, pX, pQb);
    ShowDvector(pX);

    Dvector *pResi, *pAX;
    pResi = CreateDvector(4898);
    pAX  = CreateDvector(4898);
    Lapack_Dgemv(pAX, pA, pX);
    DvectorArithmetic(pResi, 1, pAX, -1, pb );
    printf("Residue: %f\n", Dvector_2Norm(pResi));

    DeleteDmatrixList( 5, &pA, &pQ_A, &pR_A, &pR_Inv_A, &pQ_T_A);
    DeleteDvectorList( 5, &pb, &pX, &pQb, &pResi, &pAX);
}


void ReadCSV( Dmatrix *pMatrix, char * file)
    /* return: a Dmatrix
    */
{
    FILE *fp;
    int i = 0;
    double trash;
    char delimiter;
    Dmatrix *pMatrix_tem;
    fp = fopen(file, "r");
    
    if(fp == NULL)
        printf("REadCSV: open CSV file failed!\n");
    
    pMatrix_tem = CreateDmatrix(pMatrix->nDimCol, pMatrix->nDimRow);
    while(1)
    {
        int col = 1;
        while(1) 
        {
            if(col < 12 )
            {
                fscanf(fp, "%lf", &(pMatrix_tem->data[i]) );
                i++;
            }
            else{
                fscanf(fp, "%lf", &(trash) );
            }
            col++;
            delimiter = fgetc(fp);
            if(delimiter == '\n' || delimiter == EOF)
                break;
        }
        if(delimiter == EOF)
            break;
    }
    fclose(fp);
    TransposeDmatrix(pMatrix, pMatrix_tem);
    DeleteDmatrix(&(pMatrix_tem));
}

void ReadCSV_Col( Dvector *pVector, char * file, int Col)
    /* return: a Dvector
     *read the specified column of a CSV file
     * */

{
    FILE *fp;
    int i = 0;
    double trash;
    char delimiter;
    fp = fopen(file, "r");
    if(fp == NULL)
        printf("REadCSV: open CSV file failed!\n");
    
    while(1)
    {
        int col = 1;
        while(1) 
        {   
            if(col == Col)
            {
                fscanf(fp, "%lf", &(pVector->data[i]) );
                i++;
            }
            else{
                fscanf(fp, "%lf", &(trash) );
            }
            col++;
            delimiter = fgetc(fp);
            if(delimiter == '\n' || delimiter == EOF)
                break;
        }
        if(delimiter == EOF)
            break;
    }
    fclose(fp);
}
