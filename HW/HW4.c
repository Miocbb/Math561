#include<stdio.h>
#include<stdlib.h>
#include"matrix.h"
#include"vector.h"
#include"QR.h"
#include"least_square.h"

void ReadCSV( Dmatrix *pMatrix, char * file);
void ReadCSV_Col( Dvector *pVector, char * file, int Col);

int main()
{
    /*******************************
     * HW4_a: least square problem
     * AX ~ b 
     *******************************/
    printf("HW4: (a)\n");
    Dmatrix *pA;
    Dvector *pResi, *pAX;
    Dvector *pb, *pX;

    pA = CreateDmatrix(4898,11);
    pb = CreateDvector(4898);
    pX = CreateDvector(11);
    pResi = CreateDvector(4898);
    pAX   = CreateDvector(4898);
    
    /*load matrix A and vector b*/
    ReadCSV(pA, "datasets/HW4.csv");
    ReadCSV_Col(pb, "datasets/HW4.csv", 12);

    LeastSquare(pA, pX, pb);
    printf("Coefficient vector Beta:\n");
    ShowDvector(pX);

    /*show Residue*/
    Lapack_Dgemv(pAX, pA, pX);
    DvectorArithmetic(pResi, 1, pAX, -1, pb ); /*Resi = AX - b*/
    printf("Residue: %f\n", Dvector_2Norm(pResi));
    DeleteDvector(&pX);

    /*************************************** 
     * HW4_b: modified least square problem
     * AX+C ~ b
     *****************************************/
    printf("\nHW4: (b)\n");
    Dvector *pC;
    pC = CreateDvector(4898);
    pX = CreateDvector(12);
    for(int i=0; i < pC->dim; i++)
        pC->data[i]  = 1;

    DmatrixExpansionByColVect(pA, pC); /* A := A(left) + C(right)*/
    LeastSquare(pA, pX, pb);
    printf("Coefficient vector Beta:\n");
    ShowDvector(pX);
    printf("Constant vector C: %f\n", pX->data[11]);

    /*Show Residue*/
    Lapack_Dgemv(pAX, pA, pX);
    DvectorArithmetic(pResi, 1, pAX, -1, pb ); /* Resi = AX -b */
    printf("Residue: %f\n", Dvector_2Norm(pResi));

    /*clean up*/
    DeleteDmatrixList(1, &pA);
    DeleteDvectorList(2, &pb, &pX);
}


void ReadCSV( Dmatrix *pMatrix, char * file)
    /* return: a Dmatrix*/
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
