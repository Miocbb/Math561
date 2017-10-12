#ifndef _VECTOR
#define _VECTOR

typedef struct
{
    int     dim;
    double *data;
}Dvector;

Dvector *CreateDvector(int dim);
void DeleteDvector(Dvector **pVect);
void DeleteDvectorList(int count, Dvector **pX, ...);

void RandomDvector(Dvector *pVect);
void InitDvector(Dvector *pVect, double *pArray);
void ShowDvector(Dvector *pVect);
void DvectorArithmetic(Dvector *pVectOut, double alpha, Dvector *pVect1,
    double beta, Dvector *pVect2);
void DvectorScalar(Dvector *pVectOut, double alpha, Dvector *pVect);

double DvectorIP(Dvector *pVect1, Dvector *pVect2);
double Dvector_2Norm(Dvector *pVect);

void CopyDvector(Dvector *pVectOut, Dvector *pVect);
#endif
