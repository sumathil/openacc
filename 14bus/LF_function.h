// Function declaration //

#pragma once

#include "LF_include.h"

// Function Declaration //

void bustype(int *PV,int *PQ,int *type_bus);
void initializevector(double *vect, int n);
void printvector(double *vect,int n);
void initializematrix( double **mat, int nrows, int ncols);
void printmatrix( double **mat, int nrows, int ncols);
void Jacobian(double **J, double **J11, double **J12, double **J21, double **J22, int n);
//void InverseGE(double **Jinv, double **J, int n);

void cg(double *X1,double *M,double **J,double *busvoltage,double *busangle,int *type_bus);

void matvec(double *outvec, double **inmat, double *invec, int n);
void pol2deg(double *out,double *in,int n);
