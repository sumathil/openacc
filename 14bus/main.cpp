#include <stdio.h>
#include <iostream>
#include <complex.h>
#include <iomanip>
#include <string>
#include <fstream>
#include <cmath>
using namespace std;

#define num_buses 14
#define linefile 22

typedef struct
{
	double real[linefile];
	double imag[linefile];
} rect;

typedef struct
{
	double real[num_buses][num_buses];
	double imag[num_buses][num_buses];
}rect_matrix;



