#include "LF_function.h"


void bustype(int *PV, int *PQ, int *type_bus)
{
	int j = 0;
	int k = 0;
	for (int i = 0; i < num_buses; i++)
	{
		if (type_bus[i] == 2 || type_bus[i] == 1)
		{
			PV[j] = i;
			j++;
		}
		if (type_bus[i] == 3)
		{
			PQ[k] = i;
			k++;
		}
	}
}


void initializevector( double *vect, int n)
{
	for (int i=0;i<n;i++)
	{
		vect[i]=0.0f;
	}

}

void printvector(double *vect,int n)
{
	for(int i=0;i<n;i++)
	{
		cout<< vect[i]<<endl;
	}


}


void initializematrix( double **mat, int nrows, int ncols)
{
	for (int i=0;i<nrows;i++)
	{
		for(int j=0;j<ncols;j++)
		{
			mat[i][j]=0.0f;
		}
	}

}


void printmatrix( double **mat, int nrows, int ncols)
{
	for (int i=0;i<nrows;i++)
	{
		for(int j=0;j<ncols;j++)
		{
			cout<<mat[i][j]<<"\t";
		}
		cout<<endl;
	}

}

void Jacobian(double **J, double **J11, double **J12, double **J21, double **J22, int n)
{
	for(int i=0;i<n-1;i++)
	{
		for (int j=0;j<n-1;j++)
		{
			J[i][j] = J11[i][j];
		}
	}
	for(int i=0;i<n-1;i++)
	{
		for (int j=n-1;j<22;j++)
		{
			J[i][j] = J12[i][j-13];
		}
	}
	for(int i=n-1;i<22;i++)
	{
		for(int j=0;j<n-1;j++)
		{
			J[i][j] = J21[i-13][j];
		}

	}
	for(int i=n-1;i<22;i++)
	{
		for(int j=n-1;j<22;j++)
		{
			J[i][j] = J22[i-13][j-13];
		}

	}

}


void InverseGE(double **Jinv, double **J, int n)
{
	double **temp = new double *[44];
	double d;
	for (int tt=0;tt<44;tt++)
	{
		temp[tt]= new double[44];
	}
	initializematrix(temp,44,44);
	//Assigning the matrix

	for (int i=0;i<n;i++)
	{
		for (int j=0;j<n;j++)
		{
			temp[i][j]=J[i][j];
		}
	}
	// Gauss Elimination method	
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < 2 * n; j++)
		{
			if (j == (i + n))
			{
				temp[i][j] = 1;
			}
		}
	}
	for (int i = n-1; i > 0; i--)
	{
		if (temp[i - 1][0] < temp[i][0])
			for (int j = 0; j <= n * 2; j++)
			{
				d = temp[i][j];
				temp[i][j] = temp[i - 1][j];
				temp[i - 1][j] = d;
			}
	}
	/********** reducing to diagonal  matrix ***********/

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n * 2; j++)
			if (j != i)
			{
				d = temp[j][i] / temp[i][i];
				for (int k = 0; k < n * 2; k++)
					temp[j][k] -= temp[i][k] * d;
			}

	}


	for (int i = 0; i < n; i++)
	{
		d = temp[i][i];
		for (int j = 0; j < n * 2; j++)
			temp[i][j] = temp[i][j] / d;
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = n ; j < n * 2; j++)
		{
			//cout << temp[i][j] << "    ";
			Jinv[i][j-n]=temp[i][j];
		}
		//cout << endl;
	}




	for(int tt=0;tt<44;tt++)
	{
		delete[] temp[tt];
	}

	delete[] temp;
}

void cg(double *X1,double *M,double **J,double *busvoltage,double *busangle,int *type_bus)
{
	double eps = 1e-12;
	double *r = new double[num_buses+8];
	double alpha = 0.0f;
	double num = 0.0f;
	double *temp2 = new double[num_buses+8];
	double *p = new double[num_buses+8];
	double *x = new double[num_buses+8];
	double *temp = new double[num_buses+8];
	double den = 0.0f;
	double tol = 0.0f;
	double lambda;
	int kk;
	double *tbestimated = new double[22];
	int iter = 10000;
	initializevector(tbestimated,22);
	for(int i =1;i<num_buses;i++)
	{
		tbestimated[i-1] = busangle[i];
	}
	kk=1;
	for (int j=1;j<num_buses;j++)
	{
		if(type_bus[j]==3)
		{
			tbestimated[kk+12]=busvoltage[j];
			kk++;
		}
	}
	//printvector(tbestimated,22);
	matvec(temp, J, tbestimated,22);
	//printvector(temp,22);
#pragma acc parallel
	for (int i = 0; i < num_buses+8; i++)
	{
		r[i] = M[i] - temp[i];
		p[i] = r[i];
		x[i] = tbestimated[i];
	}

	for (int k = 0; k < iter; k++)
	{
		num = 0.0f;
		den = 0.0f;
		matvec(temp2, J, p,22);
#pragma acc parallel
		for (int i = 0; i < num_buses+8; i++)
		{
			num += p[i] * r[i];
			den += p[i] * temp2[i];
		}
		alpha = num / den;
#pragma acc parallel
		for (int i = 0; i < num_buses+8; i++)
		{
			X1[i] = x[i] + alpha*p[i];

		}
		matvec(temp, J, X1,22);
		tol = 0.0f;
#pragma acc parallel
		for (int i = 0; i < num_buses+8; i++)
		{
			r[i] = 0.0f;
			r[i] = M[i] - temp[i];
			tol += r[i] * r[i];
		}
		tol = abs(tol);
		if (tol < eps)
		{
			cout << "Number of iterations to converge for CG is  " << k << endl;
			break;

		}
		den = 0.0f;
		num = 0.0f;
		matvec(temp2, J, p,22);
#pragma acc parallel
		for (int i = 0; i < num_buses+8; i++)
		{
			num -= r[i] * temp2[i];
			den += p[i] * temp2[i];
		}
		lambda = num / den;
#pragma acc parallel
		for (int i = 0; i < num_buses+8; i++)
		{
			p[i] = r[i] + lambda*p[i];
			x[i] = X1[i];
		}
	//	printvector(X1,22);
	//	cout<<k<<endl;	
	}
	delete[] r;
	delete[] p;
	delete[] temp;
	delete[] temp2;
	delete[] x;
	delete[] tbestimated;

}


void matvec(double *outvec, double **inmat, double *invec,int n)
{
	for (int i = 0; i < n; i++)
	{
		outvec[i] = 0.0f;
		for (int j = 0; j < n; j++)
		{
			outvec[i] += inmat[i][j] * invec[j];
		}
	}
}


void pol2deg(double *out,double *in,int n)
{
 for(int i =0;i<n;i++)
	{
	out[i] = (in[i]*180)/pi;
	}

}

