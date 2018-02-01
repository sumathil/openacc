#include "LF_function.h"

int main()
{
	
	int *fb = new int[linefile];
	int *tb = new int[linefile];
	double *R = new double[linefile];
	double *X = new double[linefile];
	double *B_shunt = new double[linefile];
	double *tap = new double[linefile];
	int *bus_name = new int[num_buses];
	int *type_bus = new int[num_buses];
	int *PV = new int[5];
	int *PQ = new int[9];
	int kk,mm,nn;
	double *busvoltage = new double[num_buses];
	double *busangle = new double[num_buses];
	double *Real_powerg = new double[num_buses];
	double *Reactive_powerg = new double[num_buses];
	double *Real_powerl = new double[num_buses];
	double *Reactive_powerl = new double[num_buses];
	double *Qmax = new double[num_buses];
	double *Qmin = new double[num_buses];
	double *P = new double[num_buses];
	double *Q = new double[num_buses];
	double *Pspe = new double[num_buses];
	double *Qspe = new double[num_buses];
	double *dP = new double[num_buses];
	double *dQ = new double[num_buses];
	double *dPm = new double[13];
	double *dQm = new double[9];
	double *M = new double[22];
	double **J11 = new double *[num_buses-1];
	double **J12 = new double *[num_buses-1];
	double **J21 = new double *[9];
	double **J22 = new double *[9];
	double **J = new double *[22];
	double **JinvG = new double *[22];
	double *X1 = new double[22];
	double *dth = new double[num_buses-1];
	double *dV = new double[9];
	double epsilon = 1e-4;
	double Qg,Tol;
	for (int i=0;i<22;i++)
	{
		JinvG[i] = new double[22];
	}

	for (int i=0;i<22;i++)
	{
		J[i] = new double[22];
	}

	for (int i=0;i<num_buses-1;i++)
	{
		J11[i] = new double[num_buses-1];
	}
	for (int i=0;i<num_buses-1;i++)
	{
		J12[i] = new double[9];
	}

	for (int i=0;i<9;i++)
	{
		J21[i] = new double[num_buses-1];
	}

	for (int i=0;i<9;i++)
	{
		J22[i] = new double[9];
	}

	rect *y_line = new rect;
	rect_matrix *Y = new rect_matrix;
	ifstream linedata("linedata.txt");
	for (int i = 0; i < linefile; i++)
	{
		if (linedata.eof())
		{
			cout << "end of file reached" << endl;
			break;
		}
		linedata >> fb[i] >> tb[i] >> R[i] >> X[i] >> B_shunt[i] >> tap[i];
	}

	ifstream busdata("busdata.txt");
	for (int i = 0; i < num_buses; i++)
	{
		if (busdata.eof())
		{
			cout << "end of file reached" << endl;
			break;
		}
		busdata >> bus_name[i] >> type_bus[i] >> busvoltage[i] >> busangle[i] >> Real_powerg[i] >> Reactive_powerg[i] >> Real_powerl[i] >> Reactive_powerl[i] >> Qmin[i] >> Qmax[i];
		Real_powerg[i] = Real_powerg[i] / 100;
		Real_powerl[i] = Real_powerl[i] / 100;
		Reactive_powerg[i] = Reactive_powerg[i] / 100;
		Reactive_powerl[i] = Reactive_powerl[i] / 100;
		Qmin[i] = Qmin[i] / 100;
		Qmax[i] = Qmax[i] / 100;
		Pspe[i] = Real_powerg[i] - Real_powerl[i];
		Qspe[i] = Reactive_powerg[i] - Reactive_powerl[i];

	}
		
	for (int i = 0; i < linefile; i++)
	{
		y_line->real[i] = R[i] / (R[i] * R[i] + X[i] * X[i]);
		y_line->imag[i] = -X[i]/ (R[i] * R[i] + X[i] * X[i]);
	}


	for (int i = 0; i < num_buses; i++)
	{
		for (int j = 0; j < num_buses; j++)
		{
			Y->real[i][j] = 0.0f;
			Y->imag[i][j] = 0.0f;
		}
	}
	initializevector(dP,num_buses);
	initializevector(dQ,num_buses);
	for (int k = 0; k < linefile; k++)
	{
		Y->real[fb[k]-1][tb[k]-1] = Y->real[fb[k]-1][tb[k]-1] - y_line->real[k] / tap[k];
		Y->imag[fb[k]-1][tb[k]-1] = Y->imag[fb[k]-1][tb[k]-1] - y_line->imag[k] / tap[k];
		Y->real[tb[k]-1][fb[k]-1] = Y->real[fb[k]-1][tb[k]-1];
		Y->imag[tb[k]-1][fb[k]-1] = Y->imag[fb[k]-1][tb[k]-1];
	}

	for (int m = 0; m < num_buses; m++)
	{
		for (int n = 0; n < linefile; n++)
		{
			if ((fb[n]-1) == m)
			{
				Y->real[m][m] = Y->real[m][m] + y_line->real[n] / (tap[n]*tap[n]);
				Y->imag[m][m] = Y->imag[m][m] + y_line->imag[n] / (tap[n]*tap[n])+B_shunt[n];
			}
			else if ((tb[n]-1) == m)
			{
				Y->real[m][m] = Y->real[m][m] + y_line->real[n];
				Y->imag[m][m] = Y->imag[m][m] + y_line->imag[n] + B_shunt[n];
			}
		}
	}

	/*for (int i = 0; i < num_buses; i++)
	  {
	  for (int j = 0; j < num_buses; j++)
	  {
	  cout << Y->real[i][j] << "+" << Y->imag[i][j] << '\t';
	  }
	  cout << endl;

	  }*/

	bustype(PV,PQ,type_bus);


	/*for (int i = 0; i< 5; i++)
	  {
	  cout << "PV buses are " << PV[i] << endl;
	  }*/
	/*for (int i = 0; i< 9; i++)

	  cout << "PQ buses are " << PQ[i] << endl;
	  }*/

	// Calculating the Power flow values based on the conditions


//while(Tol < eps)

	for(int it=1;it<100;it++)
	{
	for (int i = 0; i < num_buses; i++)
	{
		P[i]=0.0f;
		Q[i]=0.0f;
		for (int k = 0; k < num_buses; k++)
		{
			P[i] = P[i] + busvoltage[i] * busvoltage[k] * (Y->real[i][k] * cos(busangle[i] - busangle[k]) + Y->imag[i][k] * sin(busangle[i] - busangle[k]));
			Q[i] = Q[i] + busvoltage[i] * busvoltage[k] * (Y->real[i][k] * sin(busangle[i] - busangle[k]) - Y->imag[i][k] * cos(busangle[i] - busangle[k]));
		}
	}
/*cout<<"for iteration "<<it<< endl<<endl;;
cout<<endl;
printvector(Q,num_buses);*/


	for(int i=0;i<num_buses;i++)
	{
		//cout<<P[i]<<" = Real Power" << "\t" << Q[i] <<" =Reactive Power"<<endl;
	}

	// Checking Q-limit violations
	if(it <= 7 && it>2)
	{
	for (int ut=1;ut<num_buses;ut++)
	{
		if(type_bus[ut] == 2)
		{
		Qg = Q[ut]+Reactive_powerl[ut];
		//cout<< Qg <<"for the bus "<<i<<endl;
		
		if(Qg<Qmin[ut])
			{
			busvoltage[ut]=busvoltage[ut]+0.01;
			}
		else if(Qg>Qmax[ut])
			{
                        busvoltage[ut]=busvoltage[ut]-0.01;
                        }	
		}
	}
	
	}	
/*cout<<"for iteration "<<it<< endl<<endl;;
cout<<endl;
printvector(busvoltage,num_buses);*/	
	// Calculating the difference in the power
	
	for (int i =0;i<num_buses;i++)
	{
		dP[i] = Pspe[i]-P[i];
		dQ[i] = Qspe[i]-Q[i];
		
	//	cout<< dP[i] <<"\t"<< dQ[i]<<endl;
	}


	// Calculating the mismatch vector


	initializevector(dQm,9);
	kk=0;
	for (int i=0;i<num_buses;i++)
	{
		if(type_bus[i]==3)
		{
			dQm[kk]=dQ[i];
			kk++;
		}		
	}


	initializevector(dPm,13);
	for (int i=1;i<num_buses;i++)
	{
		dPm[i-1]=dP[i];
	}


//printvector(dQm,9);
//cout<<"Second mismatch"<<endl;
//printvector(dPm,13);


initializevector(M,22);
for(int i=0;i<num_buses-1;i++)
{
	M[i]= dPm[i];
}
for (int i=num_buses-1;i<22;i++)
{
	M[i]=dQm[i-(num_buses-1)];
}
//printvector(M,22);

//Calculate the Jacobian Matrix

initializematrix(J11,num_buses-1,num_buses-1);

for (int i=0;i<num_buses-1;i++)
{
	mm=i+1;
	for (int k=0;k<num_buses-1;k++)
	{	
		nn=k+1;
		if (nn==mm)
		{
			for(nn = 0;nn<num_buses;nn++)
			{
				J11[i][k] = J11[i][k] + busvoltage[mm]*busvoltage[nn] *(-Y->real[mm][nn]*sin(busangle[mm]-busangle[nn])+Y->imag[mm][nn]*cos(busangle[mm]-busangle[nn]));
			}
			J11[i][k] = J11[i][k] - busvoltage[mm]*busvoltage[mm]*Y->imag[mm][mm];
		}	
		else
		{
			J11[i][k] = busvoltage[mm]*busvoltage[nn] *(Y->real[mm][nn]*sin(busangle[mm]-busangle[nn])-Y->imag[mm][nn]*cos(busangle[mm]-busangle[nn]));
		}

	}	
}

initializematrix(J12,num_buses-1,9);
mm=0;
nn=0;
for (int i=0;i<num_buses-1;i++)
{
	mm = i+1;
	for (int k=0;k<9;k++)
	{
		nn = PQ[k];
		if (nn==mm)
		{
			for (nn=0;nn<num_buses;nn++)
			{
				J12[i][k] = J12[i][k] + busvoltage[mm]*busvoltage[nn] *(Y->real[mm][nn]*cos(busangle[mm]-busangle[nn])+Y->imag[mm][nn]*sin(busangle[mm]-busangle[nn]));
			}
			J12[i][k] = J12[i][k] + busvoltage[mm]*Y->real[mm][mm];
		}
		else
		{
			J12[i][k] = busvoltage[mm]*(Y->real[mm][nn]*cos(busangle[mm]-busangle[nn])+Y->imag[mm][nn]*sin(busangle[mm]-busangle[nn]));
		}

	}


}

initializematrix(J21,9,num_buses-1);
mm =0;
nn =0;
for (int i=0;i<9;i++)
{
	mm = PQ[i];
	for (int k=0;k<num_buses-1;k++)
	{
		nn = k+1;
		if (nn==mm)
		{
			for (nn=0;nn<num_buses;nn++)
			{
				J21[i][k] = J21[i][k] + busvoltage[mm]*busvoltage[nn]*(Y->real[mm][nn]*cos(busangle[mm]-busangle[nn])+Y->imag[mm][nn]*sin(busangle[mm]-busangle[nn]));
			}	
			J21[i][k] = J21[i][k] - busvoltage[mm]*busvoltage[mm]*Y->real[mm][mm];
		}
		else
		{
			J21[i][k] =  busvoltage[mm]*busvoltage[nn]*(-Y->real[mm][nn]*cos(busangle[mm]-busangle[nn])-Y->imag[mm][nn]*sin(busangle[mm]-busangle[nn]));
		}	
	}
}

initializematrix(J22,9,9);
mm=0;
nn=0;
for (int i=0;i<9;i++)
{
	mm=PQ[i];
	for(int k=0;k<9;k++)
	{
		nn =PQ[k];
		if(nn==mm)
		{
			for(nn = 0;nn<num_buses;nn++)
			{
				J22[i][k] = J22[i][k] + busvoltage[nn] *(Y->real[mm][nn]*sin(busangle[mm]-busangle[nn])-Y->imag[mm][nn]*cos(busangle[mm]-busangle[nn]));
			}
			J22[i][k] = J22[i][k] - busvoltage[mm]*Y->imag[mm][mm];
		}
		else
		{
			J22[i][k] = busvoltage[mm]*(Y->real[mm][nn]*sin(busangle[mm]-busangle[nn])-Y->imag[mm][nn]*cos(busangle[mm]-busangle[nn]));
		}
	}

}

//

/*cout<< "Jacobian Part1"<<endl;
  printmatrix(J11,num_buses-1,num_buses-1);
  cout<<endl;

  cout<< "Jacobian Part2"<<endl;
  printmatrix(J12,num_buses-1,9);
  cout<<endl;

  cout<< "Jacobian Part3"<<endl;
  printmatrix(J21,9,num_buses-1);
  cout<<endl;*/
/*
  cout<< "Jacobian Part4"<<endl;
  printmatrix(J22,9,9);
  cout<<endl;*/
//cout<<"For interation= "<<it<<endl<<endl<<endli;

initializematrix(J,22,22);
Jacobian(J,J11,J12,J21,J22,num_buses);
//printmatrix(J,22,22);


//cout<<"For interation= "<<it<<endl<<endl<<endl;
//Inverse through GE method

//initializematrix(JinvG,22,22);
//InverseGE(JinvG,J,22);
//printmatrix(JinvG,22,22);

initializevector(X1,22);
//matvec(X1,JinvG,M,22);
//printvector(X1,22);


// Conjugate Gradient Method


cg(X1, M, J, busvoltage,busangle,type_bus);


for(int i=0;i<num_buses-1;i++)
{
	dth[i] = X1[i];
}
for(int j=0;j<9;j++)
{
	dV[j] = X1[(num_buses-1)+j];
}
//cout<<"for iteration "<<it<< endl<<endl;
//printvector(dth,13);
//printvector(dV,9);
for(int i=0;i<num_buses;i++)
{
busangle[i+1] = dth[i]+busangle[i+1];
}
//printvector(busangle,num_buses);
kk=0;
for (int j=1;j<num_buses;j++)
{
if(type_bus[j] == 3)
{
	busvoltage[j]=dV[kk]+busvoltage[j];
	kk++;
} 
}

/*cout<<"for iteration "<<it<< endl<<endl;
cout<<endl;
printvector(busvoltage,num_buses);


printvector(busangle,num_buses);*/

Tol =abs(M[0]);
//cout<<Tol<<endl;
for (int i=1;i<22;i++)
{
if(Tol<abs(M[i]))
{
	Tol=abs(M[i]);
}
}
//cout<<Tol<<endl;
if (Tol < epsilon)
{
	cout<<"Voltage Vector"<<endl;
	printvector(busvoltage,num_buses);
	pol2deg(busangle,busangle,14);
	cout<<"Angle Vector"<<endl;
 	printvector(busangle,num_buses);
	cout << "Number of iterations to converge is " << it << endl;
	break;
	
}
//cout<<it<< endl;

}


delete[] fb;
delete[] tb;
delete[] R;
delete[] X;
delete[] B_shunt;             
delete[] tap;
delete[] y_line;
delete[] Y;
delete[] bus_name;
delete[] type_bus;
delete[] busvoltage;
delete[] busangle;
delete[] Real_powerg;
delete[] Reactive_powerg;
delete[] Real_powerl;
delete[] Reactive_powerl;
delete[] Qmax;
delete[] Qmin;
delete[] P;
delete[] Q;
delete[] Pspe;
delete[] Qspe;
delete[] PV;
delete[] PQ;
delete[] dP;
delete[] dQ;
delete[] dPm;
delete[] dQm;
delete[] M;

for (int i=0;i<num_buses-1;i++)
{
	delete[] J11[i];
}
for (int i=0;i<num_buses-1;i++)
{
	delete[] J12[i];
}

for (int i=0;i<9;i++)
{
	delete[] J21[i];
}

for (int i=0;i<9;i++)
{
	delete[] J22[i];
}

delete[] J11;
delete[] J12;
delete[] J21;
delete[] J22;
for (int i=0;i<22;i++)
{
	delete[] J[i];
}
delete[] J;

for(int i=0;i<22;i++)
{
delete[] JinvG[i];
}
delete[] JinvG;
delete[] X1;
delete[] dth;
delete[] dV;


}	



