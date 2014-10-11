#ifndef __INDIVIDUAL_H_
#define __INDIVIDUAL_H_

#include "global.h"
#include "common.h"


#define PI  3.1415926535897932384626433832795

class CIndividual{
public:
	CIndividual();
	virtual ~CIndividual();

	vector <double> x_var;
	vector <double> y_obj;
	int    rank;
	double    crowd_distance;

	void   rnd_init();
	void   obj_eval();
	void   show_objective();
	void   show_variable();

    bool   operator<(const CIndividual &ind2);
	bool   operator<<(const CIndividual &ind2);
    bool   operator==(const CIndividual &ind2);
    void   operator=(const CIndividual &ind2);
};

CIndividual::CIndividual()
{
	x_var = vector<double>(nvar, 0);
    y_obj = vector<double>(nobj, 0);
	rank = max_int;
	crowd_distance = 0;
}

CIndividual::~CIndividual()
{
//	x_var.swap(vector<double>());
//	y_obj.swap(vector<double>());
}

void CIndividual::rnd_init()
{
    for(int n=0;n<nvar;n++)
        x_var[n] = lowBound[n] + rnd_uni(&rnd_uni_init)*(uppBound[n] - lowBound[n]);    

}

void CIndividual::obj_eval()
{
/*
	if(!strcmp("UF1", strTestInstance))  CEC09_F1(y_obj, x_var);
	if(!strcmp("UF2", strTestInstance))  CEC09_F2(y_obj, x_var);
	if(!strcmp("UF3", strTestInstance))  CEC09_F3(y_obj, x_var);
	if(!strcmp("UF4", strTestInstance))  CEC09_F4(y_obj, x_var);
	if(!strcmp("UF5", strTestInstance))  CEC09_F5(y_obj, x_var);
	if(!strcmp("UF6", strTestInstance))  CEC09_F6(y_obj, x_var);
	if(!strcmp("UF7", strTestInstance))  CEC09_F7(y_obj, x_var);
	if(!strcmp("UF8", strTestInstance))  CEC09_F8(y_obj, x_var);
	if(!strcmp("UF9", strTestInstance))  CEC09_F9(y_obj, x_var);
	if(!strcmp("UF10", strTestInstance)) CEC09_F10(y_obj, x_var);
*/
	int nx = x_var.size();
	vector<double> x=x_var;
	int k = position_parameters;

	if(!strcmp("UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		y_obj[0] = x[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
		return;
	}
	if(!strcmp("UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		y_obj[0] = x[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 1.0 - sqrt(x[0]) + 2.0 * sum2 / (double)count2;
		return;
	}
	if(!strcmp("UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		y_obj[0] = x[0]				+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] = 1.0 - sqrt(x[0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	}
	if(!strcmp("UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		y_obj[0] = x[0]				+ 2.0*sum1 / (double)count1;
		y_obj[1] = 1.0 - x[0]*x[0]	+ 2.0*sum2 / (double)count2;
		return;
	}
	if(!strcmp("UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
		y_obj[0] = x[0]	      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = 1.0 - x[0] + hj + 2.0*sum2 / (double)count2;
		return;
	}
	if(!strcmp("UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2,prod1, prod2, yj, hj, N, E,pj;
		
		sum1   = sum2   = 0.0;
		prod1  = prod2  = 1.0;
		count1 = count2 = 0;
		N = 2.0; E = 0.1;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		hj = 2*(0.5/N + E)*sin(2.0*N*PI*x[0]);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = x[0]	      + hj + 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] = 1.0 - x[0] + hj + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	}
	if(!strcmp("UF7", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x[0],0.2);
		y_obj[0] = yj	    + 2.0*sum1 / (double)count1;
		y_obj[1] = 1.0 - yj + 2.0*sum2 / (double)count2;
		return;
	}
	if(!strcmp("UF8", strTestInstance))
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj;
		
		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			if(j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += yj*yj;
				count2++;
			}
			else
			{
				sum3  += yj*yj;
				count3++;
			}
		}
		y_obj[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
		y_obj[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
		y_obj[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
		return;
	}
	if(!strcmp("UF9", strTestInstance))
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, E;
		
		E = 0.1;
		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			if(j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += yj*yj;
				count2++;
			}
			else
			{
				sum3  += yj*yj;
				count3++;
			}
		}
		yj = (0.5+E)*(1.0-4.0*(2.0*x[0]-1.0)*(2.0*x[0]-1.0));
		if(yj<0.0) yj = 0.0;
		y_obj[0] = 0.5*(yj + 2*x[0])*x[1]		+ 2.0*sum1 / (double)count1;
		y_obj[1] = 0.5*(yj - 2*x[0] + 2.0)*x[1] + 2.0*sum2 / (double)count2;
		y_obj[2] = 1.0 - x[1]                   + 2.0*sum3 / (double)count3;
		return;
	}
	if(!strcmp("UF10", strTestInstance))
	{
		unsigned int j, count1, count2, count3;
		double sum1, sum2, sum3, yj, hj;
		
		sum1   = sum2   = sum3   = 0.0;
		count1 = count2 = count3 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - 2.0*x[1]*sin(2.0*PI*x[0]+j*PI/nx);
			hj = 4.0*yj*yj - cos(8.0*PI*yj) + 1.0;
			if(j % 3 == 1) 
			{
				sum1  += hj;
				count1++;
			} 
			else if(j % 3 == 2) 
			{
				sum2  += hj;
				count2++;
			}
			else
			{
				sum3  += hj;
				count3++;
			}
		}
		y_obj[0] = cos(0.5*PI*x[0])*cos(0.5*PI*x[1]) + 2.0*sum1 / (double)count1;
		y_obj[1] = cos(0.5*PI*x[0])*sin(0.5*PI*x[1]) + 2.0*sum2 / (double)count2;
		y_obj[2] = sin(0.5*PI*x[0])                  + 2.0*sum3 / (double)count3;
		return;
	}

	if(!strcmp("SCH", strTestInstance))
	{
		y_obj[1] = (x_var[0] - 5) * (x_var[0] - 5);

		if (x_var[0] >= -5 && x_var[0] <= 1)	y_obj[0] = -x_var[0];
		else if (x_var[0] > 1 && x_var[0] < 3)	y_obj[0] = -2 + x_var[0];
		else if (x_var[0] >= 3 && x_var[0] <= 4)	y_obj[0] = 4 - x_var[0];
		else if (x_var[0] > 4 && x_var[0] <= 10)	y_obj[0] = -4 + x_var[0];
		return;
	}

	if(!strcmp("DEB", strTestInstance))
	{
		y_obj[0] = x_var[0];
		y_obj[1] = (1+10*x_var[1])*(1-(x_var[0]/(1+10*x_var[1]))*(x_var[0]/(1+10*x_var[1]))-x_var[0]*sin(8*PI*x_var[0])/(1+10*x_var[1]));
		return;
	}

	if(!strcmp("KUR", strTestInstance))
	{
		int i;
		double sum1=0,sum2=0;
		for (i=0; i < nvar - 1; i++) sum1 = sum1 -10 *exp(-0.2*sqrt(x_var[i] * x_var[i]+x_var[i+1] * x_var[i+1]));
		for (i=0; i < nvar; i++)	sum2 = sum2 + pow(fabs(x_var[i]),0.8)+5*sin(pow(x_var[i],3));
		y_obj[0] = sum1;
		y_obj[1] = sum2;	
		return;
	}

	if(!strcmp("ZDT1",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);

		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(y_obj[0]/g));
		return;
	}

	if(!strcmp("ZDT2",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - pow(y_obj[0]/g,2));
		return;
	}

	if(!strcmp("ZDT3",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);

		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(x_var[0]/g) - x_var[0]*sin(10*PI*x_var[0])/g);
		return;
	}

	if(!strcmp("ZDT4",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
		{
			double x = x_var[n];
			g+= x*x - 10*cos(4*PI*x);
		}
		g = 1 + 10*(nvar-1) + g;
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1- sqrt(y_obj[0]/g));
		return;
	}

	if(!strcmp("ZDT6",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n]/(nvar - 1);
		g = 1 + 9*pow(g,0.25) ;

		y_obj[0] = 1 - exp(-4*x_var[0])*pow(sin(6*PI*x_var[0]),6);
		y_obj[1] = g*(1- pow(y_obj[0]/g,2));
		return;
	}

	if(!strcmp("DTLZ1",strTestInstance))
	{
		double g = 0;
		for(int n=2; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);
		y_obj[0] = 0.5*(1 + g)*x_var[0]*x_var[1];
		y_obj[1] = 0.5*(1 + g)*x_var[0]*(1 - x_var[1]);
		y_obj[2] = 0.5*(1 + g)*(1 - x_var[0]);
		return;
	}

	if(!strcmp("DTLZ2",strTestInstance))
	{
		double g = 0;		
		for(int n=2; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}
		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[2] = (1 + g)*sin(x_var[0]*PI/2);
		return;
	}

	if(!strcmp("DTLZ3",strTestInstance))
	{
		double g = 0;
		for(int n=2; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);

		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[2] = (1 + g)*sin(x_var[0]*PI/2);
		return;
	}

	if(!strcmp("DTLZ4",strTestInstance))
	{
		double g = 0;
		double alpha = 100;
		for(int n=2; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}

		y_obj[0] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*cos(pow(x_var[1],alpha)*PI/2);
		y_obj[1] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*sin(pow(x_var[1],alpha)*PI/2);
		y_obj[2] = (1 + g)*sin(pow(x_var[0],alpha)*PI/2);
		return;
	}

	if(!strcmp("DTLZ6",strTestInstance))
	{
		int n;
		double g = 0,h=0;
		for(n=2; n<nvar;n++)	g = g + x_var[n];
		g = 1 + 9*g /(nvar + 1 - nobj);
		for (n =0; n < nobj-1; n++) h = h + x_var[n]/(1+g)*(1+sin(3*PI*x_var[n]));
		h = nobj - h;

		y_obj[0] = x_var[0];
		y_obj[1] = x_var[1];
		y_obj[2] = (1+g) * h;
		return;
	}

	if(!strcmp("R2_DTLZ2_M5", strTestInstance))
	{
		
	}
	if(!strcmp("R2_DTLZ3_M5", strTestInstance))
	{
		
	}
	if(!strcmp("WFG1_M5", strTestInstance))
	{
		
	}
	if(!strcmp("CF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6*PI*x[0] + j*PI/nvar);
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		y_obj[0] = x[0]				+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] = 1.0 - sqrt(x[0]) + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	}
		if(!strcmp("MOP1", strTestInstance))
	{//Schaffer	 	
		y_obj[0] = x[0]*x[0];
		y_obj[1] = (x[0]-2)*(x[0]-2);
		return;
	}
	if(!strcmp("MOP2", strTestInstance))
	{//Fonsenca and Fleming
		int j;
		double sum1 = 0, sum2 = 0;
		for (j = 0; j < nx; j++)
		{
			sum1 += (x[j] - 1/sqrt(nx*1.0)) * (x[j] - 1/sqrt(nx*1.0));
			sum2 += (x[j] + 1/sqrt(nx*1.0)) * (x[j] + 1/sqrt(nx*1.0));
		}
		y_obj[0] = 1 - exp(-sum1);
		y_obj[1] = 1 - exp(-sum2);
		return;
	}
	if(!strcmp("MOP3", strTestInstance))
	{//Poloni
		double A1 = 0.5*sin(1.0) -2*cos(1.0)+sin(2.0)-1.5*cos(2.0);
		double A2 = 1.5*sin(1.0) -cos(1.0)+2*sin(2.0)-0.5*cos(2.0);
		double B1 = 0.5*sin(x[0.0]) -2*cos(x[0])+sin(x[1])-1.5*cos(x[1]);
		double B2 = 1.5*sin(x[0]) -cos(x[0])+2*sin(x[1])-0.5*cos(x[1]);
		
		y_obj[0] = 1 + (A1 - B1)*(A1 - B1) + (A2 - B2) * (A2 - B2);
		y_obj[1] = (x[0] + 3)*(x[0] + 3)+ (x[1] + 1)*(x[1] + 1);
		return;
	}
	if(!strcmp("MOP4", strTestInstance))
	{//Kursawe	
		double sum1=0,sum2=0;
		int i;
		for (i=0; i < nvar - 1; i++) sum1 = sum1 -10 *exp(-0.2*sqrt(x_var[i] * x_var[i]+x_var[i+1] * x_var[i+1]));
		for (i=0; i < nvar; i++)	sum2 = sum2 + pow(fabs(x_var[i]),0.8)+5*sin(pow(x_var[i],3));
		y_obj[0] = sum1;
		y_obj[1] = sum2;	
		return;
	}
	if(!strcmp("MOP5", strTestInstance))
	{//Viennet		
		y_obj[0] = 0.5*(x[0]*x[0]+x[1]*x[1])+sin(x[0]*x[0]+x[1]*x[1]);
		y_obj[1] = (3*x[0] - 2*x[1]+4)*(3*x[0] - 2*x[1]+4)/8+(x[0]-x[1]+1)*(x[0]-x[1]+1)/27+15;
		y_obj[2] = 1/(x[0]*x[0]+x[1]*x[1]+1)-1.1*exp(-x[0]*x[0] -x[1]*x[1]);
		return;
	}
	if(!strcmp("MOP6", strTestInstance))
	{//Deb	
		y_obj[0] = x_var[0];
		y_obj[1] = (1+10*x_var[1])*(1-(x_var[0]/(1+10*x_var[1]))*(x_var[0]/(1+10*x_var[1]))-x_var[0]*sin(8*PI*x_var[0])/(1+10*x_var[1]));
		return;
	}
	if(!strcmp("MOP7", strTestInstance))
	{//Viennet			
		y_obj[0] = (x[0]-2)*(x[0]-2)/2 + (x[1]+1)*(x[1]+1)/13 + 3;
		y_obj[1] = (x[0]+x[1]-3)*(x[0]+x[1]-3)/36 + (-x[0]+x[1]+2)*(-x[0]+x[1]+2)/8 -17;
		y_obj[2] = (x[0]+2*x[1]-1)*(x[0]+2*x[1]-1)/175 + (-x[0]+2*x[1])*(-x[0]+2*x[1])/17 -13;
		return;
	}
	if(!strcmp("WFG1", strTestInstance))
	{
		vector<double> y = WFG_normalise_z( x_var );
		
		y = WFG1_t1(y, k);
		y = WFG1_t2(y, k);
		//y = WFG1_t3(y);
		y = WFG1_t4(y, k);
		
		y_obj = WFG1_shape(y);

		return;
	}
	if(!strcmp("WFG2", strTestInstance))
	{
		vector< double > y = WFG_normalise_z( x_var );
		
		y = WFG1_t1( y, k);
		y = WFG2_t2( y, k);
		y = WFG2_t3( y, k);
		
		y_obj = WFG2_shape( y );

		return;
	}
	if(!strcmp("WFG3", strTestInstance))
	{
		vector< double > y = WFG_normalise_z( x_var );
		
		y = WFG1_t1( y, k);
		y = WFG2_t2( y, k);
		y = WFG2_t3( y, k);
		
		y_obj = WFG3_shape( y );

		return;
	}
	if(!strcmp("WFG4", strTestInstance))
	{
		vector< double > y = WFG_normalise_z( x_var );
		
		y = WFG4_t1( y);
		y = WFG2_t3( y, k);		
		
		y_obj = WFG4_shape( y );

		return;
	}
	if(!strcmp("WFG5", strTestInstance))
	{
		vector< double > y = WFG_normalise_z( x_var );

		y = WFG5_t1( y );
		y = WFG2_t3( y, k);
		
		y_obj = WFG4_shape( y );
		return;
	}
	if(!strcmp("WFG6", strTestInstance))
	{
		vector< double > y = WFG_normalise_z( x_var );

		y = WFG1_t1( y, k);
		y = WFG6_t2( y, k);
		
		y_obj = WFG4_shape( y );
		return;
	}
	if(!strcmp("WFG7", strTestInstance))
	{
		vector< double > y = WFG_normalise_z( x_var );

		y = WFG7_t1( y, k );
		y = WFG1_t1( y, k);
		y = WFG2_t3( y, k);
		
		y_obj = WFG4_shape( y );
		return;
	}
	if(!strcmp("WFG8", strTestInstance))
	{
		vector< double > y = WFG_normalise_z( x_var );

		y = WFG8_t1( y, k );
		y = WFG1_t1( y, k);
		y = WFG2_t3( y, k);
		
		y_obj = WFG4_shape( y );
		return;
	}
	if(!strcmp("WFG9", strTestInstance))
	{
		vector< double > y = WFG_normalise_z( x_var );

		y = WFG9_t1( y);
		y = WFG9_t2( y, k);
		y = WFG6_t2( y, k);
		
		y_obj = WFG4_shape( y );
		return;
	} 
/*	if(!strcmp("SCH_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]) + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}

		y_obj[1] = (x_var[0] - 5) * (x_var[0] - 5);

		if (x_var[0] >= -5 && x_var[0] <= 1)	y_obj[0] = -x_var[0];
		else if (x_var[0] > 1 && x_var[0] < 3)	y_obj[0] = -2 + x_var[0];
		else if (x_var[0] >= 3 && x_var[0] <= 4)	y_obj[0] = 4 - x_var[0];
		else if (x_var[0] > 4 && x_var[0] <= 10)	y_obj[0] = -4 + x_var[0];

		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} */
	if(!strcmp("SCH_VAR_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]) + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}

		double K=0.16;	double A = 7.4;	double B = 1;	double D = 0;	double E = -B;
		double C = uppBound[0]-sqrt((B-D)/K);
		y_obj[1] = K*(x_var[0] - uppBound[0]) * (x_var[0] - uppBound[0]);

		if (x_var[0] >= lowBound[0] && x_var[0] <= A)	y_obj[0] = E+(B-E)/(A-lowBound[0])*(x_var[0]-lowBound[0]);
		else if (x_var[0] > A && x_var[0] <= C)	y_obj[0] = B+(B-D)/(A-C)* (x_var[0]-A);
		else if (x_var[0] > C && x_var[0] <= uppBound[0])	y_obj[0] = D+(D-B)/(C-uppBound[0])*( x_var[0]-C);

		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 
	if(!strcmp("SCH_VAR_UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		double K=0.16;	double A = 7.4;	double B = 1;	double D = 0;	double E = -B;
		double C = uppBound[0]-sqrt((B-D)/K);
		y_obj[1] = K*(x_var[0] - uppBound[0]) * (x_var[0] - uppBound[0]);

		if (x_var[0] >= lowBound[0] && x_var[0] <= A)	y_obj[0] = E+(B-E)/(A-lowBound[0])*(x_var[0]-lowBound[0]);
		else if (x_var[0] > A && x_var[0] <= C)	y_obj[0] = B+(B-D)/(A-C)* (x_var[0]-A);
		else if (x_var[0] > C && x_var[0] <= uppBound[0])	y_obj[0] = D+(D-B)/(C-uppBound[0])*( x_var[0]-C);

		y_obj[0] = y_obj[0]					+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;
		return;
	} 
	if(!strcmp("SCH_VAR_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		double K=0.16;	double A = 7.4;	double B = 1;	double D = 0;	double E = -B;
		double C = uppBound[0]-sqrt((B-D)/K);
		y_obj[1] = K*(x_var[0] - uppBound[0]) * (x_var[0] - uppBound[0]);

		if (x_var[0] >= lowBound[0] && x_var[0] <= A)	y_obj[0] = E+(B-E)/(A-lowBound[0])*(x_var[0]-lowBound[0]);
		else if (x_var[0] > A && x_var[0] <= C)	y_obj[0] = B+(B-D)/(A-C)* (x_var[0]-A);
		else if (x_var[0] > C && x_var[0] <= uppBound[0])	y_obj[0] = D+(D-B)/(C-uppBound[0])*( x_var[0]-C);
		y_obj[0] = y_obj[0]			+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] =  y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	} 
	if(!strcmp("SCH_VAR_UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}	
		double K=0.16;	double A = 7.4;	double B = 1;	double D = 0;	double E = -B;
		double C = uppBound[0]-sqrt((B-D)/K);
		y_obj[1] = K*(x_var[0] - uppBound[0]) * (x_var[0] - uppBound[0]);

		if (x_var[0] >= lowBound[0] && x_var[0] <= A)	y_obj[0] = E+(B-E)/(A-lowBound[0])*(x_var[0]-lowBound[0]);
		else if (x_var[0] > A && x_var[0] <= C)	y_obj[0] = B+(B-D)/(A-C)* (x_var[0]-A);
		else if (x_var[0] > C && x_var[0] <= uppBound[0])	y_obj[0] = D+(D-B)/(C-uppBound[0])*( x_var[0]-C);

		y_obj[0] = y_obj[0]				+ 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1]	+ 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("SCH_VAR_UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, Emxitong;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; Emxitong = 0.1;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		double K=0.16;	double A = 7.4;	double B = 1;	double D = 0;	double E = -B;
		double C = uppBound[0]-sqrt((B-D)/K);
		y_obj[1] = K*(x_var[0] - uppBound[0]) * (x_var[0] - uppBound[0]);

		if (x_var[0] >= lowBound[0] && x_var[0] <= A)	y_obj[0] = E+(B-E)/(A-lowBound[0])*(x_var[0]-lowBound[0]);
		else if (x_var[0] > A && x_var[0] <= C)	y_obj[0] = B+(B-D)/(A-C)* (x_var[0]-A);
		else if (x_var[0] > C && x_var[0] <= uppBound[0])	y_obj[0] = D+(D-B)/(C-uppBound[0])*( x_var[0]-C);

		hj = (0.5/N + Emxitong)*fabs(sin(2.0*N*PI*x[0]));
		y_obj[0] =  y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;		
		return;
	} 
	if(!strcmp("SCH_VAR_UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, Emxitong;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 2.0; Emxitong = 0.1;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		double K=0.16;	double A = 7.4;	double B = 1;	double D = 0;	double E = -B;
		double C = uppBound[0]-sqrt((B-D)/K);
		y_obj[1] = K*(x_var[0] - uppBound[0]) * (x_var[0] - uppBound[0]);

		if (x_var[0] >= lowBound[0] && x_var[0] <= A)	y_obj[0] = E+(B-E)/(A-lowBound[0])*(x_var[0]-lowBound[0]);
		else if (x_var[0] > A && x_var[0] <= C)	y_obj[0] = B+(B-D)/(A-C)* (x_var[0]-A);
		else if (x_var[0] > C && x_var[0] <= uppBound[0])	y_obj[0] = D+(D-B)/(C-uppBound[0])*( x_var[0]-C);

		hj = 2*(0.5/N + Emxitong)*sin(2.0*N*PI*x[0]);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("SCH_VAR_UF7", strTestInstance))
	{
	unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow((x[0]-lowBound[0])/(uppBound[0]-lowBound[0]),0.2)*(uppBound[0]-lowBound[0])+lowBound[0];
		double K=0.16;	double A = 7.4;	double B = 1;	double D = 0;	double E = -B;
		double C = uppBound[0]-sqrt((B-D)/K);
		y_obj[1] = K*(yj - uppBound[0]) * (yj - uppBound[0]);

		if (yj >= lowBound[0] && yj <= A)	y_obj[0] = E+(B-E)/(A-lowBound[0])*(yj-lowBound[0]);
		else if (yj > A && yj <= C)	y_obj[0] = B+(B-D)/(A-C)* (yj-A);
		else if (yj > C && yj <= uppBound[0])	y_obj[0] = D+(D-B)/(C-uppBound[0])*(yj-C);

		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 
	if(!strcmp("DEB_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]) + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		y_obj[0] = x_var[0];
		y_obj[1] = (1+10*x_var[1])*(1-(x_var[0]/(1+10*x_var[1]))*(x_var[0]/(1+10*x_var[1]))-x_var[0]*sin(8*PI*x_var[0])/(1+10*x_var[1]));
		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 	
	if(!strcmp("DEB_UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 3; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		y_obj[0] = x_var[0];
		y_obj[1] = (1+10*x_var[1])*(1-(x_var[0]/(1+10*x_var[1]))*(x_var[0]/(1+10*x_var[1]))-x_var[0]*sin(8*PI*x_var[0])/(1+10*x_var[1]));
		y_obj[0] = y_obj[0]					+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;
		return;
	} 
	if(!strcmp("DEB_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		y_obj[0] = x_var[0];
		y_obj[1] = (1+10*x_var[1])*(1-(x_var[0]/(1+10*x_var[1]))*(x_var[0]/(1+10*x_var[1]))-x_var[0]*sin(8*PI*x_var[0])/(1+10*x_var[1]));
		y_obj[0] = y_obj[0]			+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] =  y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	} 
	if(!strcmp("DEB_UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}	
		y_obj[0] = x_var[0];
		y_obj[1] = (1+10*x_var[1])*(1-(x_var[0]/(1+10*x_var[1]))*(x_var[0]/(1+10*x_var[1]))-x_var[0]*sin(8*PI*x_var[0])/(1+10*x_var[1]));
		y_obj[0] = y_obj[0]				+ 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1]	+ 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("DEB_UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		y_obj[0] = x_var[0];
		y_obj[1] = (1+10*x_var[1])*(1-(x_var[0]/(1+10*x_var[1]))*(x_var[0]/(1+10*x_var[1]))-x_var[0]*sin(8*PI*x_var[0])/(1+10*x_var[1]));

		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
		y_obj[0] =  y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;		
		return;
	} 
	if(!strcmp("DEB_UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 2.0; E = 0.1;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		y_obj[0] = x_var[0];
		y_obj[1] = (1+10*x_var[1])*(1-(x_var[0]/(1+10*x_var[1]))*(x_var[0]/(1+10*x_var[1]))-x_var[0]*sin(8*PI*x_var[0])/(1+10*x_var[1]));

		hj = 2*(0.5/N + E)*sin(2.0*N*PI*x[0]);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("DEB_UF7", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x_var[0],0.2);	
		y_obj[0] = x_var[0];
		y_obj[1] = (1+10*x_var[1])*(1-(x_var[0]/(1+10*x_var[1]))*(x_var[0]/(1+10*x_var[1]))-x_var[0]*sin(8*PI*x_var[0])/(1+10*x_var[1]));
		//y_obj[0] = yj;
		//y_obj[1] = (1+10*x_var[1])*(1-(yj/(1+10*x_var[1]))*(yj/(1+10*x_var[1]))-yj*sin(8*PI*yj)/(1+10*x_var[1]));


		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 
	if(!strcmp("KUR_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		double x1=(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]);
		for(j = 4; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x1 + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		int i;
		y_obj[0]   =0.0, y_obj[1]  = 0.0;
		for (i=0; i < 2; i++) y_obj[0] = y_obj[0] -10 *exp(-0.2*sqrt(x_var[i] * x_var[i]+x_var[i+1] * x_var[i+1]));
		for (i=0; i < 3; i++)	y_obj[1] = y_obj[1] + pow(fabs(x_var[i]),0.8)+5*sin(pow(x_var[i],3));
		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 
	if(!strcmp("KUR_UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		double x1=(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]);
		for(j = 4; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*x1*(x1*cos(4.0*(6.0*PI*x1+j*PI/nx))+2.0)*cos(6.0*PI*x1+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x1*(x1*cos(4.0*(6.0*PI*x1+j*PI/nx))+2.0)*sin(6.0*PI*x1+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		int i;
		y_obj[0]   =0.0, y_obj[1]  = 0.0;
		for (i=0; i < 2; i++) y_obj[0] = y_obj[0] -10 *exp(-0.2*sqrt(x_var[i] * x_var[i]+x_var[i+1] * x_var[i+1]));
		for (i=0; i < 3; i++)	y_obj[1] = y_obj[1] + pow(fabs(x_var[i]),0.8)+5*sin(pow(x_var[i],3));

		y_obj[0] = y_obj[0]					+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;
		return;
	} 
	if(!strcmp("KUR_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		double x1=(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]);
		for(j = 4; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x1,0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		int i;
		y_obj[0]   =0.0, y_obj[1]  = 0.0;
		for (i=0; i < 2; i++) y_obj[0] = y_obj[0] -10 *exp(-0.2*sqrt(x_var[i] * x_var[i]+x_var[i+1] * x_var[i+1]));
		for (i=0; i < 3; i++)	y_obj[1] = y_obj[1] + pow(fabs(x_var[i]),0.8)+5*sin(pow(x_var[i],3));

		y_obj[0] = y_obj[0]			+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] =  y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	} 
	if(!strcmp("KUR_UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		double x1=(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]);
		for(j = 4; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x1+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}	
		int i;
		y_obj[0]   =0.0, y_obj[1]  = 0.0;
		for (i=0; i < 2; i++) y_obj[0] = y_obj[0] -10 *exp(-0.2*sqrt(x_var[i] * x_var[i]+x_var[i+1] * x_var[i+1]));
		for (i=0; i < 3; i++)	y_obj[1] = y_obj[1] + pow(fabs(x_var[i]),0.8)+5*sin(pow(x_var[i],3));

		y_obj[0] = y_obj[0]				+ 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1]	+ 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("KUR_UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		double x1=(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]);
		for(j = 4; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x1);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		int i;
		y_obj[0]   =0.0, y_obj[1]  = 0.0;
		for (i=0; i < 2; i++) y_obj[0] = y_obj[0] -10 *exp(-0.2*sqrt(x_var[i] * x_var[i]+x_var[i+1] * x_var[i+1]));
		for (i=0; i < 3; i++)	y_obj[1] = y_obj[1] + pow(fabs(x_var[i]),0.8)+5*sin(pow(x_var[i],3));

		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x1));
		y_obj[0] =  y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;		
		return;
	} 
	if(!strcmp("KUR_UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 2.0; E = 0.1;
		double x1=(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]);
		for(j = 4; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x1+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		int i;
		y_obj[0]   =0.0, y_obj[1]  = 0.0;
		for (i=0; i < 2; i++) y_obj[0] = y_obj[0] -10 *exp(-0.2*sqrt(x_var[i] * x_var[i]+x_var[i+1] * x_var[i+1]));
		for (i=0; i < 3; i++)	y_obj[1] = y_obj[1] + pow(fabs(x_var[i]),0.8)+5*sin(pow(x_var[i],3));

		hj = 2*(0.5/N + E)*sin(2.0*N*PI*x1);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("KUR_UF7", strTestInstance))
	{
	unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		double x1=(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]);
		for(j = 4; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x1+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		//yj = pow((x[0]-lowBound[0])/(uppBound[0]-lowBound[0]),0.2)*(uppBound[0]-lowBound[0])+lowBound[0];
		int i;
		y_obj[0]   =0.0, y_obj[1]  = 0.0;
		for (i=0; i < 2; i++) y_obj[0] = y_obj[0] -10 *exp(-0.2*sqrt(x_var[i] * x_var[i]+x_var[i+1] * x_var[i+1]));
		for (i=0; i < 3; i++)	y_obj[1] = y_obj[1] + pow(fabs(x_var[i]),0.8)+5*sin(pow(x_var[i],3));


		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 
	/*if(!strcmp("KUR_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum11, sum21, prod1, prod2, yj, pj;
		
		sum11   = sum21   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 4; j <= nx; j++) 
		{
			yj = x[j-1]-pow((x[0]-lowBound[0])/(uppBound[0]-lowBound[0]),0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum21  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum11 += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		int i;
		double sum1   =0.0, sum2   = 0.0;
		for (i=0; i < 2; i++) sum1 = sum1 -10 *exp(-0.2*sqrt(x_var[i] * x_var[i]+x_var[i+1] * x_var[i+1]));
		for (i=0; i < 3; i++)	sum2 = sum2 + pow(fabs(x_var[i]),0.8)+5*sin(pow(x_var[i],3));
		y_obj[0] = sum1;
		y_obj[1] = sum2;
		y_obj[0] = y_obj[0]				+ 2.0*(4.0*sum11 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] = y_obj[1] + 2.0*(4.0*sum21 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	}*/
	if(!strcmp("ZDT3_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		
		y_obj[0] = x_var[0];
		y_obj[1] = 1-sqrt(x_var[0])-x_var[0]*sin(10*PI*x_var[0]);
		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 	
	if(!strcmp("ZDT3_UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		y_obj[0] = x_var[0];
		y_obj[1] = 1-sqrt(x_var[0])-x_var[0]*sin(10*PI*x_var[0]);
		y_obj[0] = y_obj[0]					+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;
		return;
	} 
	if(!strcmp("ZDT3_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		y_obj[0] = x_var[0];
		y_obj[1] = 1-sqrt(x_var[0])-x_var[0]*sin(10*PI*x_var[0]);
		y_obj[0] = y_obj[0]			+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] =  y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	} 
	if(!strcmp("ZDT3_UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}	
		y_obj[0] = x_var[0];
		y_obj[1] = 1-sqrt(x_var[0])-x_var[0]*sin(10*PI*x_var[0]);
		y_obj[0] = y_obj[0]				+ 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1]	+ 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("ZDT3_UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		y_obj[0] = x_var[0];
		y_obj[1] = 1-sqrt(x_var[0])-x_var[0]*sin(10*PI*x_var[0]);

		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
		y_obj[0] =  y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;		
		return;
	} 
	if(!strcmp("ZDT3_UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 2.0; E = 0.1;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		y_obj[0] = x_var[0];
		y_obj[1] = 1-sqrt(x_var[0])-x_var[0]*sin(10*PI*x_var[0]);
		hj = 2*(0.5/N + E)*sin(2.0*N*PI*x[0]);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("ZDT3_UF7", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x_var[0],0.2);	
		y_obj[0] = x_var[0];
		y_obj[1] = 1-sqrt(x_var[0])-x_var[0]*sin(10*PI*x_var[0]);


		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 
	if(!strcmp("DTLZ6_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*(x[0]-lowBound[0])/(uppBound[0]-lowBound[0]) + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
double g = 0,h=0;
		int n;		
		g = 1 ;
		for (n =0; n < nobj-1; n++) h = h + x_var[n]/(1+g)*(1+sin(3*PI*x_var[n]));
		h = nobj - h;
		if (nobj == 2)
		{
			y_obj[0] = x_var[0];
			y_obj[1] = (1+g) * h;
			
			y_obj[0] = y_obj[0]			+ 2.0 * sum1 / (double)count1;
			y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;
		}
		else if (nobj == 3)
		{}
		else
		{}
		
		return;
	} 	
	if(!strcmp("DTLZ6_UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
double g = 0,h=0;
		int n;		
		g = 1 ;
		for (n =0; n < nobj-1; n++) h = h + x_var[n]/(1+g)*(1+sin(3*PI*x_var[n]));
		h = nobj - h;
		if (nobj == 2)
		{
			y_obj[0] = x_var[0];
			y_obj[1] = (1+g) * h;
			
			y_obj[0] = y_obj[0]			+ 2.0 * sum1 / (double)count1;
			y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;
		}
		else if (nobj == 3)
		{}
		else
		{}

		return;
	} 
	if(!strcmp("DTLZ6_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
double g = 0,h=0;
		int n;		
		g = 1 ;
		for (n =0; n < nobj-1; n++) h = h + x_var[n]/(1+g)*(1+sin(3*PI*x_var[n]));
		h = nobj - h;
		if (nobj == 2)
		{
			y_obj[0] = x_var[0];
			y_obj[1] = (1+g) * h;
			
			y_obj[0] = y_obj[0]			+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
			y_obj[1] = y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		}
		else if (nobj == 3)
		{}
		else
		{}
		return;
	} 
	if(!strcmp("DTLZ6_UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}	
double g = 0,h=0;
		int n;		
		g = 1 ;
		for (n =0; n < nobj-1; n++) h = h + x_var[n]/(1+g)*(1+sin(3*PI*x_var[n]));
		h = nobj - h;
		if (nobj == 2)
		{
			y_obj[0] = x_var[0];
			y_obj[1] = (1+g) * h;
			
			y_obj[0] = y_obj[0]			+ 2.0*sum1 / (double)count1;
			y_obj[1] = y_obj[1] + 2.0*sum2 / (double)count2;
		}
		else if (nobj == 3)
		{}
		else
		{}
		return;
	} 
	if(!strcmp("DTLZ6_UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
double g = 0,h=0;
		int n;		
		g = 1 ;
		for (n =0; n < nobj-1; n++) h = h + x_var[n]/(1+g)*(1+sin(3*PI*x_var[n]));
		h = nobj - h;
		if (nobj == 2)
		{
			y_obj[0] = x_var[0];
			y_obj[1] = (1+g) * h;
			
			hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
			y_obj[0] =  y_obj[0]      + hj + 2.0*sum1 / (double)count1;
			y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;		
		}
		else if (nobj == 3)
		{}
		else
		{}
		
		return;
	} 
	if(!strcmp("DTLZ6_UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 2.0; E = 0.1;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
double g = 0,h=0;
		int n;		
		g = 1 ;
		for (n =0; n < nobj-1; n++) h = h + x_var[n]/(1+g)*(1+sin(3*PI*x_var[n]));
		h = nobj - h;
		if (nobj == 2)
		{
			y_obj[0] = x_var[0];
			y_obj[1] = (1+g) * h;
			
			hj = 2*(0.5/N + E)*sin(2.0*N*PI*x[0]);
			if(hj<0.0) hj = 0.0;
			y_obj[0] = y_obj[0]      + hj + 2.0*sum1 / (double)count1;
			y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;		
		}
		else if (nobj == 3)
		{}
		else
		{}
		
		return;
	} 
	if(!strcmp("DTLZ6_UF7", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 2; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x_var[0],0.2);	
double g = 0,h=0;
		int n;		
		g = 1 ;
		for (n =0; n < nobj-1; n++) h = h + x_var[n]/(1+g)*(1+sin(3*PI*x_var[n]));
		h = nobj - h;
		if (nobj == 2)
		{
			y_obj[0] = x_var[0];
			y_obj[1] = (1+g) * h;
			
			y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
			y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		}
		else if (nobj == 3)
		{}
		else
		{}
			
		
		return;
	} 
	if(!strcmp("WFG1_VAR_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		double A=3,	alfha=3;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 );
		y_obj[1] = 5 * pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ), alfha );
		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 	
	if(!strcmp("WFG1_VAR_UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		double A=3,	alfha=3;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 );
		y_obj[1] = 5 * pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ) , alfha );
		y_obj[0] = y_obj[0]					+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;
		return;
	} 
	if(!strcmp("WFG1_VAR_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		double	A=3,	alfha=3;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 );
		y_obj[1] = 5 * pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ) , alfha );
		y_obj[0] = y_obj[0]			+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] =  y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	} 
	if(!strcmp("WFG1_VAR_UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}	
		double A=3,	alfha=3;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 );
		y_obj[1] = 5 * pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ) , alfha );
		y_obj[0] = y_obj[0]				+ 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1]	+ 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("WFG1_VAR_UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		double	A=3,	alfha=3;
		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 ) + hj;
		y_obj[1] = 5 * (pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ) , alfha ) + hj);

	
		y_obj[0] =  y_obj[0]    + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1]  + 2.0*sum2 / (double)count2;		
		return;
	} 
	if(!strcmp("WFG1_VAR_UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 2.0; E = 0.1;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		double 	A=3,	alfha=3;
		hj = 2*(0.5/N + E)*sin(2.0*N*PI*x[0]);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 )+ hj;
		y_obj[1] = 5 * (pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ) , alfha )+ hj);


		y_obj[0] = y_obj[0]       + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1]  + 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("WFG1_VAR_UF7", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x_var[0],0.2);	
		double 	A=3,	alfha=3;
		y_obj[0] = 1 - cos( yj * PI / 2 );
		y_obj[1] = 5 * pow( 1 - yj - cos( 2 * A * PI * yj + PI / 2 ) / ( 2 * A * PI ) , alfha );

		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;	
		
		return;
	} 
	if(!strcmp("WFG2_VAR_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		double A=5,	alfha=1,	belta = 1.5;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 );
		y_obj[1] = 5 * ( 1 - pow( x_var[0], alfha ) * cos( A * pow( x_var[0], belta ) * PI ) * cos( A * pow( x_var[0], belta) * PI) );

		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 	
	if(!strcmp("WFG2_VAR_UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		double 	A=5,	alfha=1,	belta = 1.5;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 );
		y_obj[1] = 5 * ( 1 - pow( x_var[0], alfha ) * cos( A * pow( x_var[0], belta ) * PI ) * cos( A * pow( x_var[0], belta) * PI) );

		y_obj[0] = y_obj[0]					+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;
		return;
	} 
	if(!strcmp("WFG2_VAR_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		double	A=5,	alfha=1,	belta = 1.5;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 );
		y_obj[1] = 5 * ( 1 - pow( x_var[0], alfha ) * cos( A * pow( x_var[0], belta ) * PI ) * cos( A * pow( x_var[0], belta) * PI) );

		y_obj[0] = y_obj[0]			+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] =  y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	} 
	if(!strcmp("WFG2_VAR_UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj+1; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}	
		double 	A=5,	alfha=1,	belta = 1.5;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 );
		y_obj[1] = 5 * ( 1 - pow( x_var[0], alfha ) * cos( A * pow( x_var[0], belta ) * PI ) * cos( A * pow( x_var[0], belta) * PI) );

		y_obj[0] = y_obj[0]				+ 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1]	+ 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("WFG2_VAR_UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for(j = nobj+1; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		double 	A=5,	alfha=1,	belta = 1.5;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 );
		y_obj[1] = 5 * ( 1 - pow( x_var[0], alfha ) * cos( A * pow( x_var[0], belta ) * PI ) * cos( A * pow( x_var[0], belta) * PI) );

		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
		y_obj[0] =  y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;		
		return;
	} 
	if(!strcmp("WFG2_VAR_UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 2.0; E = 0.1;
		for(j = nobj+1; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		double 	A=5,	alfha=1,	belta = 1.5;
		y_obj[0] = 1 - cos( x_var[0] * PI / 2 );
		y_obj[1] = 5 * ( 1 - pow( x_var[0], alfha ) * cos( A * pow( x_var[0], belta ) * PI ) * cos( A * pow( x_var[0], belta) * PI) );

		hj = 2*(0.5/N + E)*sin(2.0*N*PI*x[0]);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = y_obj[1] + hj + 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("WFG2_VAR_UF7", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj+1; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x_var[0],0.2);	
		double 	A=5,	alfha=1,	belta = 1.5;
		y_obj[0] = 1 - cos( yj * PI / 2 );
		double temp = cos( A * pow( yj, belta ) * PI );
		y_obj[1] = 5 * ( 1 - pow( yj, alfha ) * temp * temp );

		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = y_obj[1] + 2.0 * sum2 / (double)count2;	
		
		return;
	} 
	if(!strcmp("SQRT_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		double A=5,	alfha=5;
		y_obj[0] = x_var[0];
		y_obj[1] = 1-pow(x_var[0], 0.1);
		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 	
	if(!strcmp("SQRT_UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		double A=5,	alfha=5;
		y_obj[0] = x_var[0];
		y_obj[1] = 1-pow(x_var[0], 0.1);
		y_obj[0] = y_obj[0]					+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + 2.0 * sum2 / (double)count2;
		return;
	} 
	if(!strcmp("SQRT_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		double	A=5,	alfha=5;
		y_obj[0] = x_var[0];
		y_obj[1] = 1-pow(x_var[0], 0.1);
		y_obj[0] = y_obj[0]			+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] = 5* y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	} 
	if(!strcmp("SQRT_UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}	
		double A=5,	alfha=5;
		y_obj[0] = x_var[0];
		y_obj[1] = 1-pow(x_var[0], 0.1);
		y_obj[0] = y_obj[0]				+ 2.0*sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1]	+ 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("SQRT_UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		double	A=5,	alfha=5;
		y_obj[0] = x_var[0];
		y_obj[1] = 1-pow(x_var[0], 0.1);
		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
		y_obj[0] =  y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + hj + 2.0*sum2 / (double)count2;		
		return;
	} 
	if(!strcmp("SQRT_UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 2.0; E = 0.1;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		double 	A=5,	alfha=5;
		y_obj[0] = x_var[0];
		y_obj[1] = 1-pow(x_var[0], 0.1);

		hj = 2*(0.5/N + E)*sin(2.0*N*PI*x[0]);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + hj + 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("SQRT_UF7", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x_var[0],0.2);	
		double 	A=5,	alfha=5;
		y_obj[0] = yj;
		y_obj[1] = 1-pow(yj, 0.1);

		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + 2.0 * sum2 / (double)count2;	
		
		return;
	} 
	if(!strcmp("FONSICA_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*(x[0]+1)/2 + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		double A=5,	alfha=5;
		y_obj[0] = ( exp( -2.0*(x_var[0]-1)*(x_var[0]-1) ) - exp(-8.0) ) / (1.0-exp(-8.0));
		y_obj[1] = ( exp( -2.0*(x_var[0]+1)*(x_var[0]+1) ) - exp(-8.0) ) / (1-exp(-8.0));
		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 	
	if(!strcmp("FONSICA_UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*(x[0]+1)/2*((x[0]+1)/2*cos(4.0*(6.0*PI*(x[0]+1)/2+j*PI/nx))+2.0)*cos(6.0*PI*(x[0]+1)/2+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*(x[0]+1)/2*((x[0]+1)/2*cos(4.0*(6.0*PI*(x[0]+1)/2+j*PI/nx))+2.0)*sin(6.0*PI*(x[0]+1)/2+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		double A=5,	alfha=5;
		y_obj[0] = ( exp( -2.0*(x_var[0]-1)*(x_var[0]-1) ) - exp(-8.0) ) / (1.0-exp(-8.0));
		y_obj[1] = ( exp( -2.0*(x_var[0]+1)*(x_var[0]+1) ) - exp(-8.0) ) / (1-exp(-8.0));
		y_obj[0] = y_obj[0]					+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + 2.0 * sum2 / (double)count2;
		return;
	} 
	if(!strcmp("FONSICA_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-pow((x[0]+1)/2,0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		double	A=5,	alfha=5;
			y_obj[0] = ( exp( -2.0*(x_var[0]-1)*(x_var[0]-1) ) - exp(-8.0) ) / (1.0-exp(-8.0));
		y_obj[1] = ( exp( -2.0*(x_var[0]+1)*(x_var[0]+1) ) - exp(-8.0) ) / (1-exp(-8.0));
		y_obj[0] = y_obj[0]			+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] =  5*y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	} 
	if(!strcmp("FONSICA_UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*(x[0]+1)/2+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}	
		double A=5,	alfha=5;
		y_obj[0] = ( exp( -2.0*(x_var[0]-1)*(x_var[0]-1) ) - exp(-8.0) ) / (1.0-exp(-8.0));
		y_obj[1] = ( exp( -2.0*(x_var[0]+1)*(x_var[0]+1) ) - exp(-8.0) ) / (1-exp(-8.0));
		y_obj[0] = y_obj[0]				+ 2.0*sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1]	+ 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("FONSICA_UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*(x[0]+1)/2);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		double	A=5,	alfha=5;
			y_obj[0] = ( exp( -2.0*(x_var[0]-1)*(x_var[0]-1) ) - exp(-8.0) ) / (1.0-exp(-8.0));
		y_obj[1] = ( exp( -2.0*(x_var[0]+1)*(x_var[0]+1) ) - exp(-8.0) ) / (1-exp(-8.0));

		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
		y_obj[0] =  y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + hj + 2.0*sum2 / (double)count2;		
		return;
	} 
	if(!strcmp("FONSICA_UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 2.0; E = 0.1;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*(x[0]+1)/2+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}

		y_obj[0] = ( exp( -2.0*(x_var[0]-1)*(x_var[0]-1) ) - exp(-8.0) ) / (1.0-exp(-8.0));
		y_obj[1] = ( exp( -2.0*(x_var[0]+1)*(x_var[0]+1) ) - exp(-8.0) ) / (1-exp(-8.0));

		hj = 2*(0.5/N + E)*sin(2.0*N*PI*x[0]);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + hj + 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("FONSICA_UF7", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{			
			yj = x[j-1] - sin(6.0*PI*((x[0]+1)/2+1)/2+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x_var[0],0.2);	

		y_obj[0] = ( exp( -2.0*(yj-1)*(yj-1) ) - exp(-8.0) ) / (1-exp(-8.0));
		y_obj[1] = ( exp( -2.0*(yj+1)*(yj+1) ) - exp(-8.0) ) / (1-exp(-8.0));

		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + 2.0 * sum2 / (double)count2;	
		
		return;
	} 
	if(!strcmp("LZDTZ3_UF1", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		double A=5,	alfha=5;
		y_obj[0] = x_var[0];
		if (x_var[0] <= 0.05)
			y_obj[1] = 1.0 - 19.0 * x_var[0];
		else
			y_obj[1] = 1.0/19.0 - x_var[0]/19.0;
		
		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + 2.0 * sum2 / (double)count2;		
		
		return;
	} 	
	if(!strcmp("LZDTZ3_UF2", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			if(j % 2 == 1) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			} 
			else 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
		}
		double A=5,	alfha=5;
		y_obj[0] = x_var[0];
		if (x_var[0] <= 0.05)
			y_obj[1] = 1.0 - 19.0 * x_var[0];
		else
			y_obj[1] = 1.0/19.0 - x_var[0]/19.0;

		y_obj[0] = y_obj[0]					+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + 2.0 * sum2 / (double)count2;
		return;
	} 
	if(!strcmp("LZDTZ3_UF3", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, prod1, prod2, yj, pj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		prod1  = prod2  = 1.0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
		}
		double	A=5,	alfha=5;
		y_obj[0] = x_var[0];
		if (x_var[0] <= 0.05)
			y_obj[1] = 1.0 - 19.0 * x_var[0];
		else
			y_obj[1] = 1.0/19.0 - x_var[0]/19.0;

		y_obj[0] = y_obj[0]			+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] =  5*y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		return;
	} 
	if(!strcmp("LZDTZ3_UF4", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}	
		double A=5,	alfha=5;
		y_obj[0] = x_var[0];
		if (x_var[0] <= 0.05)
			y_obj[1] = 1.0 - 19.0 * x_var[0];
		else
			y_obj[1] = 1.0/19.0 - x_var[0]/19.0;

		y_obj[0] = y_obj[0]				+ 2.0*sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1]	+ 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("LZDTZ3_UF5", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 10.0; E = 0.1;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 2 == 0) 
			{
				sum2  += hj;
				count2++;
			} 
			else 
			{
				sum1  += hj;
				count1++;
			}
		}
		double	A=5,	alfha=5;
		y_obj[0] = x_var[0];
		if (x_var[0] <= 0.05)
			y_obj[1] = 1.0 - 19.0 * x_var[0];
		else
			y_obj[1] = 1.0/19.0 - x_var[0]/19.0;

		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]));
		y_obj[0] =  y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + hj + 2.0*sum2 / (double)count2;		
		return;
	} 
	if(!strcmp("LZDTZ3_UF6", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj, hj, N, E;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		N = 2.0; E = 0.1;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		double 	A=5,	alfha=5;
		y_obj[0] = x_var[0];
		if (x_var[0] <= 0.05)
			y_obj[1] = 1.0 - 19.0 * x_var[0];
		else
			y_obj[1] = 1.0/19.0 - x_var[0]/19.0;

		hj = 2*(0.5/N + E)*sin(2.0*N*PI*x[0]);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + hj + 2.0*sum2 / (double)count2;
		return;
	} 
	if(!strcmp("LZDTZ3_UF7", strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 2 == 0) 
			{
				sum2  += yj*yj;
				count2++;
			} 
			else 
			{
				sum1  += yj*yj;
				count1++;
			}
		}
		yj = pow(x_var[0],0.2);	
		double 	A=5,	alfha=5;
		y_obj[0] = x_var[0];
		if (yj <= 0.05)
			y_obj[1] = 1.0 - 19.0 * yj;
		else
			y_obj[1] = 1.0/19.0 - yj/19.0;
	
		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 5*y_obj[1] + 2.0 * sum2 / (double)count2;	
		
		return;
	} 
	if(!strcmp("WFG1_UF1", strTestInstance))
	{
		unsigned int j, count3, count1, count2;
		double sum3, sum1, sum2, yj;
		
		sum1   = sum2 = sum3  = 0.0;
		count1 = count2 = count3 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 3 == 0) 
			{
				sum3 += yj;
				count3++;
			} 
			else if (j % 3 == 1) 
			{
				sum1 += yj;
				count1++;
			}
			else
			{
				sum2 += yj;
				count2++;
			}
		}

		double A=5,	alfha=5;

		y_obj[0] = (1 - cos( x_var[0] * PI / 2 ))*(1 - cos( x_var[1] * PI / 2 ));
		y_obj[1] = (1 - cos( x_var[0] * PI / 2 ))*(1 - sin( x_var[1] * PI / 2 ));
		y_obj[2] = 5 * pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ), alfha );

		y_obj[0] = y_obj[0]	+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 2*y_obj[1]	+ 2.0 * sum2 / (double)count2;
		y_obj[2] = 4*y_obj[2] + 2.0 * sum3 / (double)count3;		
		
		return;
	} 	
	if(!strcmp("WFG1_UF2", strTestInstance))
	{
		unsigned int j, count3, count1, count2;
		double sum3, sum1, sum2, yj;
		
		sum1   = sum2 = sum3  = 0.0;
		count1 = count2 = count3 = 0;
	
		for(j = nobj; j <= nx; j++) 
		{
			if(j % 3 == 0) 
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*cos(6.0*PI*x[0]+j*PI/nx);
				sum3 += yj*yj;
				count3++;
			} 
			else if (j % 3 == 1)
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum1 += yj*yj;
				count1++;
			}
			else
			{
				yj = x[j-1]-0.3*x[0]*(x[0]*cos(4.0*(6.0*PI*x[0]+j*PI/nx))+2.0)*sin(6.0*PI*x[0]+j*PI/nx);
				sum2 += yj*yj;
				count2++;
			}
		}
		double A=5,	alfha=5;
		y_obj[0] = (1 - cos( x_var[0] * PI / 2 ))*(1 - cos( x_var[1] * PI / 2 ));
		y_obj[1] = (1 - cos( x_var[0] * PI / 2 ))*(1 - sin( x_var[1] * PI / 2 ));
		y_obj[2] = 5 * pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ), alfha );

		y_obj[0] = y_obj[0]	+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 2*y_obj[1]	+ 2.0 * sum2 / (double)count2;
		y_obj[2] = 4*y_obj[2] + 2.0 * sum3 / (double)count3;

		return;
	} 
	if(!strcmp("WFG1_UF3", strTestInstance))
	{
		unsigned int j, count3, count1, count2;
		double sum3, sum1, sum2, prod3, prod1, prod2, yj, pj;
		
		sum3 = sum1   = sum2   = 0.0;
		count3 = count1 = count2 = 0;
		prod3 = prod1  = prod2  = 1.0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-pow(x[0],0.5*(1.0+3.0*(j-2.0)/(nx-2.0)));
			pj = cos(20.0*yj*PI/sqrt(j+0.0));
			if (j % 3 == 0) 
			{
				sum3  += yj*yj;
				prod3 *= pj;
				count3++;
			} 
			else if (j % 3 == 1)
			{
				sum1  += yj*yj;
				prod1 *= pj;
				count1++;
			}
			else 
			{
				sum2  += yj*yj;
				prod2 *= pj;
				count2++;
			}
		}
		double	A=5,	alfha=5;
		y_obj[0] = (1 - cos( x_var[0] * PI / 2 ))*(1 - cos( x_var[1] * PI / 2 ));
		y_obj[1] = (1 - cos( x_var[0] * PI / 2 ))*(1 - sin( x_var[1] * PI / 2 ));
		y_obj[2] = 5 * pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ), alfha );

		y_obj[0] = y_obj[0]	+ 2.0*(4.0*sum1 - 2.0*prod1 + 2.0) / (double)count1;
		y_obj[1] = 2*y_obj[1] + 2.0*(4.0*sum2 - 2.0*prod2 + 2.0) / (double)count2;
		y_obj[2] = 4*y_obj[2] + 2.0*(4.0*sum3 - 2.0*prod3 + 2.0) / (double)count3;

		return;
	} 
	if(!strcmp("WFG1_UF4", strTestInstance))
	{
		unsigned int j, count3, count1, count2;
		double sum3, sum1, sum2, yj, hj;
		
		sum1   = sum2 = sum3  = 0.0;
		count1 = count2 = count3 = 0;

		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1]-sin(6.0*PI*x[0]+j*PI/nx);
			hj = fabs(yj)/(1.0+exp(2.0*fabs(yj)));
			if (j % 3 == 0) 
			{
				sum3  += hj;
				count3++;
			} 
			else if (j % 3 == 1)
			{
				sum1  += hj;
				count1++;
			}
			else
			{
				sum2  += hj;
				count2++;
			}
		}	
		double A=5,	alfha=5;
		y_obj[0] = (1 - cos( x_var[0] * PI / 2 ))*(1 - cos( x_var[1] * PI / 2 ));
		y_obj[1] = (1 - cos( x_var[0] * PI / 2 ))*(1 - sin( x_var[1] * PI / 2 ));
		y_obj[2] = 5 * pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ), alfha );

		y_obj[0] = y_obj[0]	+ 2.0*sum1 / (double)count1;
		y_obj[1] = 2*y_obj[1]	+ 2.0*sum2 / (double)count2;
		y_obj[2] = 4*y_obj[2]	+ 2.0*sum3 / (double)count3;
		return;
	} 
	if(!strcmp("WFG1_UF5", strTestInstance))
	{
		unsigned int j, count1, count2,count3;
		double sum1, sum2,sum3, yj, hj, N, E;
		
		sum1   = sum2  =sum3 = 0.0;
		count1 = count2 = count3 = 0;
		N = 10.0; E = 0.1;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]);
			hj = 2.0*yj*yj - cos(4.0*PI*yj) + 1.0;
			if (j % 3 == 0) 
			{
				sum3  += hj;
				count3++;
			} 
			else if (j % 3 == 1) 
			{
				sum1  += hj;
				count1++;
			}
			else
			{
				sum2  += hj;
				count2++;
			}
		}
		double	A=5,	alfha=5;
		y_obj[0] = (1 - cos( x_var[0] * PI / 2 ))*(1 - cos( x_var[1] * PI / 2 ));
		y_obj[1] = (1 - cos( x_var[0] * PI / 2 ))*(1 - sin( x_var[1] * PI / 2 ));
		y_obj[2] = 5 * pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ), alfha );

		hj = (0.5/N + E)*fabs(sin(2.0*N*PI*x[0]))*fabs(sin(2.0*N*PI*x[1]));
		y_obj[0] =  y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = 2*y_obj[1] + hj + 2.0*sum2 / (double)count2;
		y_obj[2] = 4*y_obj[2] + hj + 2.0*sum3 / (double)count3;
		return;
	} 
	if(!strcmp("WFG1_UF6", strTestInstance))
	{
		unsigned int j, count1, count2,count3;
		double sum1, sum2,sum3, yj, hj, N, E;
		
		sum1   = sum2 =sum3  = 0.0;
		count1 = count2=count3 = 0;
		N = 2.0; E = 0.1;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 3 == 0) 
			{
				sum3  += yj*yj;
				count3++;
			} 
			else if (j % 3 == 1)
			{
				sum1  += yj*yj;
				count1++;
			}
			else
			{
				sum2  += yj*yj;
				count2++;
			}
		}
		double 	A=5,	alfha=5;
		y_obj[0] = (1 - cos( x_var[0] * PI / 2 ))*(1 - cos( x_var[1] * PI / 2 ));
		y_obj[1] = (1 - cos( x_var[0] * PI / 2 ))*(1 - sin( x_var[1] * PI / 2 ));
		y_obj[2] = 5 * pow( 1 - x_var[0] - cos( 2 * A * PI * x_var[0] + PI / 2 ) / ( 2 * A * PI ), alfha );

		hj = 2*(0.5/N + E)*sin(2.0*N*PI*x_var[0])*sin(2.0*N*PI*x_var[1]);
		if(hj<0.0) hj = 0.0;
		y_obj[0] = y_obj[0]      + hj + 2.0*sum1 / (double)count1;
		y_obj[1] = 2*y_obj[1] + hj + 2.0*sum2 / (double)count2;
		y_obj[2] = 4*y_obj[2] + hj + 2.0*sum3 / (double)count3;
		return;
	} 
	if(!strcmp("WFG1_UF7", strTestInstance))
	{
		unsigned int j, count1, count2,count3;
		double sum1, sum2,sum3, yj;
		
		sum1   = sum2  =sum3 = 0.0;
		count1 = count2 = count3 = 0;
		for(j = nobj; j <= nx; j++) 
		{
			yj = x[j-1] - sin(6.0*PI*x[0]+j*PI/nx);
			if (j % 3 == 0) 
			{
				sum3  += yj*yj;
				count3++;
			} 
			else if (j % 3 == 1) 
			{
				sum1  += yj*yj;
				count1++;
			}
			else
			{
				sum2  += yj*yj;
				count2++;
			}
		}
		double y1 = pow(x_var[0],0.2);	
		double y2 = pow(x_var[1],0.2);
		double 	A=5,	alfha=5;
		y_obj[0] = (1 - cos( y1 * PI / 2 ))*(1 - cos( y2 * PI / 2 ));
		y_obj[1] = (1 - cos( y1 * PI / 2 ))*(1 - sin( y2 * PI / 2 ));
		y_obj[2] = 5 * pow( 1 - y1 - cos( 2 * A * PI * y1 + PI / 2 ) / ( 2 * A * PI ), alfha );

		y_obj[0] = y_obj[0]				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = 2*y_obj[1] + 2.0 * sum2 / (double)count2;	
		y_obj[2] = 4*y_obj[2] + 2.0 * sum3 / (double)count3;	
		
		return;
	}
		if(!strcmp("AB", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	
	if(!strcmp("ABC1", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(10);
		temp_object.push_back(50);
		y_obj[1] = dist_vector(temp_object, x_var);	
		temp_object.clear();

		temp_object.clear();
		temp_object.push_back(90);
		temp_object.push_back(50);
		y_obj[2] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 

	if(!strcmp("ABC2", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(48);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		temp_object.clear();

		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(10);
		y_obj[2] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 

	if(!strcmp("ABCDNONCON1", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(45);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(15);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(55);
		temp_object.push_back(15);
		y_obj[3] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 

	if(!strcmp("ABCDNONCON2", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(70);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(35);
		temp_object.push_back(30);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(40);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(65);
		temp_object.push_back(30);
		y_obj[3] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABCD13", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(48);
		temp_object.push_back(11);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(11);
		y_obj[3] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 

	if(!strcmp("ABCD22", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(10);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(90);
		y_obj[3] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABCDEF1", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(48);
		temp_object.push_back(11);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(11);
		y_obj[3] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(48);
		temp_object.push_back(89);
		y_obj[4] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(89);
		y_obj[5] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABCDEF2", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(45);
		temp_object.push_back(15);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(55);
		temp_object.push_back(15);
		y_obj[3] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(45);
		temp_object.push_back(85);
		y_obj[4] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(55);
		temp_object.push_back(85);
		y_obj[5] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABCDEF3", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(30);
		temp_object.push_back(20);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(70);
		temp_object.push_back(20);
		y_obj[3] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(30);
		temp_object.push_back(80);
		y_obj[4] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(70);
		temp_object.push_back(80);
		y_obj[5] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABCD1", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(70);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(68);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(28);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(30);
		y_obj[3] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	/*
	if(!strcmp("AB", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABC", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(10);
		temp_object.push_back(50);
		y_obj[1] = dist_vector(temp_object, x_var);	
		temp_object.clear();

		temp_object.clear();
		temp_object.push_back(90);
		temp_object.push_back(50);
		y_obj[2] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 

	if(!strcmp("ABCDNONCON1", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(45);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(15);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(55);
		temp_object.push_back(15);
		y_obj[3] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 

	if(!strcmp("ABCDNONCON2", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(70);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(35);
		temp_object.push_back(30);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(40);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(65);
		temp_object.push_back(30);
		y_obj[3] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 

	if(!strcmp("ABCD13", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(48);
		temp_object.push_back(11);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(11);
		y_obj[3] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABCD1", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(70);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(68);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(28);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(30);
		y_obj[3] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	}
	if(!strcmp("ABCDEF1", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(48);
		temp_object.push_back(88);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(48);
		temp_object.push_back(12);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[3] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(12);
		y_obj[4] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(52);
		temp_object.push_back(88);
		y_obj[5] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	}
	if(!strcmp("ABCDEF2", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(70);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(45);
		temp_object.push_back(65);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(45);
		temp_object.push_back(35);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(30);
		y_obj[3] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(55);
		temp_object.push_back(35);
		y_obj[4] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(55);
		temp_object.push_back(65);
		y_obj[5] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABCDEF3", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(70);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(20);
		temp_object.push_back(65);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(20);
		temp_object.push_back(35);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(30);
		y_obj[3] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(80);
		temp_object.push_back(35);
		y_obj[4] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(80);
		temp_object.push_back(65);
		y_obj[5] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} */
	/*
	if(!strcmp("ABCD22", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(51);
		temp_object.push_back(90);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(51);
		temp_object.push_back(10);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[3] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABCDEF1", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(49);
		temp_object.push_back(89);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(49);
		temp_object.push_back(11);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[3] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(51);
		temp_object.push_back(11);
		y_obj[4] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(51);
		temp_object.push_back(89);
		y_obj[5] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABCDEF2", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(45);
		temp_object.push_back(15);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(55);
		temp_object.push_back(15);
		y_obj[3] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(45);
		temp_object.push_back(85);
		y_obj[4] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(55);
		temp_object.push_back(85);
		y_obj[5] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	if(!strcmp("ABCDEF3", strTestInstance))
	{	
		vector<double> temp_object;
		//coordinate of object A
		temp_object.push_back(50);
		temp_object.push_back(90);
		y_obj[0] = dist_vector(temp_object, x_var);

		temp_object.clear();
		temp_object.push_back(20);
		temp_object.push_back(85);
		y_obj[1] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(20);
		temp_object.push_back(15);
		y_obj[2] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(50);
		temp_object.push_back(10);
		y_obj[3] = dist_vector(temp_object, x_var);	
		
		temp_object.clear();
		temp_object.push_back(80);
		temp_object.push_back(15);
		y_obj[4] = dist_vector(temp_object, x_var);	

		temp_object.clear();
		temp_object.push_back(80);
		temp_object.push_back(85);
		y_obj[5] = dist_vector(temp_object, x_var);	
		temp_object.clear();
		
		return;
	} 
	*/
	if(!strcmp("ZDT1_4",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);

		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(y_obj[0]/g));
		y_obj[2] = 1*y_obj[0] +0.1*y_obj[1];
		y_obj[3] = 0.1*y_obj[0] +1*y_obj[1];
		return;
	}

	if(!strcmp("ZDT2_4",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - pow(y_obj[0]/g,2));
		y_obj[2] = 1*y_obj[0] +0.1*y_obj[1];
		y_obj[3] = 0.1*y_obj[0] +1*y_obj[1];
		return;
	}

	if(!strcmp("ZDT3_4",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);

		y_obj[0] = x_var[0];
		y_obj[1] = g*(1 - sqrt(x_var[0]/g) - x_var[0]*sin(10*PI*x_var[0])/g);
		y_obj[2] = 1*y_obj[0] +0.1*y_obj[1];
		y_obj[3] = 0.1*y_obj[0] +1*y_obj[1];
		return;
	}

	if(!strcmp("ZDT4_4",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
		{
			double x = x_var[n];
			g+= x*x - 10*cos(4*PI*x);
		}
		g = 1 + 10*(nvar-1) + g;
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1- sqrt(y_obj[0]/g));
		y_obj[2] = 1*y_obj[0] +0.1*y_obj[1];
		y_obj[3] = 0.1*y_obj[0] +1*y_obj[1];
		return;
	}

	if(!strcmp("ZDT6_4",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n]/(nvar - 1);
		g = 1 + 9*pow(g,0.25) ;

		y_obj[0] = 1 - exp(-4*x_var[0])*pow(sin(6*PI*x_var[0]),6);
		y_obj[1] = g*(1- pow(y_obj[0]/g,2));
		y_obj[2] = 1*y_obj[0] +0.1*y_obj[1];
		y_obj[3] = 0.1*y_obj[0] +1*y_obj[1];
		return;
	}

	if(!strcmp("DTLZ1_6",strTestInstance))
	{
		double g = 0;
		for(int n=2; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - 3 + g);
		y_obj[0] = 0.5*(1 + g)*x_var[0]*x_var[1];
		y_obj[1] = 0.5*(1 + g)*x_var[0]*(1 - x_var[1]);
		y_obj[2] = 0.5*(1 + g)*(1 - x_var[0]);
		y_obj[3] = 10 * y_obj[0] + 1*y_obj[1]+0.1*y_obj[2];
		y_obj[4] = 0.1 * y_obj[0] +10*y_obj[1]+1*y_obj[2];
		y_obj[5] = 1 * y_obj[0] +0.1*y_obj[1]+10*y_obj[2];
		return;
	}

	if(!strcmp("DTLZ2_6",strTestInstance))
	{
		double g = 0;		
		for(int n=2; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}
		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[2] = (1 + g)*sin(x_var[0]*PI/2);
		y_obj[3] = 10 * y_obj[0] + 1*y_obj[1]+0.1*y_obj[2];
		y_obj[4] = 0.1 * y_obj[0] +10*y_obj[1]+1*y_obj[2];
		y_obj[5] = 1 * y_obj[0] +0.1*y_obj[1]+10*y_obj[2];
		return;
	}

	if(!strcmp("DTLZ3_6",strTestInstance))
	{
		double g = 0;
		for(int n=2; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - 3 + g);

		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[2] = (1 + g)*sin(x_var[0]*PI/2);
		y_obj[3] = 10 * y_obj[0] + 1*y_obj[1]+0.1*y_obj[2];
		y_obj[4] = 0.1 * y_obj[0] +10*y_obj[1]+1*y_obj[2];
		y_obj[5] = 1 * y_obj[0] +0.1*y_obj[1]+10*y_obj[2];
		return;
	}

	if(!strcmp("DTLZ4_6",strTestInstance))
	{
		double g = 0;
		double alpha = 100;
		for(int n=2; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}

		y_obj[0] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*cos(pow(x_var[1],alpha)*PI/2);
		y_obj[1] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*sin(pow(x_var[1],alpha)*PI/2);
		y_obj[2] = (1 + g)*sin(pow(x_var[0],alpha)*PI/2);
		y_obj[3] = 10 * y_obj[0] + 1*y_obj[1]+0.1*y_obj[2];
		y_obj[4] = 0.1 * y_obj[0] +10*y_obj[1]+1*y_obj[2];
		y_obj[5] = 1 * y_obj[0] +0.1*y_obj[1]+10*y_obj[2];
		return;
	}

	if(!strcmp("DTLZ6_6",strTestInstance))
	{
		double g = 0,h=0;
		int n;
		for(n=2; n<nvar;n++)	g = g + x_var[n];
		g = 1 + 9*g /(nvar + 1 - 3);
		for (n =0; n < 3-1; n++) h = h + x_var[n]/(1+g)*(1+sin(3*PI*x_var[n]));
		h = 3 - h;

		y_obj[0] = x_var[0];
		y_obj[1] = x_var[1];
		y_obj[2] = (1+g) * h;
		y_obj[3] = 10 * y_obj[0] + 1*y_obj[1]+0.1*y_obj[2];
		y_obj[4] = 0.1 * y_obj[0] +10*y_obj[1]+1*y_obj[2];
		y_obj[5] = 1 * y_obj[0] +0.1*y_obj[1]+10*y_obj[2];
		return;
	}

	if(!strcmp("TDY1",strTestInstance))
	{
		unsigned int j, count1, count2;
		double sum1, sum2, yj;
		
		sum1   = sum2   = 0.0;
		count1 = count2 = 0;
		for(j = 3; j <= nx; j++) 
		{
			yj = x_var[j-1] - sin(6.0*PI*x_var[0] + j*PI/nx);
			yj = yj * yj;
			if(j % 2 == 0) 
			{
				sum2 += yj;
				count2++;
			} 
			else 
			{
				sum1 += yj;
				count1++;
			}
		}
		double h = (x_var[0] + x_var[1]) / 2;
		double tmp = fabs( cos(h * PI / 0.6) );
		y_obj[0] = -h - tmp				+ 2.0 * sum1 / (double)count1;
		y_obj[1] = -1.0 + h - tmp		+ 2.0 * sum2 / (double)count2;
	}

	if(!strcmp("TDY2",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nx;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nx-1);
		y_obj[0] = 2*x_var[0]+4*cos(8*PI*x_var[0])+g*g;
		y_obj[1] = 8*exp(-2*x_var[0]*x_var[0])+g*g;
	}

	if(!strcmp("TDY3",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nx;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nx-1);

		y_obj[0] = x_var[0] + 2*sin(8*PI*x_var[0])+sin(4*PI*x_var[0])+g*g;
		y_obj[1] = 8*exp(-2*x_var[0]*x_var[0])+4*exp(-2*(x_var[0]-0.2)*(x_var[0]-0.2))+g*g;
	}

	if(!strcmp("TDY4",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nx;n++)
			g+= x[n];
		g = 1 + 9*g/(nx-1);

		y_obj[0] = -2*x[0]+g*g;
		y_obj[1] = 2*x[0]+32*sin(2*PI*x[0])*sin(2*PI*x[0])+g*g;
	}

	if(!strcmp("TDY5",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nx;n++)
			g+= x[n];
		g = 1 + 9*g/(nx-1);
		
		double y = pow(x_var[0],0.7);
		y_obj[0] = x_var[0]+g*g;
		y_obj[1] = exp(-y)+sin(2*PI*y)+g*g;
	}

	if(!strcmp("TDY6",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nx;n++)
			g+= x[n];
		g = 1 + 9*g/(nx-1);
		
		double y = log(x_var[0]);

		y_obj[0] = x_var[0] + g*g;
		y_obj[1] = -0.5*y+sin(4*PI*y)+g*g;
	}

	if(!strcmp("DTLZ1_2_6",strTestInstance))
	{
		double g = 0;
		for(int n=nobj-1; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);
		double cc = 2;
		double temp2=1.0/cc*(1+2*g*x_var[2])/(1+g);
		double temp3=1.0/cc*(1+2*g*x_var[3])/(1+g);
		double temp4=1.0/cc*(1+2*g*x_var[4])/(1+g);
		y_obj[0] = 0.5*(1 + g)*x_var[0]*x_var[1]*temp2*temp3*temp4;
		y_obj[1] = 0.5*(1 + g)*x_var[0]*x_var[1]*temp2*temp3*(1-temp4);
		y_obj[2] = 0.5*(1 + g)*x_var[0]*x_var[1]*temp2*(1-temp3);
		y_obj[3] = 0.5*(1 + g)*x_var[0]*x_var[1]*(1-temp2);
		y_obj[4] = 0.5*(1 + g)*x_var[0]*(1-x_var[1]);
		y_obj[5] = 0.5*(1 + g)*(1-x_var[0]);
		return;
	}

	if(!strcmp("DTLZ2_2_6",strTestInstance))
	{
		double g = 0;		
		for(int n=nobj-1; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}
		double temp2=PI/4*(1+2*g*x_var[2])/(1+g);
		double temp3=PI/4*(1+2*g*x_var[3])/(1+g);
		double temp4=PI/4*(1+2*g*x_var[4])/(1+g);
		int I=3;
		y_obj[0] = pow(2,(nobj-I)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*cos(temp3)*cos(temp4);
		y_obj[1] = pow(2,(nobj-I)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*cos(temp3)*sin(temp4);
		y_obj[2] = pow(2,(nobj-I-1)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*sin(temp3);
		y_obj[3] = pow(2,(nobj-I-2)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(temp2);
		y_obj[4] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[5] = (1 + g)*sin(x_var[0]*PI/2);
		return;
	}

	if(!strcmp("DTLZ3_2_6",strTestInstance))
	{
		double g = 0;
		for(int n=nobj-1; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);

		double temp2=PI/4*(1+2*g*x_var[2])/(1+g);
		double temp3=PI/4*(1+2*g*x_var[3])/(1+g);
		double temp4=PI/4*(1+2*g*x_var[4])/(1+g);
		int I=3;
		y_obj[0] = pow(2,(nobj-I)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*cos(temp3)*cos(temp4);
		y_obj[1] = pow(2,(nobj-I)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*cos(temp3)*sin(temp4);
		y_obj[2] = pow(2,(nobj-I-1)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*sin(temp3);
		y_obj[3] = pow(2,(nobj-I-2)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(temp2);
		y_obj[4] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[5] = (1 + g)*sin(x_var[0]*PI/2);
		return;
	}

	if(!strcmp("DTLZ4_2_6",strTestInstance))
	{
		double g = 0;
		double alpha = 100;
		for(int n=nobj-1; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}

		double temp2=PI/4*(1+2*g*pow(x_var[2],alpha))/(1+g);
		double temp3=PI/4*(1+2*g*pow(x_var[3],alpha))/(1+g);
		double temp4=PI/4*(1+2*g*pow(x_var[4],alpha))/(1+g);
		int I=3;
		y_obj[0] = pow(2,(nobj-I)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*cos(temp3)*cos(temp4);
		y_obj[1] = pow(2,(nobj-I)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*cos(temp3)*sin(temp4);
		y_obj[2] = pow(2,(nobj-I-1)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*sin(temp3);
		y_obj[3] = pow(2,(nobj-I-2)*1.0/2) * (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(temp2);
		y_obj[4] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[5] = (1 + g)*sin(x_var[0]*PI/2);

		return;
	}
	if(!strcmp("DTLZ1_4_4",strTestInstance))
	{
		double g = 0;
		for(int n=nobj-1; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);
		y_obj[0] = 0.5*(1 + g)*x_var[0]*x_var[1]*x_var[2];
		y_obj[1] = 0.5*(1 + g)*x_var[0]*x_var[1]*(1 - x_var[2]);
		y_obj[2] = 0.5*(1 + g)*x_var[0]*(1 - x_var[1]);
		y_obj[3] = 0.5*(1 + g)*(1-x_var[0]);
		return;
	}

	if(!strcmp("DTLZ2_4_4",strTestInstance))
	{
		double g = 0;		
		for(int n=nobj-1; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}
		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(x_var[2]*PI/2);
		y_obj[2] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[3] = (1 + g)*sin(x_var[0]*PI/2);
		return;
	}

	if(!strcmp("DTLZ3_4_4",strTestInstance))
	{
		double g = 0;
		for(int n=nobj-1; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);

		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(x_var[2]*PI/2);
		y_obj[2] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[3] = (1 + g)*sin(x_var[0]*PI/2);
		return;
	}

	if(!strcmp("DTLZ4_4_4",strTestInstance))
	{
		double g = 0;
		double alpha = 100;
		for(int n=nobj-1; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}

		y_obj[0] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*cos(pow(x_var[1],alpha)*PI/2)*cos(pow(x_var[2],alpha)*PI/2);
		y_obj[1] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*cos(pow(x_var[1],alpha)*PI/2)*sin(pow(x_var[2],alpha)*PI/2);
		y_obj[2] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*sin(pow(x_var[1],alpha)*PI/2);
		y_obj[3] = (1 + g)*sin(pow(x_var[0],alpha)*PI/2);
		return;
	}
	if(!strcmp("DTLZ1_5_5",strTestInstance))
	{
		double g = 0;
		for(int n=nobj-1; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);

		y_obj[0] = 0.5*(1 + g)*x_var[0]*x_var[1]*x_var[2]*x_var[3];
		y_obj[1] = 0.5*(1 + g)*x_var[0]*x_var[1]*x_var[2]*(1-x_var[3]);
		y_obj[2] = 0.5*(1 + g)*x_var[0]*x_var[1]*(1 - x_var[2]);
		y_obj[3] = 0.5*(1 + g)*x_var[0]*(1 - x_var[1]);
		y_obj[4] = 0.5*(1 + g)*(1-x_var[0]);
		return;
	}

	if(!strcmp("DTLZ2_5_5",strTestInstance))
	{
		double g = 0;		
		for(int n=nobj-1; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}

		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*cos(x_var[3]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*sin(x_var[3]*PI/2);
		y_obj[2] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(x_var[2]*PI/2);
		y_obj[3] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[4] = (1 + g)*sin(x_var[0]*PI/2);

		return;
	}

	if(!strcmp("DTLZ3_5_5",strTestInstance))
	{
		double g = 0;
		for(int n=nobj-1; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);

		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*cos(x_var[3]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*sin(x_var[3]*PI/2);
		y_obj[2] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(x_var[2]*PI/2);
		y_obj[3] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[4] = (1 + g)*sin(x_var[0]*PI/2);
		return;
	}

	if(!strcmp("DTLZ4_5_5",strTestInstance))
	{
		double g = 0;
		double alpha = 100;
		for(int n=nobj-1; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}

		y_obj[0] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*cos(pow(x_var[1],alpha)*PI/2)*cos(pow(x_var[2],alpha)*PI/2)*cos(pow(x_var[3],alpha)*PI/2);
		y_obj[1] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*cos(pow(x_var[1],alpha)*PI/2)*cos(pow(x_var[2],alpha)*PI/2)*sin(pow(x_var[3],alpha)*PI/2);
		y_obj[2] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*cos(pow(x_var[1],alpha)*PI/2)*sin(pow(x_var[2],alpha)*PI/2);
		y_obj[3] = (1 + g)*cos(pow(x_var[0],alpha)*PI/2)*sin(pow(x_var[1],alpha)*PI/2);
		y_obj[4] = (1 + g)*sin(pow(x_var[0],alpha)*PI/2);

		return;
	}

	if(!strcmp("DTLZ1_6_6",strTestInstance))
	{
		double g = 0;
		for(int n=nobj-1; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);

		y_obj[0] = 0.5*(1 + g)*x_var[0]*x_var[1]*x_var[2]*x_var[3]*x_var[4];
		y_obj[1] = 0.5*(1 + g)*x_var[0]*x_var[1]*x_var[2]*x_var[3]*(1-x_var[4]);
		y_obj[2] = 0.5*(1 + g)*x_var[0]*x_var[1]*x_var[2]*(1-x_var[3]);
		y_obj[3] = 0.5*(1 + g)*x_var[0]*x_var[1]*(1 - x_var[2]);
		y_obj[4] = 0.5*(1 + g)*x_var[0]*(1 - x_var[1]);
		y_obj[5] = 0.5*(1 + g)*(1-x_var[0]);
		return;
	}

	if(!strcmp("DTLZ2_6_6",strTestInstance))
	{
		double g = 0;		
		for(int n=nobj-1; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}
		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*cos(x_var[3]*PI/2)*cos(x_var[4]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*cos(x_var[3]*PI/2)*sin(x_var[4]*PI/2);
		y_obj[2] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*sin(x_var[3]*PI/2);
		y_obj[3] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(x_var[2]*PI/2);
		y_obj[4] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[5] = (1 + g)*sin(x_var[0]*PI/2);
		return;
	}

	if(!strcmp("DTLZ3_6_6",strTestInstance))
	{
		double g = 0;
		for(int n=nobj-1; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);

		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*cos(x_var[3]*PI/2)*cos(x_var[4]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*cos(x_var[3]*PI/2)*sin(x_var[4]*PI/2);
		y_obj[2] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*sin(x_var[3]*PI/2);
		y_obj[3] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(x_var[2]*PI/2);
		y_obj[4] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[5] = (1 + g)*sin(x_var[0]*PI/2);

		return;
	}

	if(!strcmp("DTLZ4_6_6",strTestInstance))
	{
		double g = 0;
		double alpha = 100;
		for(int n=nobj-1; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;;
		}

		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*cos(x_var[3]*PI/2)*cos(x_var[4]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*cos(x_var[3]*PI/2)*sin(x_var[4]*PI/2);
		y_obj[2] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(x_var[2]*PI/2)*sin(x_var[3]*PI/2);
		y_obj[3] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(x_var[2]*PI/2);
		y_obj[4] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[5] = (1 + g)*sin(x_var[0]*PI/2);

		return;
	}

	if(!strcmp("DTLZ1_3_6",strTestInstance))
	{
		double g = 0;
		for(int n=nobj-1; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);
		double temp2=0.5*(1+2*g*x_var[2])/(1+g);
		double temp3=0.5*(1+2*g*x_var[3])/(1+g);
		double temp4=0.5*(1+2*g*x_var[4])/(1+g);
		y_obj[0] = 0.5*(1 + g)*x_var[0]*x_var[1]*temp2*temp3*temp4;
		y_obj[1] = 0.5*(1 + g)*x_var[0]*x_var[1]*temp2*temp3*(1-temp4);
		y_obj[2] = 0.5*(1 + g)*x_var[0]*x_var[1]*temp2*(1-temp3);
		y_obj[3] = 0.5*(1 + g)*x_var[0]*x_var[1]*(1-temp2);
		y_obj[4] = 0.5*(1 + g)*x_var[0]*(1-x_var[1]);
		y_obj[5] = 0.5*(1 + g)*(1-x_var[0]);
		return;
	}

	if(!strcmp("DTLZ3_3_6",strTestInstance))
	{
		double g = 0;
		for(int n=nobj-1; n<nvar;n++)				
			g = g + pow(x_var[n]-0.5,2) - cos(20*PI*(x_var[n] - 0.5));
		g = 100*(nvar + 1 - nobj + g);

		double temp2=PI/4*(1+2*g*x_var[2])/(1+g);
		double temp3=PI/4*(1+2*g*x_var[3])/(1+g);
		double temp4=PI/4*(1+2*g*x_var[4])/(1+g);
		int I=3;
		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*cos(temp3)*cos(temp4);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*cos(temp3)*sin(temp4);
		y_obj[2] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*sin(temp3);
		y_obj[3] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(temp2);
		y_obj[4] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[5] = (1 + g)*sin(x_var[0]*PI/2);
		return;
	}

	if(!strcmp("DTLZ5_3_10",strTestInstance))
	{
		//double g = 0;		
		//for(int n=nobj-1; n<nvar;n++)				
		//{
		//	double x = (x_var[n] - 0.5);
		//	g = g + x*x;;
		//}
		//double temp2=PI/4*(1+2*g*x_var[2])/(1+g);
		//double temp3=PI/4*(1+2*g*x_var[3])/(1+g);
		//double temp4=PI/4*(1+2*g*x_var[4])/(1+g);
		//int I=3;
		//y_obj[0] = (1 + 100*g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*cos(temp3)*cos(temp4);
		//y_obj[1] = (1 + 100*g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*cos(temp3)*sin(temp4);
		//y_obj[2] = (1 + 100*g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*cos(temp2)*sin(temp3);
		//y_obj[3] = (1 + 100*g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2)*sin(temp2);
		//y_obj[4] = (1 + 100*g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		//y_obj[5] = (1 + 100*g)*sin(x_var[0]*PI/2);

		int n,l;
		double g = 0;		
		for(n=nobj-1; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;
		}
		int I=3,M=10;
		double theta[11]; 
		for (n=1; n<=I-1; n++)  theta[n] = PI/2*x_var[n-1];
		for (n=I; n<=M-1; n++)	theta[n] = PI/4*(1+2*g*x_var[n-1])/(1+g);
		double multi = 1.0;
		for (n=1;n<=M;n++)
		{
			multi = 1.0;
			for (l=1;l<=M-n;l++)	multi = multi*cos(theta[l]);
			if (n>1)				multi = multi*sin(theta[M-n+1]);
			y_obj[n-1] = (1 + 100*g) * multi;
		}

		return;
	}

		if(!strcmp("ZDT1_Scale",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);

		y_obj[0] = x_var[0];
		y_obj[1] = 10*g*(1 - sqrt(y_obj[0]/g));
		return;
	}

	if(!strcmp("DTLZ2_Scale",strTestInstance))
	{
		double g = 0;		
		for(int n=2; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;
		}
		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2);
		y_obj[1] = 10*(1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[2] = 100*(1 + g)*sin(x_var[0]*PI/2);
		return;
	}


	if(!strcmp("ZDT1_HSLT",strTestInstance))
	{
		double g = 0;
		for(int n=1;n<nvar;n++)
			g+= x_var[n];
		g = 1 + 9*g/(nvar-1);


		double N = 2.8;
		y_obj[0] = x_var[0];
		y_obj[1] = g*(1-pow(1 - pow(1-x_var[0],N),1.0/N));

		return;
	}

	if(!strcmp("DTLZ2_HSLT",strTestInstance))
	{
		double g = 0;		
		for(int n=2; n<nvar;n++)				
		{
			double x = (x_var[n] - 0.5);
			g = g + x*x;
		}
		y_obj[0] = (1 + g)*cos(x_var[0]*PI/2)*cos(x_var[1]*PI/2);
		y_obj[1] = (1 + g)*cos(x_var[0]*PI/2)*sin(x_var[1]*PI/2);
		y_obj[2] = (1 + g)*sin(x_var[0]*PI/2);


		y_obj[0] = y_obj[0]*y_obj[0]*y_obj[0]*y_obj[0];
		y_obj[1] = y_obj[1]*y_obj[1]*y_obj[1]*y_obj[1];
		y_obj[2] = y_obj[2]*y_obj[2];

		return;
	}

}


void CIndividual::show_objective()
{
    for(int n=0; n<nobj; n++)
		printf("%f ",y_obj[n]);
	printf("\n");
}

void CIndividual::show_variable()
{
    for(int n=0; n<nvar; n++)
		printf("%f ",x_var[n]);
	printf("\n");
}

void CIndividual::operator=(const CIndividual &ind2)
{
    x_var = ind2.x_var;
	y_obj = ind2.y_obj;
	rank  = ind2.rank;
	crowd_distance = ind2.crowd_distance;
}

bool CIndividual::operator<(const CIndividual &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}


bool CIndividual::operator<<(const CIndividual &ind2)
{
	bool dominated = true;
    for(int n=0; n<nobj; n++)
	{
		if(ind2.y_obj[n]<y_obj[n]  - 0.0001) return false;
	}
	if(ind2.y_obj==y_obj) return false;
	return dominated;
}

bool CIndividual::operator==(const CIndividual &ind2)
{
	if(ind2.y_obj==y_obj) return true;
	else return false;
}


// defining subproblem 

class CSubproblem  
{
public:
	CSubproblem();
	virtual ~CSubproblem();

	void show();

	CIndividual     indiv;     // best solution
	CIndividual     saved;     // last solution
	vector <double> namda;     // weight vector
	vector <int>    table;     // neighbourhood table

	double          fitness;

    void  operator=(const CSubproblem &sub2);
};

CSubproblem::CSubproblem()
{
    namda = vector<double>(nobj, 0);
}

CSubproblem::~CSubproblem()
{
	namda.clear();
	table.clear();
}

void CSubproblem::show()
{
   for(int n=0; n<namda.size(); n++)
   {
       printf("%f ",namda[n]);
   }
   printf("\n");
}

void CSubproblem::operator=(const CSubproblem &sub2)
{
    indiv  = sub2.indiv;
	table  = sub2.table;
	namda  = sub2.namda;
}


#endif
