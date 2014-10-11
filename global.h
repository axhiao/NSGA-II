#ifndef __GLOBAL_H_
#define __GLOBAL_H_

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <memory.h>
#include <vector>
#include <cassert>
#include <algorithm>

using namespace std;

#include "random.h"

//------------- Parameters for max number ------------------
#define max_int 2147483647				//	the max sign integer number
#define inf     1e10					//	the approminate max double number
#define PI		3.1415926535897932384626433832795

//------------- Parameters in test instance ------------------
int     nvar,  nobj;                    //  the number of variables and objectives
vector<double>  lowBound,   uppBound;   //  lower and upper bounds of variables
char    strTestInstance[256];			//	the string of test case
int		num_max_evulation;				//	the max number of function evulation
int		pop;							//	the size of population 
int		num_neighbour = 10;
int		num_store = 300;
double	emxitong = 1e-3;

//------------- Parameters in random number ------------------
int     seed    = 177;					//	the initialization seed
long    rnd_uni_init;					//	

//------------- Parameters in evolution operator ------------------
char	crossoverType[1024];
char	mutationType[1024];
int		etax    = 20, 	etam    = 20;   // distribution indexes of crossover and mutation
double  realx,  realm,  realb = 0.9;    // crossover, mutation, selection probabilities


int num_of_random =10;
int position_parameters = 2;
int distance_parameters = 2;
int inst=1;

//---------------------------Parameters in nondonimate sorting---------------
vector< vector<int> > front;
vector< vector< int > > donimate;
vector<int>	num_donimated;


void getVarBoundy(void)
{
	int i =0;
	if(!strcmp("UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);	
		return;
	}
	if(!strcmp("UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("UF8", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=2;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("UF9", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=2;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("UF10", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=2;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if (!strcmp("SCH", strTestInstance))
	{
		lowBound = vector<double>(nvar, -5);
		uppBound = vector<double>(nvar, 10);
		return;
	}
	if (!strcmp("DEB", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("KUR", strTestInstance))
	{
		lowBound = vector<double>(nvar, -5);
		uppBound = vector<double>(nvar, 5);
		return;
	}
	if (!strcmp("ZDT1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("ZDT2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("ZDT3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("ZDT4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=2; i<=nvar; i++)
		{
			lowBound[i-1] = -5;
			uppBound[i-1] = 5;
		}
		return;
	}
	if (!strcmp("ZDT6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	
	if(!strcmp("MOP1", strTestInstance))
	{
		lowBound = vector<double>(nvar, -100000);
		uppBound = vector<double>(nvar, 100000);
		return;
	}
	if(!strcmp("MOP2", strTestInstance))
	{
		lowBound = vector<double>(nvar, -4);
		uppBound = vector<double>(nvar, 4);
		return;
	}if(!strcmp("MOP3", strTestInstance))
	{
		lowBound = vector<double>(nvar, -PI);
		uppBound = vector<double>(nvar, PI);
		return;
	}
	if(!strcmp("MOP4", strTestInstance))
	{
		lowBound = vector<double>(nvar, -5);
		uppBound = vector<double>(nvar, 5);
		return;
	}if(!strcmp("MOP5", strTestInstance))
	{
		lowBound = vector<double>(nvar, -30);
		uppBound = vector<double>(nvar, 30);
		return;
	}if(!strcmp("MOP6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if(!strcmp("MOP7", strTestInstance))
	{
		lowBound = vector<double>(nvar, -400);
		uppBound = vector<double>(nvar, 400);
		return;
	}
	
	if(!strcmp("CF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("WFG1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG7", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG8", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
	if(!strcmp("WFG9", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i < nvar; i++)		uppBound[i] = 2*(i + 1);
		return;
	}
		if(!strcmp("SCH_VAR_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 10;
		return;
	}
	if(!strcmp("SCH_VAR_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 10;
		return;
	}
	if(!strcmp("SCH_VAR_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 10;
		return;
	}
	if(!strcmp("SCH_VAR_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, -2);
		uppBound = vector<double>(nvar, 2);
		lowBound[0] = 0;
		uppBound[0] = 10;
		return;
	}
	if(!strcmp("SCH_VAR_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 10;
		return;
	}
	if(!strcmp("SCH_VAR_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 10;
		return;
	}
	if(!strcmp("SCH_VAR_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 10;
		return;
	}
	if(!strcmp("DEB_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		lowBound[1] = 0;
		return;
	}
	if(!strcmp("DEB_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		lowBound[1] = 0;
		return;
	}
	if(!strcmp("DEB_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if(!strcmp("DEB_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, -2);
		uppBound = vector<double>(nvar, 2);
		lowBound[0] = 0;uppBound[0] = 1;
		lowBound[1] = 0;uppBound[1] = 1;
		return;
	}
	if(!strcmp("DEB_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		lowBound[1] = 0;
		return;
	}
	if(!strcmp("DEB_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		lowBound[1] = 0;
		return;
	}
	if(!strcmp("DEB_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		lowBound[1] = 0;
		return;
	}
	if(!strcmp("KUR_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i = 1; i <= 3; i++)
		{
			lowBound[i-1] = -5;
			uppBound[i-1] = 5;
		}
		return;
	}
	if(!strcmp("KUR_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i = 1; i <= 3; i++)
		{
			lowBound[i-1] = -5;
			uppBound[i-1] = 5;
		}
		return;
	}
	if(!strcmp("KUR_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i = 1; i <= 3; i++)
		{
			lowBound[i-1] = -5;
			uppBound[i-1] = 5;
		}
		return;
	}
	if(!strcmp("KUR_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, -2);
		uppBound = vector<double>(nvar, 2);
		for (i = 1; i <= 3; i++)
		{
			lowBound[i-1] = -5;
			uppBound[i-1] = 5;
		}
		return;
	}
	if(!strcmp("KUR_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i = 1; i <= 3; i++)
		{
			lowBound[i-1] = -5;
			uppBound[i-1] = 5;
		}
		return;
	}
	if(!strcmp("KUR_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i = 1; i <= 3; i++)
		{
			lowBound[i-1] = -5;
			uppBound[i-1] = 5;
		}
		return;
	}
	if(!strcmp("KUR_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i = 1; i <= 3; i++)
		{
			lowBound[i-1] = -5;
			uppBound[i-1] = 5;
		}
		return;
	}
	if(!strcmp("ZDT3_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;		
		return;
	}
	if(!strcmp("ZDT3_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;		
		return;
	}
	if(!strcmp("ZDT3_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;		
		return;
	}
	if(!strcmp("ZDT3_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, -2);
		uppBound = vector<double>(nvar, 2);
		lowBound[0] = 0;
		uppBound[0] = 1;		
		return;
	}
	if(!strcmp("ZDT3_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;		
		return;
	}
	if(!strcmp("ZDT3_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;		
		return;
	}
	if(!strcmp("ZDT3_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;		
		return;
	}
	if(!strcmp("DTLZ6_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;
		return;
	}
	if(!strcmp("DTLZ6_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;
		return;
	}
	if(!strcmp("DTLZ6_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;
		return;
	}
	if(!strcmp("DTLZ6_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, -2);
		uppBound = vector<double>(nvar, 2);
		lowBound[0] = 0;
		uppBound[0] = 1;
		return;
	}
	if(!strcmp("DTLZ6_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;
		return;
	}
	if(!strcmp("DTLZ6_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;
		return;
	}
	if(!strcmp("DTLZ6_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;
		return;
	}
	if(!strcmp("WFG1_VAR_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG1_VAR_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG1_VAR_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);

		return;
	}
	if(!strcmp("WFG1_VAR_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, -2);
		uppBound = vector<double>(nvar, 2);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG1_VAR_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG1_VAR_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG1_VAR_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG2_VAR_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG2_VAR_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG2_VAR_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);

		return;
	}
	if(!strcmp("WFG2_VAR_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, -2);
		uppBound = vector<double>(nvar, 2);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG2_VAR_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG2_VAR_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("WFG2_VAR_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		lowBound[0] = 0;
		uppBound[0] = 1;

		return;
	}
	if(!strcmp("SQRT_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("SQRT_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("SQRT_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);	
		return;
	}
	if(!strcmp("SQRT_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("SQRT_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("SQRT_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("SQRT_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("FONSICA_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("FONSICA_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("FONSICA_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);	
		lowBound[0] = -1;
		uppBound[0] = 1;
		return;
	}
	if(!strcmp("FONSICA_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("FONSICA_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("FONSICA_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("FONSICA_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, -1);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("LZDTZ3_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("LZDTZ3_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("LZDTZ3_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);	
		return;
	}
	if(!strcmp("LZDTZ3_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("LZDTZ3_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("LZDTZ3_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("LZDTZ3_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("WFG1_UF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("WFG1_UF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("WFG1_UF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);	
		return;
	}
	if(!strcmp("WFG1_UF4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -2;
			uppBound[i] = 2;
		}
		return;
	}
	if(!strcmp("WFG1_UF5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("WFG1_UF6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("WFG1_UF7", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=1;i<nvar;i++)
		{
			lowBound[i] = -1;
			uppBound[i] = 1;
		}
		return;
	}
	if(!strcmp("AB", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABC", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABCD13", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABCD1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABCDNONCON1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABCDNONCON2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABCDEF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABCDEF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABCDEF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	/*
	if(!strcmp("ABCD22", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABCDEF1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABCDEF2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}
	if(!strcmp("ABCDEF3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 100);

		return;
	}*/
	if (!strcmp("ZDT1_4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("ZDT2_4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("ZDT3_4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("ZDT4_4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		for (i=2; i<=nvar; i++)
		{
			lowBound[i-1] = -5;
			uppBound[i-1] = 5;
		}
		return;
	}
	if (!strcmp("ZDT6_4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ1_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ2_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ3_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ4_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ6_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("TDY1", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1.5);
		return;
	}
	if (!strcmp("TDY2", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1.46);
		return;
	}
	if (!strcmp("TDY3", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		return;
	}
	if (!strcmp("TDY4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 2);
		return;
	}
	if (!strcmp("TDY5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 3);
		lowBound[0] = 0.6;
		uppBound[0] = 4.6;
		return;
	}
	if (!strcmp("TDY6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 3);
		lowBound[0] = 0.7;
		uppBound[0] = 4.6;
		return;
	}
	if (!strcmp("DTLZ1_2_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ2_2_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ3_2_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ4_2_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ1_4_4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ2_4_4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ3_4_4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ4_4_4", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ1_5_5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ2_5_5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ3_5_5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ4_5_5", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ1_6_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ2_6_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ3_6_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ4_6_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ1_3_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ3_3_6", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ5_3_10", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("ZDT1_Scale", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ2_Scale", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("ZDT1_HSLT", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
	if (!strcmp("DTLZ2_HSLT", strTestInstance))
	{
		lowBound = vector<double>(nvar, 0);
		uppBound = vector<double>(nvar, 1);
		return;
	}
}
#endif