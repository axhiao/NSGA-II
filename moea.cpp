/*==========================================================================
//  C++ Implementation of NSGAII for Contest Multiobjective Problems
//
//  Author: Xiaoliang Ma
//
//  See the details of NSGAII and test problems in the following papers
//
//  1) K. Deb, S. Agrawal, A. Pratap, and T. Meyarivan, ¡°A fast and elitist multiobjective genetic algorithm:
       NSGA-II,¡± IEEE Trans. Evol.Comput., vol. 6, no. 2, pp. 182¨C197, Apr. 2002.
//
//  If you have any questions about the codes, please contact
//
//  Xiaoliang Ma at   maxiaoliang@yeah.net
//
//  Date: 21/DEC/2009
//
// ===========================================================================*/
#include "algorithm.h"

void main()
{	
	std::ifstream readf("TestInstance.txt");//	TestInstance.txt store the 
											//	parameter setting of test instance
	strcpy(crossoverType, "SBX");			//	the type of crossover operater
	
	strcpy(mutationType, "ploynomial");		//	theoperater
	//--- set needed test case index from begin to end --- type of mutation 
	int fun_start_num = 23;					//	the start index of needed test case
	int fun_end_num =23;					//	the end index of needed test case
	int num_run = 30;						//	the number of run for each test case
	//--- skip the front test instance ---
	int seq;
	for (int i =1;i<fun_start_num;i++)
	{	//--- the parameter setting of test instance ---
		readf>>seq;
		readf>>strTestInstance;				//	the name of test case
		readf>>nvar;						//	the number of variable for the test case
		readf>>nobj;						//	the number of objetive for the test case
	}
	//--- running from the fun_start_num case to the fun_end_num case ---
	for(inst=fun_start_num; inst<=fun_end_num; inst++)
	{	//--- the parameter setting of test instance ---
		readf>>seq;
		readf>>strTestInstance;				//	the name of test case
		readf>>nvar;						//	the number of variable for the test case
		readf>>nobj;						//	the number of objetive for the test case

		getVarBoundy();						//	get  the domain of x_var

		if (inst >= 11 && inst <= 23||inst >= 34 && inst <= 40)
		{
			if (nobj == 2)
			{
				pop=100;
				num_max_evulation = 50000;
				strcpy(crossoverType,"SBX");
				//pops = 100;
			}
			else if (nobj == 3)
			{
				pop=300;
				num_max_evulation = 75000;//  DTLZ6:23
				strcpy(crossoverType,"SBX");
				//pops = 100;
			}
		}
		else if (inst >= 1 && inst <= 10 || inst >= 24 && inst <= 33|| inst >= 41 && inst <= 133)
		{
			if (nobj == 2)
			{
				pop = 300;
				num_neighbour = int(pop * 0.1);
				num_max_evulation = 60000;
				strcpy(crossoverType,"DE");
				//pops = 600;
			}
			else if (nobj == 3)
			{
				pop = 600;
				num_neighbour = int(pop * 0.1);
				num_max_evulation = 300000;
				strcpy(crossoverType,"DE");
				//pops = 1000;
			}
		}
		else if (inst >= 134 && inst <= 144)
		{
			if (nobj == 2)
			{
				pop=50;
				num_max_evulation = 30000;
				strcpy(crossoverType,"SBX");
				//pops = 100;
			}
			else if (nobj == 3)
			{
				pop=100;
				num_max_evulation = 30000;
				strcpy(crossoverType,"SBX");
			}
			else if (nobj == 4)
			{
				pop=120;
				num_max_evulation = 30000;
				strcpy(crossoverType,"SBX");
				//pops = 100;
			}
			else if (nobj == 6)
			{
				pop=126;
				num_max_evulation = 30000;
				strcpy(crossoverType,"SBX");
			}
		}
		else if (inst >= 145 && inst <= 145)
		{//ZDT_4
			pop=220;
			num_max_evulation = 200000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst >= 146 && inst <= 150)
		{//DTLZ_5
			pop=252;
			num_max_evulation = 200000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst >= 151 && inst <= 160)
		{//ZDT_3
			num_max_evulation = 50000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst >= 161 && inst <= 166)
		{
			pop = 100;
			num_max_evulation = 30000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst >= 167 && inst <= 170)
		{
			pop = 252;
			num_max_evulation = 200000;
			strcpy(crossoverType,"SBX");
		}	
		else if (inst >= 171 && inst <= 174)
		{
			pop = 220;
			num_max_evulation = 130000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst >= 175 && inst <= 178)
		{
			pop = 495;
			num_max_evulation = 180000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst >= 179 && inst <= 182)
		{
			pop = 792;
			num_max_evulation = 250000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst >= 183 && inst <= 184)
		{
			pop = 252;
			num_max_evulation = 350000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst == 185)
		{
			pop = 220;
			num_max_evulation = 350000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst == 186)
		{
			pop = 100;
			num_max_evulation = 50000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst == 187)
		{
			pop = 300;
			num_max_evulation = 100000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst == 188)
		{
			pop = 100;
			num_max_evulation = 50000;
			strcpy(crossoverType,"SBX");
		}
		else if (inst == 189)
		{
			pop = 300;
			num_max_evulation = 100000;
			strcpy(crossoverType,"SBX");
		}


	
		CMOEA MOEA(pop);							//

		printf("\n -- Instance: %s, %d variables %d objectives \n\n", strTestInstance, nvar, nobj);		
		for(int run=1; run<=num_run; run++)
		{
			printf("\n -- %d-th run  -- \n", run);			
			MOEA.exec_emo(run);				//
		}

		//lowBound.clear();
		//uppBound.clear();
	}
	//--- release the store resource in heap ---
	lowBound.clear();
	uppBound.clear();
}
