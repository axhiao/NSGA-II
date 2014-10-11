#ifndef _PROBLEM_H
#define _PROBLEM_H

#include "cec09.h"/*

// *********************** CEC 2009 ************************************

void CEC09_F1(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF1(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F2(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF2(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F3(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	CEC09::UF3(&(*(X.begin())), &(*(F.begin())), X.size());
}
void CEC09_F4(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF4(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF5(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F6(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF6(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F7(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(2); 
	std::vector<double> XX(X);
	for(unsigned int i=1; i<X.size(); i++) XX[i] = -1.0 + 2.0*X[i];
	CEC09::UF7(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F8(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF8(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F9(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF9(&(*(XX.begin())), &(*(F.begin())), X.size());
}
void CEC09_F10(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(3); 
	std::vector<double> XX(X);
	for(unsigned int i=2; i<X.size(); i++) XX[i] = -2.0 + 4.0*X[i];
	CEC09::UF10(&(*(XX.begin())), &(*(F.begin())), X.size());
}

// ---------------- 5 objective test instances ------------------

void CEC09_R2_DTLZ2_M5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(5); 
	std::vector<double> XX(X);
	
	double low30[] = {-1.773,	-1.846,	-1.053,	-2.370,	-1.603,	-1.878,	-1.677,	-0.935,	-1.891,	-0.964,	-0.885,	-1.690,	-2.235,	-1.541,	-0.720,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000};
	double upp30[] = {1.403,	1.562,	2.009,	0.976,	1.490,	1.334,	1.074,	2.354,	1.462,	2.372,	2.267,	1.309,	0.842,	1.665,	2.476,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000};

	for(unsigned int i=0; i<X.size(); i++) XX[i] = low30[i] + X[i]*(upp30[i] - low30[i]);
	CEC09::R2_DTLZ2_M5(&(*(XX.begin())), &(*(F.begin())), X.size(), nobj);
}

void CEC09_R2_DTLZ3_M5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(5); 
	std::vector<double> XX(X);
	
	double low30[] = {-1.773,	-1.846,	-1.053,	-2.370,	-1.603,	-1.878,	-1.677,	-0.935,	-1.891,	-0.964,	-0.885,	-1.690,	-2.235,	-1.541,	-0.720,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000,	0.000};
	double upp30[] = {1.403,	1.562,	2.009,	0.976,	1.490,	1.334,	1.074,	2.354,	1.462,	2.372,	2.267,	1.309,	0.842,	1.665,	2.476,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000,	1.000};

	for(unsigned int i=0; i<X.size(); i++) XX[i] = low30[i] + X[i]*(upp30[i] - low30[i]);
	CEC09::R2_DTLZ3_M5(&(*(XX.begin())), &(*(F.begin())), X.size(), nobj);
}

void CEC09_WFG1_M5(std::vector< double >& F, std::vector< double >& X)
{
	F.resize(5); 
	std::vector<double> XX(X);
	CEC09::WFG1_M5(&(*(XX.begin())), &(*(F.begin())), X.size(), nobj);
}
*/
#endif


  