#ifndef __COMMON_H_
#define __COMMON_H_

#include "global.h"

void minfastsort(vector<double> &x, vector<int> &idx, int n, int m)
{
    for(int i=0; i<m; i++)
	{
	    for(int j=i+1; j<n; j++)
			if(x[i]>x[j])
			{
			    double temp = x[i];
				x[i]        = x[j];
				x[j]        = temp;
				int id      = idx[i];
				idx[i]      = idx[j];
				idx[j]      = id;
			}
	}
}

void quick_sort(vector<double> &dist, vector<int> &seq, int num_relest)
{
	int begin=0,stop = dist.size()-1,quick_index = stop+1;
	int start=0,end = dist.size()-1;
	double temp;
	int temp_int;
	while (quick_index > num_neighbour + num_relest || quick_index < num_neighbour)
	{
		if (quick_index > num_neighbour + num_relest) 
		{
			start=begin;
			end = quick_index - 1;
			quick_index=end;
			stop = end;
		}
		else
		{
			start=quick_index +1;
			begin = start;
			end = stop;
			quick_index=end;
		}		
		while (start < end)
		{		
			while (start < end && dist[start] <= dist[quick_index]) start++;
			if (start < end)
			{
				temp = dist[start];
				dist[start] = dist[quick_index];
				dist[quick_index] = temp;
				temp_int = seq[start];
				seq[start] = seq[quick_index];
				seq[quick_index] = temp_int;
				
				quick_index = start;	
			}
				

			while (start < end && dist[end] >= dist[quick_index]) end--;
			if (start < end)
			{
				temp = dist[end];
				dist[end] = dist[quick_index];
				dist[quick_index] = temp;
				temp_int = seq[end];
				seq[end] = seq[quick_index];
				seq[quick_index] = temp_int;

				quick_index = end;
			}
		}
		quick_index = start;	
	}
}

double dist(double *pf,double *pf_find)
{
	double temp  = 0;
	for (int i=0; i<nobj; i++)
	{
		temp += (pf[i]-pf_find[i])*(pf[i]-pf_find[i]);
	}
	return temp;
}

int find_point(double **pf,vector<bool> flag,double **pf_find,int num,int pop)
{
	int i,j;
	double max_distance = -1e9;
	int max_index = -1;
	for (i=0; i<pop; i++)
	{
		if (flag[i] == true)
		{
			double distance_point = 1e10;
			for (j=0; j<num; j++)
			{
				double distance = dist(pf[i],pf_find[j]);
				if (distance <distance_point)
				{
					distance_point = distance;
				}
			}
			if (distance_point > max_distance)
			{
				max_distance = distance_point;
				max_index = i;
			}
		}
	}
	return max_index;
}

int donimate_judge(vector<double> pf1, vector<double> pf2)
{
	//return 1: point j donimate i;return -1: point i donimate j
	double big=0,small=0;
	for (int i=0; i<nobj; i++)
	{
		if (pf1[i]>=pf2[i])	big++;
		if (pf1[i]<=pf2[i]) small++;		
	}
	if (small == nobj)	return -1;
	else if (big == nobj) return 1;
	else return 0;
}

double dist_vector(vector <double> &vec1, vector <double> &vec2)
{
	int dim = vec1.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum+=(vec1[n] - vec2[n])*(vec1[n] - vec2[n]);
	return sqrt(sum);
}

double norm(vector <double> vec)
{
	int dim = vec.size();
    double sum = 0;
	for(int n=0; n<dim; n++)
	    sum += vec[n] * vec[n];
	return sqrt(sum);
}

double innerproduct(vector <double> x, vector <double> y)
{
	int dim = x.size();
	double sum = 0;
	for (int n = 0; n < dim; n++)	sum = sum + x[n] * y[n];
	return sum;
}

bool   dominate(vector<double> &u, vector<double> &v, double epsilon)
{
    int dim = u.size();
	for(int i=0; i<dim; i++)
	{
	    if(u[i]<v[i]-epsilon)
		{
		    return false;
		}
	}
	return true;
}

void load_data(char filename[1024], vector< vector<double> > &data, int dim)
{
	std::ifstream readf(filename);
	vector<double> vect = vector<double>(dim, 0);
	while(!readf.eof())
	{
        for(int i=0; i<dim; i++)
		{
			readf>>vect[i];
			//printf(" %f ", vect[i]);
		}
		data.push_back(vect);
		//vect.clear();    // bug here. 
	}
	vect.clear();
	readf.close();    
}

//----------------------------------------------//
/* the following function designed for  WFG*/
//----------------------------------------------//
// the following function designed for shape for object-function//
double Convex(vector<double> x_var, int m)
{
	double result = 1;
	int i;
	
	for (i = 1; i <= nobj - m; i++)		result *= 1 - cos(x_var[i-1]*PI/2); 
	if (m == 1)			return result;
	else				return result *= 1-sin(x_var[nobj - m + 1 - 1]*PI/2);		
}

double Linear(vector<double> x_var, int m)
{
	double result = 1;
	int i;
	
	for (i = 1; i <= nobj - m; i++)		result *= x_var[i-1]; 
	if (m == 1)			return result;
	else				return result *= 1-x_var[nobj - m + 1 - 1];
}

double Concave(vector<double> x_var, int m)
{
	double result = 1;
	int i;
	
	for (i = 1; i <= nobj - m; i++)		result *= sin(x_var[i-1]*PI/2); 
	if (m == 1)			return result;
	else				return result *= cos(x_var[nobj - m + 1 - 1]*PI/2);		
}

double Mixed(vector<double> x_var, double alfa, int A)
{//alfa > 0,A subject to {1,2,3,...}
	return pow(1 - x_var[0] - (cos(2*A*PI + PI/2)) / (2*A*PI), alfa);	
}

double Disc(vector<double> x_var, double alfa, double belta, int A)
{//alfa > 0, belta > 0, A subject to {1,2,3,...}
	return 1 - pow(x_var[0], alfa) * cos(A * pow(x_var[0], belta) * PI)* cos(A * pow(x_var[0], belta) * PI);	
}

//----------------------------------------------//
// the following function designed for Transformation functions//
double b_poly(double y, double alfa)
{//alfa > 0, alfa not equal to 1
	return pow(y, alfa);	
}

double b_flat(double y, double A, double B, double C)
{//A,B,C subject to [0,1]
	double temp1,temp2;
	if (floor(y - B) <= 0)	temp1 = floor(y - B);
	else					temp1 = 0;
	if (floor(C - y) <= 0)	temp2 = floor(C - y);
	else					temp2 = 0;
	return A + temp1 * A* (B - y) / B - temp2 * (1 - A) * (y - C) / (1 - C);	
}

double b_param(double y, double u_y, double A, double B, double C)
{//A subject to [0,1], 0 < B < C	
	return pow(y, B + (C - B) * (A - (1 - 2 * u_y) * fabs(floor(0.5 - u_y) + A) ));	
}

double s_linear(double y, double A)
{//A subject to [0,1]	
	return fabs(y - A) / fabs(floor(A - y) + A);	
}

double s_multi(double y, int A, double B, double C)
{//A subject to {1,2,...}, B>=0, (4*A+2)*PI>=4*B, C subject to [0,1]	
	return (1 + cos((4*A + 2) * PI * (0.5 - fabs(y-C)/(2*(floor(C - y) + C)))) + 4 * B * fabs(y-C)/(2*(floor(C-y)+C)) *fabs(y-C)/(2*(floor(C-y)+C)) ) / (B + 2);	
}

double s_decept(double y, double A, double B, double C)
{//A subject to [0,1], 0<=B<<1, 0<=C<<1, A-B>0,A+B<1	
	return 1+(fabs(y-A)-B)*((floor(y-A+B)*(1-C+(A-B)/B))/(A-B)+(floor(A+B-y)*(1-C+(1-A-B)/B))/(1-A-B)+1.0/B);
}

double r_sum(vector<double> y, vector<double> w)
{//|w| = |y|, w1,w2,...w|y| > 0
	double total_w = 0;
	double weight_sum = 0;
	int i;
	for (i = 0; i < w.size(); i++)
	{
		total_w += w[i];
		weight_sum += w[i] * y[i];
	}
	return weight_sum/total_w;	
}

double r_nonsep(vector<double> y, double A)
{//A subject to {1,...,|y|}, |y| mod A = 0
	double numerator = 0.0;
	
	for( int j = 0; j < y.size(); j++ )
	{
		numerator += y[j];

		for( int k = 0; k <= A-2; k++ )
		{
			numerator += fabs( y[j] - y[( j+k+1 ) % y.size()] );
		}
	}

	const double tmp = ceil( A/2.0 );
	const double denominator = y.size()*tmp*( 1.0 + 2.0*A - 2.0*tmp )/A;

	return  numerator / denominator;
}

vector<double>  WFG_normalise_z(vector<double> x)
{	
	int i;
	vector<double> z(nvar, 0);
	for (i = 1; i <= x.size(); i++)		z[i-1] = x[i-1]*1.0 / (2*i);	
	return z;
}

vector<double> WFG1_t1(vector<double> y_old, int k)
{	
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i <= k; i++)	y_new[i-1] = y_old[i-1];
	for (i = k + 1; i <= nvar; i++)		y_new[i-1] = s_linear(y_old[i-1],0.35);
	
	return y_new;
}

vector<double> WFG1_t2(vector<double> y_old, int k)
{	
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i <= k; i++)	y_new[i-1] = y_old[i-1];
	for (i = k + 1; i <= nvar; i++)		y_new[i-1] = b_flat(y_old[i-1],0.8,0.75,0.85);
	
	return y_new;
}

vector<double> WFG1_t3(vector<double> y_old)
{	
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i <= nvar; i++)	y_new[i-1] = b_poly(y_old[i-1],0.02);	
	
	return y_new;
}

vector<double> WFG1_t3_var(vector<double> y_old)
{	
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i <= nvar; i++)	y_new[i-1] = b_poly(y_old[i-1],0.8);	
	
	return y_new;
}

vector<double> WFG1_t4(vector<double> y_old, int k)
{	
	int i,j;
	vector<double> y_new(nobj, 0);

	vector<double> temp_y;
	vector<double> w;

	for (i = 1; i <= nobj - 1; i++)	
	{
		for (j = (i-1) *k /(nobj -1)+1; j <= i * k /(nobj -1); j++)
		{
			temp_y.push_back(y_old[j-1]);
			w.push_back(2*j);
		}
		y_new[i-1] = r_sum(temp_y, w);
		
		temp_y.clear();
		w.clear();
	}
	for (j = k+1; j < nvar; j++)
	{
		temp_y.push_back(y_old[j-1]);
		w.push_back(2*j);
	}
	y_new[nobj - 1] = r_sum(temp_y, w);

	return y_new;
}

vector<double> WFG1_shape(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nobj, 0);
	vector<double> x(nobj, 0);
	
	double max;
	if (y_old[nobj - 1] < 1)	max = 1;
	else						max = y_old[nobj - 1];
	//turn to x
	for (i = 1; i < nobj; i++)	x[i-1] = max*(y_old[i-1] - 0.5)+0.5;
	x[nobj - 1] = y_old[nobj - 1];
	//shape
	for (i = 1; i < nobj; i++)	y_new[i-1] = x[nobj - 1] + 2*i *Convex(x,i);
	y_new[nobj - 1] = x[nobj - 1] + 2*nobj *Mixed(x,1,5);

	return y_new;
}
vector<double> WFG1_VAR_shape(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nobj, 0);
	vector<double> x(nobj, 0);
	
	double max;
	if (y_old[nobj - 1] < 1)	max = 1;
	else						max = y_old[nobj - 1];
	//turn to x
	for (i = 1; i < nobj; i++)	x[i-1] = max*(y_old[i-1] - 0.5)+0.5;
	x[nobj - 1] = y_old[nobj - 1];
	//shape
	for (i = 1; i < nobj; i++)	y_new[i-1] = x[nobj - 1] + 2*i *Convex(x,i);
	y_new[nobj - 1] = x[nobj - 1] + 2*nobj *Mixed(x,7,8);

	return y_new;
}


vector<double> WFG2_t2(vector<double> y_old, int k)
{	
	int i,j;
	vector<double> y_new;
	int l = nvar - k;
	
	for (i = 1; i <= k ; i++) y_new.push_back( y_old[i-1] );
	for (i = k + 1;  i <= k+l/2; i++)
	{
		vector<double> temp;
		temp.push_back(y_old[k+2*(i-k) - 1-1]);
		temp.push_back(y_old[k+2*(i-k)-1]);
		y_new.push_back(r_nonsep(temp, 2));

		temp.clear();
	}

	return y_new;
}

vector<double> WFG2_t3(vector<double> y_old, int k)
{	
	int i,j;
	int l = nvar - k;
	vector<double> y_new(nobj, 0);

	vector<double> temp_y;
	vector<double> w;

	for (i = 1; i <= nobj - 1; i++)	
	{
		for (j = (i-1) *k /(nobj -1)+1; j <= i * k /(nobj -1); j++)
		{
			temp_y.push_back(y_old[j-1]);
			w.push_back(1);
		}
		y_new[i-1] = r_sum(temp_y, w);
		
		temp_y.clear();
		w.clear();
	}
	for (j = k+1; j <= k + l/2; j++)
	{
		temp_y.push_back(y_old[j-1]);
		w.push_back(1);
	}
	y_new[nobj - 1] = r_sum(temp_y, w);

	return y_new;
}

vector<double> WFG2_shape(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nobj, 0);
	vector<double> x(nobj, 0);
	
	double max;
	if (y_old[nobj - 1] < 1)	max = 1;
	else						max = y_old[nobj - 1];
	//turn to x
	for (i = 1; i < nobj; i++)	x[i-1] = max*(y_old[i-1] - 0.5)+0.5;
	x[nobj - 1] = y_old[nobj - 1];
	//shape
	for (i = 1; i < nobj; i++)	y_new[i-1] = x[nobj - 1] + 2*i *Convex(x,i);
	y_new[nobj - 1] = x[nobj - 1] + 2*nobj *Disc(x,1,1,5);

	return y_new;
}

vector<double> WFG3_shape(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nobj, 0);
	vector<double> x(nobj, 0);
	vector<double> A(nobj, 0);
	A[0] = 1;
	
	//turn to x
	for (i = 1; i < nobj; i++)
	{
		double max;
		if (y_old[nobj - 1] < A[i-1])	max = 1;
		else						max = y_old[nobj - 1];
		x[i-1] = max*(y_old[i-1] - 0.5)+0.5;
	}
	x[nobj - 1] = y_old[nobj - 1];
	//shape
	for (i = 1; i <= nobj; i++)	y_new[i-1] = x[nobj - 1] + 2*i *Linear(x,i);

	return y_new;
}

vector<double> WFG4_t1(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i < nvar; i++) y_new[i - 1] = s_multi(y_old[i - 1], 30, 10, 0.35);

	return y_new;
}

vector<double> WFG4_shape(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nobj, 0);
	vector<double> x(nobj, 0);
	vector<double> A(nobj, 1);
	
	//turn to x
	for (i = 1; i < nobj; i++)
	{
		double max;
		if (y_old[nobj - 1] < A[i-1])	max = A[i-1];
		else							max = y_old[nobj - 1];
		x[i-1] = max*(y_old[i-1] - 0.5)+0.5;
	}
	x[nobj - 1] = y_old[nobj - 1];
	//shape
	for (i = 1; i <= nobj; i++)	y_new[i-1] = x[nobj - 1] + 2*i *Concave(x,i);

	return y_new;
}

vector<double> WFG5_t1(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i < nvar; i++) y_new[i - 1] = s_decept(y_old[i - 1], 0.35, 0.001, 0.05);

	return y_new;
}

vector<double> WFG6_t2(vector<double> y_old, int k)
{
	int i,j;
	vector<double> y_new(nobj, 0);

	vector<double> temp_y;

	for (i = 1; i <= nobj - 1; i++)	
	{
		for (j = (i-1) *k /(nobj -1)+1; j <= i * k /(nobj -1); j++)
		{
			temp_y.push_back(y_old[j-1]);
		}
		y_new[i-1] = r_nonsep(temp_y, k*1.0/(nobj - 1));
		
		temp_y.clear();
	}
	for (j = k+1; j <= nvar; j++)		temp_y.push_back(y_old[j-1]);

	y_new[nobj - 1] = r_nonsep(temp_y, nvar - k);

	return y_new;
}

vector<double> WFG7_t1(vector<double> y_old, int k)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	vector<double> temp_y;
	vector<double> w;
	for (i = 1; i <= k; i++)	
	{
		for (j = i + 1; j <= nvar; j++)
		{
			temp_y.push_back(y_old[j-1]);
			w.push_back(1);
		}
		y_new[i-1] = b_param(y_old[i-1], r_sum(temp_y, w), 0.98/49.98, 0.02, 50);
		
		temp_y.clear();
		w.clear();
	}
	for (j = k+1; j <= nvar; j++)		y_new[j - 1] = y_old[j-1];	

	return y_new;
}

vector<double> WFG8_t1(vector<double> y_old, int k)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	vector<double> temp_y;
	vector<double> w;
	
	for (j = 1; j <= k; j++)		y_new[j - 1] = y_old[j-1];
	for (i = k+1; i <= nvar; i++)	
	{
		for (j = 1; j <= i - 1; j++)
		{
			temp_y.push_back(y_old[j-1]);
			w.push_back(1);
		}
		y_new[i-1] = b_param(y_old[i-1], r_sum(temp_y, w), 0.98/49.98, 0.02, 50);
		
		temp_y.clear();
		w.clear();
	}	

	return y_new;
}

vector<double> WFG9_t1(vector<double> y_old)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	vector<double> temp_y;
	vector<double> w;
	for (i = 1; i <= nvar - 1; i++)	
	{
		for (j = i + 1; j <= nvar; j++)
		{
			temp_y.push_back(y_old[j-1]);
			w.push_back(1);
		}
		y_new[i-1] = b_param(y_old[i-1], r_sum(temp_y, w), 0.98/49.98, 0.02, 50);
		
		temp_y.clear();
		w.clear();
	}
	y_new[nvar - 1] = y_old[nvar - 1];	

	return y_new;
}

vector<double> WFG9_t2(vector<double> y_old, int k)
{
	int i,j;
	vector<double> y_new(nvar, 0);

	for (i = 1; i <= k; i++)		y_new[i - 1] = s_decept(y_old[ i - 1 ], 0.35, 0.001, 0.05);
	for (i = k + 1; i <= nvar; i++)	y_new[i - 1] = s_multi(y_old[i - 1],30, 95, 0.35);	

	return y_new;
}

#endif