#ifndef __EVOLUTION_H_
#define __EVOLUTION_H_

#include "global.h"
#include "recomb.h"
#include "common.h"
#include "individual.h"
#include <vector>

using namespace std;

class CMOEA
{
public:
	CMOEA();
	CMOEA(int num_pop);
	virtual ~CMOEA();

	void init_population();                  // initialize the population
	void evol_population();                                      // DE-based recombination
	void mate_selection(vector<int> &list, int cid, int size, int type);  // select mating parents
	// execute MOEAD
	void exec_emo(int run);
	void tour_selection(int depth, vector<int> &selected);
	void fitnessAssign_fast_nondonimate_sort(vector <CIndividual> &population);
	//void fitnessAssign_fast_nondonimate_sort(vector <CIndividual> &population_union, vector< vector<int> > &front);
	void nondonimate_crowd_compare(vector <CIndividual> &population);
	void crowd_distance_assignment(vector <CIndividual> &population, int k);
	//void selection_crowd_compare(vector <CIndividual> population_union, vector <CIndividual> &population_new);
	void crossover_sbx(vector <CIndividual> &population_select);
	void crossover_de();
	void mutation_polynomial(vector <CIndividual> &population_select);
	void evulation(vector <CIndividual> &population_select);	
	void update_store();
	void save_front(char savefilename[1024]);       // save the pareto front into files
	void save_front_select(char savefilename[1024]);       // save the pareto front into files
	void save_pos(char savefilename[1024]);
	void save_population_union(vector <CIndividual> population_union, char saveFilename[1024]);
	void save_individual(char saveFilename_x[1024],char saveFilename_y[1024],vector<double> x,vector<double> y);
	bool boundy_checkout(CIndividual child);
public:
	// algorithm parameters
    int		max_gen;       //  the maximal number of generations
	int     pops;    //  the population size 
	int     nfes;          //  the number of function evluations
	vector <CIndividual> population_evol;
	vector <CIndividual> population_new;
	vector <CIndividual> population_store;
	vector <CIndividual> population_union;



};

CMOEA::CMOEA(int num_pop)
{
	pops = num_pop;	
}

CMOEA::CMOEA()
{
	if (nobj == 2)		pops = 100;
	else if (nobj == 3)		pops = 300;
}

CMOEA::~CMOEA()
{
}

void CMOEA::init_population()
{
	int i,j;
    for(i=0; i<pops; i++)
	{
		CIndividual sub;
		// Randomize and evaluate solution
		sub.rnd_init();
		sub.obj_eval();	
		// Save in the population
		population_evol.push_back(sub);
		nfes++;
	}
	population_new = vector<CIndividual>(pops);
}

void CMOEA::fitnessAssign_fast_nondonimate_sort(vector <CIndividual> &population)
{
	int i,j,count;
	//clear;
	for (i = 0; i < front.size(); i++)		front[i].clear();
	front.clear();
	num_donimated.clear();
	for (i = 0; i < donimate.size(); i++)		donimate[i].clear();
	donimate.clear();

	for (i = 0; i < population.size(); i++)
	{
		front.push_back(vector <int> (0));
		donimate.push_back(vector< int >(0));
		num_donimated.push_back(0);
	}

	int temp;
	for (i = 0; i < population.size(); i++)
	{
		 for (j = i + 1; j < population.size(); j++)
		 {
			 temp = donimate_judge(population[i].y_obj, population[j].y_obj);
			 if (temp == -1) //i donimate j
			 {
				 donimate[i].push_back(j);		
				 num_donimated[j]++;
			 }
			 else if (temp == 1)//j donimate i
			 {
				 donimate[j].push_back(i);
				 num_donimated[i]++;
			 }
		 }
		 if (num_donimated[i] == 0)
		 {
			 population[i].rank = 1;
			 front[0].push_back(i);
		 }
	}

	count = 1;
	while (front[count-1].size() > 0)
	{
		for (i = 0; i < front[count-1].size(); i++)
		{
			for (j = 0; j < donimate[front[count-1][i]].size(); j++)
			{
				int temp_index = donimate[front[count-1][i]][j];
				num_donimated[temp_index]--;
				if (num_donimated[temp_index] == 0)
				{
					population[temp_index].rank = count + 1;
					front[count].push_back(temp_index);
				}
			}
		}
		count++;
	}
}

void CMOEA::crowd_distance_assignment(vector <CIndividual> &population, int k)
{
	int i,j,m;
	for (i = 0; i < front[k].size(); i++)	population[front[k][i]].crowd_distance = 0;
	vector<int>	temp_seq;
 	for (m = 0; m < nobj; m++)
	{
		double fmax_m,fmin_m;   		
		for (i = 0; i < front[k].size(); i++)
		{
			temp_seq.push_back(front[k][i]);
		}
		for (i = 0; i < temp_seq.size(); i++)
		{
			for (j = i + 1; j < temp_seq.size(); j++)
			{
				if (population[temp_seq[i]].y_obj[m] > population[temp_seq[j]].y_obj[m])
				{
					int temp_index = temp_seq[i];
					temp_seq[i] = temp_seq[j];
					temp_seq[j] = temp_index;
				}
			}
		}
		fmax_m = population[temp_seq[front[k].size() - 1]].y_obj[m];
		fmin_m = population[temp_seq[0]].y_obj[m];
		population[temp_seq[front[k].size() - 1]].crowd_distance += inf;
		population[temp_seq[0]].crowd_distance += inf;
		for (i = 1; i < front[k].size() - 1; i++)
			population[temp_seq[i]].crowd_distance += fabs(population[temp_seq[i+1]].y_obj[m] - population[temp_seq[i-1]].y_obj[m]) / fabs(fmax_m - fmin_m+1e-10);

		temp_seq.clear();
	} 
	temp_seq.clear();
}

void CMOEA::nondonimate_crowd_compare(vector <CIndividual> &population)
{
	int i,j,k;
	for (j = 0; j < population.size(); j++)		population[j].crowd_distance = 0;

	i = 0;
	population_evol.clear();
	while (population_evol.size() + front[i].size() < pops)
	{
		crowd_distance_assignment(population, i);
		for (j = 0; j < front[i].size(); j++)	population_evol.push_back(population[front[i][j]]);
		i++;
	}	

	crowd_distance_assignment(population, i);

	vector<int>	temp_seq;
	for (j = 0; j < front[i].size(); j++)		temp_seq.push_back(front[i][j]);

	if (population_evol.size() + front[i].size() > pops)
	{
		//get the first pops-population_evol.size();
		for (j = 0; j < pops-population_evol.size(); j++)
		{
			for (k = j + 1; k < front[i].size(); k++)
			{
				if (population[temp_seq[j]].crowd_distance < population[temp_seq[k]].crowd_distance)
				{
					int temp_index = temp_seq[j];
					temp_seq[j] = temp_seq[k];
					temp_seq[k]= temp_index;
				}
			}
		}
	}

	int num_temp = population_evol.size();
	for (j = 0; j < pops - num_temp; j++)
	{
		population_evol.push_back(population[temp_seq[j]]);
	}
	temp_seq.clear();
}

/*
void CMOEA::crossover_sbx(vector <CIndividual> &population_select)
{
	int i,j,index[2];
	for (i = 0; i < int(population_select.size()/2); i++)
	{
		for (j = 0; j < 2; j++)
		{
			rnd1 = int(population_evol.size() * rnd_uni(&rnd_uni_init));
			rnd2 = int(population_evol.size() * rnd_uni(&rnd_uni_init));
			while (rnd2 == rnd1) {rnd2 = int(population_evol.size() * rnd_uni(&rnd_uni_init));}
			if (population_evol[rnd1].rank < population_evol[rnd2].rank || 
				(population_evol[rnd1].rank == population_evol[rnd2].rank && 
				(population_evol[rnd1].crowd_rank < population_evol[rnd2].crowd_rank)))
				index[j] = rnd1;
			else
				index[j] = rnd2;
		}
		real_sbx_xoverA(population_select[index[0]], population_select[index[1]], population_new[2*i], population_new[2*i+1]);
	}
}
*/

void CMOEA::crossover_de()
{
	int i,j,index[3];
	int rnd1,rnd2;
	for (i = 0; i < pops; i++)
	{		
		for (j = 0; j < 3; j++)
		{
			rnd1 = int(population_evol.size() * rnd_uni(&rnd_uni_init));
			rnd2 = int(population_evol.size() * rnd_uni(&rnd_uni_init));
			while (rnd2 == rnd1) {rnd2 = int(population_evol.size() * rnd_uni(&rnd_uni_init));}
			if (population_evol[rnd1].rank < population_evol[rnd2].rank || 
				(population_evol[rnd1].rank == population_evol[rnd2].rank && 
				(population_evol[rnd1].crowd_distance < population_evol[rnd2].crowd_distance)))
				index[j] = rnd1;
			else
				index[j] = rnd2;
		}
		//2.2 de
		diff_evo_xoverB(population_evol[index[0]], population_evol[index[1]], population_evol[index[2]], population_new[i]);
	}
}

void CMOEA::crossover_sbx(vector <CIndividual> &population_select)
{
	for (int i = 0; i < int(population_select.size()/2); i++)	
		real_sbx_xoverA(population_select[2*i], population_select[2*i+1], population_new[2*i], population_new[2*i+1]);
}
/*
void CMOEA::crossover_de()
{
	int i,j;
	int rnd1,rnd2;
	for (i = 0; i < pops; i++)
	{
		int index[3];
		for (j = 0; j < 3; j++)
		{
			rnd1 = int(num_neighbour * rnd_uni(&rnd_uni_init));
			rnd2 = int(num_neighbour * rnd_uni(&rnd_uni_init));
			while (rnd2 == rnd1) {rnd2 = int(num_neighbour * rnd_uni(&rnd_uni_init));}
			if (population_evol[rnd1].rank < population_evol[rnd2].rank || 
			(population_evol[rnd1].rank == population_evol[rnd2].rank && 
			(population_evol[rnd1].crowd_distance < population_evol[rnd2].crowd_distance)))
				index[j] = rnd1;
			else
				index[j] = rnd2;
		}
		//2.2 de
		diff_evo_xoverB(population_evol[index[0]], population_evol[index[1]], population_evol[index[2]], population_new[i]);
	}
}
*/
void CMOEA::mutation_polynomial(vector <CIndividual> &population_select)
{
	for (int i = 1; i < population_select.size(); i++)	realmutation(population_select[i], 1.0/nvar);
}

void CMOEA::evulation(vector <CIndividual> &population_select)
{
	for (int i = 0; i < population_select.size(); i++)		population_select[i].obj_eval(); 
}

void  CMOEA::update_store()
{
	int i,j,k;
	for (i = 0; i < population_new.size(); i++)
	{
		int flag = 0;
		for (j = 0; j < population_store.size(); j++)
		{
			int ret = donimate_judge(population_store[j].y_obj, population_new[i].y_obj);			
			if ( ret == -1)
			{//population_store[j] donimate population_new[i]
				flag = -1;
				break;
			}
			else if (ret == 1)
			{//population_new[i] donimate population_store[j]				
			    population_store.erase(population_store.begin() + j);
				flag = 1;
				j--;
			}
		}
		if (flag == 1 || flag == 0)	population_store.push_back(population_new[i]);		
	}

	if (population_store.size() > num_store)	
	{
		vector <CIndividual> population_temp;
		int index = int(rnd_uni(&rnd_uni_init) * population_store.size());
		population_temp.push_back(population_store[index]);
		population_store.erase(population_store.begin() + index);
		
		while (population_temp.size() < pops)
		{
			double max_dist = 0;
			index = 0;
			for (i = 0; i < population_store.size(); i++)
			{
				//get the shortest distance form  population_store[i]  to population_evol
				double dist = 1e10;
				for (j = 0; j < population_temp.size(); j++)
				{
					double d = dist_vector(population_store[i].y_obj, population_temp[j].y_obj);
					if (d < dist)	dist = d;
				}
				if (dist > max_dist)	
				{//get 
					max_dist = dist;
					index = i;
				}
			}
			population_temp.push_back(population_store[index]);
			population_store.erase(population_store.begin() + index);
		}
		population_store.swap(vector<CIndividual>());
		for (i = 0; i < population_temp.size(); i++)	population_store.push_back(population_temp[i]);
		population_temp.clear();
	}
}

void CMOEA::evol_population()
{
	int i,j;
	//0.	combine parent and offset population	
	population_union.clear();
	for (i=0; i < population_evol.size(); i++)	population_union.push_back(population_evol[i]);
	for (i=0; i < population_new.size(); i++)	population_union.push_back(population_new[i]);	
	
	//1	fitness assign:	using fast non-donimate sort
	fitnessAssign_fast_nondonimate_sort(population_union);
	//3.1	selection:		selection the reproduced father indivivuals	
	nondonimate_crowd_compare(population_union);
	//3.2	crossover:		using SBX crossover
	if(!strcmp("SBX", crossoverType))	crossover_sbx(population_evol);
	else if (!strcmp("DE", crossoverType))	crossover_de();	
	//3.3	mutation:		using polynomial mutation
	mutation_polynomial(population_new);
	//3.4	evulation:		
	evulation(population_new);
	//update_store();
	nfes = 	nfes + population_new.size();
}


void CMOEA::exec_emo(int run)
{
	int i,j;
    char filename[1024];	
	seed = (seed + 23)%1377;					
	rnd_uni_init = -(long)seed;					//	initialization the random seed	
	nfes      = 0;								//	record the number of fuction evalution 
	init_population();							//	initialization the population
	//1	fitness assign:	using fast non-donimate sort
	population_union.clear();
	for (i=0; i < population_evol.size(); i++)	population_union.push_back(population_evol[i]);
	fitnessAssign_fast_nondonimate_sort(population_union);
	//2	selection:		selection the reproduced father indivivuals using crowd compare
 	nondonimate_crowd_compare(population_union);
	//3	crossover:		using SBX crossover
	if(!strcmp("SBX", crossoverType))	crossover_sbx(population_evol);
	else if (!strcmp("DE", crossoverType))	crossover_de();
	//4	mutation:		using polynomial mutation
	mutation_polynomial(population_new);
	//5	evulation:		
	evulation(population_new);
	nfes = 	nfes + population_new.size();

	int gen = 1;								//	record the times of population evolution
	int cur_gen = 0;
	//--- evolution ----
	while(nfes<=num_max_evulation)
	{		
		evol_population();						//	evolution population		
		gen++;
		
		if (nfes >= cur_gen*num_max_evulation/30)
		{
			sprintf(filename,"POF/%d/POF_NSGA2_%s_RUN%d.txt", cur_gen*10,strTestInstance,run);
			save_front(filename);
			sprintf(filename,"POS/%d/POS_NSGA2_%s_RUN%d.txt", cur_gen*10,strTestInstance,run);
			save_pos(filename);

			cur_gen++;
		}
	}
	sprintf(filename,"POF/POF_NSGA2_%s_RUN%d.txt", strTestInstance,run);
	save_front(filename);
	sprintf(filename,"POS/POS_NSGA2_%s_RUN%d.txt", strTestInstance,run);
	save_pos(filename);

	population_new.clear();
	population_evol.clear();
	population_store.clear();	
	population_union.clear();
	//clear;
	for (i = 0; i < front.size(); i++)		front[i].clear();
	front.clear();
	num_donimated.clear();
	for (i = 0; i < front.size(); i++)		donimate[i].clear();
	donimate.clear();
}

	
void CMOEA::save_population_union(vector <CIndividual> population_union, char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<population_union[n].y_obj[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void CMOEA::save_front(char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nobj;k++)
			fout<<population_evol[n].y_obj[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void CMOEA::save_front_select(char saveFilename[1024])
{
	int i,j;
	if (nobj == 2)
	{
		for (i = 0; i < population_store.size(); i++)
		{
			for (j = i + 1; j <population_store.size(); j++)
			{
				if (fabs(population_store[i].y_obj[0] - population_store[j].y_obj[0]) < emxitong)
				{
					if (population_store[i].y_obj[1] < population_store[j].y_obj[1])
					{
						population_store.erase(population_store.begin() + j);
						j--;
					}
					else
					{
						population_store.erase(population_store.begin() + i);
						i--;
						break;
					}
				}
				else if (fabs(population_store[i].y_obj[1] - population_store[j].y_obj[1]) < emxitong)
				{
					if (population_store[i].y_obj[0] < population_store[j].y_obj[0])
					{
						population_store.erase(population_store.begin() + j);
						j--;
					}
					else
					{
						population_store.erase(population_store.begin() + i);
						i--;
						break;
					}
				}
			}
		}
	}

    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	if (population_store.size() < pops)
	{
		for(int n=0; n<population_store.size(); n++)
		{
			for(int k=0;k<nobj;k++)
				fout<<population_store[n].y_obj[k]<<"  ";
			fout<<"\n";
		}
	}
	else
	{
		int i,j;
		vector<CIndividual> population_output;
		int index = int(rnd_uni(&rnd_uni_init) * population_store.size());
		population_output.push_back(population_store[index]);
		population_store.erase(population_store.begin() + index);
		
		while (population_output.size() < pops)
		{
			double max_dist = 0;
			index = 0;
			for (i = 0; i < population_store.size(); i++)
			{
				double dist = 1e10;
				for (j = 0; j < population_output.size(); j++)
				{
					double d = dist_vector(population_store[i].y_obj, population_output[j].y_obj);
					if (d < dist)	dist = d;
				}
				if (dist > max_dist)	
				{
					max_dist = dist;
					index = i;
				}
			}
			population_output.push_back(population_store[index]);
			population_store.erase(population_store.begin() + index);
		}
		for(int n=0; n<population_output.size(); n++)
		{
			for(int k=0;k<nobj;k++)
				fout<<population_output[n].y_obj[k]<<"  ";
			fout<<"\n";
		}
		population_output.clear();
		//population_output.swap(vector<CIndividual>());
	}
	
	fout.close();
}

void CMOEA::save_pos(char saveFilename[1024])
{
    std::fstream fout;
	fout.open(saveFilename,std::ios::out);
	for(int n=0; n<pops; n++)
	{
		for(int k=0;k<nvar;k++)
			fout<<population_evol[n].x_var[k]<<"  ";
		fout<<"\n";
	}
	fout.close();
}

void CMOEA::save_individual(char saveFilename_x[1024],char saveFilename_y[1024],vector<double> x,vector<double> y)
{
	std::fstream fout;
	fout.open(saveFilename_x,std::ios::out);
	int i,n;
	for(n=0; n<nvar; n++)
	{
		fout<<x[n]<<"  ";
	}
	fout.close();

	//std::fstream fout;
	fout.open(saveFilename_y,std::ios::out);
	for(n=0; n<nobj; n++)
	{
		fout<<y[n]<<"  ";
	}
	fout.close();
}

bool CMOEA::boundy_checkout(CIndividual child)
{
	int i;
	for (i=0; i<nvar; i++)
	{
		if (child.x_var[i]<lowBound[i]||child.x_var[i]>uppBound[i]) return true;
	}
	for (i=0; i<nobj; i++)
	{
		if (child.y_obj[i]<0) return true;
	}
	return false;
}
#endif
