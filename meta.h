#include <iostream>
#include <vector>
#include <string>
#include <omp.h>
#include "formula.h"

using namespace std;

template<int N>
class formula_generator{

public:
	static inline void  EXEC(result *p_result, int elcount, int minimum[], int maximum[], double masses[], int *current, double pre_mass, double loMass, double hiMass, int k)
	{
		int c = min(int((hiMass-pre_mass)/masses[N-1]),maximum[N-1]);
		if (N == elcount)
		{
		#pragma omp parallel for 
			for (k = minimum[N-1]; k<=c; k++)
			{
				current[N-1] = k;
				double current_mass_i= pre_mass+masses[N-1]*k;
				int *current_i = new int [10];
				for (int i=0;i<elcount;i++)
				{
					current_i[i]=current[i];
				}

				formula_generator<N-1>::EXEC(p_result, elcount, minimum, maximum, masses, current_i, current_mass_i, loMass, hiMass, k);
				delete []current_i;
			}

		}
		else
		{
			for (int i = minimum[N-1]; i<=c; i++)
			{

				current[N-1] = i;

				double current_mass_i= pre_mass+masses[N-1]*i;

				formula_generator<N-1>::EXEC(p_result, elcount, minimum, maximum, masses, current, current_mass_i, loMass, hiMass, k);

			}
		}
	}
};


template<>
class formula_generator<0>{

public:
	static inline void  EXEC(result *p_result, int elcount, int minimum[], int maximum[], double masses[], int *current, double pre_mass, double loMass, double hiMass, int k)
	{

		if (pre_mass >= loMass && pre_mass <= hiMass && p_result->len <= 1000000/*1000000*/)
		{
			#pragma omp critical (p_result)
			{
				for ( int i = 0; i < elcount; ++i)
				{
					if (i != elcount-1)
					{
						p_result->data[p_result->len*elcount+i] = current[i];
					}
					else
					{
						p_result->data[p_result->len*elcount+i] = k;
					}
				}
				p_result->mass[p_result->len] = pre_mass;
				p_result->len++;
			}
		}
	}
};