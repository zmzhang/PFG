#include <iostream>
#include <string.h>
#include <malloc.h>
#include <ctime>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "formula.h"
#include "meta.h"

#if defined __GNUC__
#pragma GCC diagnostic ignored "-Wwrite-strings"
#endif

using namespace std;

double *elmass(vector<char*> compositions, int size){
	double *_elmass;
	_elmass = new double [size];
	double temp;
	for(int i = 0; i < size; i++){
		temp = elementmass(compositions[i]);
		_elmass[i] = temp;
	}
	return _elmass;
}

double delta(double measuredMass, double countedMass, char units[]="ppm"){
	// Calculate error between measured Mass and counted Mass in specified units.
	if(strcmp(units,"ppm")==0)
		return (measuredMass - countedMass) / countedMass*1000000;
	else if(strcmp(units,"Da")==0)
		return (measuredMass - countedMass);
	else if(strcmp(units,"%")==0)
		return (measuredMass - countedMass) / countedMass*100;
	else printf("Unknown units for delta!");
	    return false;
}



void PFG(result *p_result, int elcount, int minimum[], int maximum[], double masses[], int *current, double pre_mass,double loMass, double hiMass)
	//call the template metaprogramming
{
	int k = 0;
	switch (elcount)
	{
		case 1:  formula_generator<1>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 2:  formula_generator<2>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 3:  formula_generator<3>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 4:  formula_generator<4>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 5:  formula_generator<5>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 6:  formula_generator<6>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 7:  formula_generator<7>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 8:  formula_generator<8>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 9:  formula_generator<9>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 10: formula_generator<10>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 11: formula_generator<11>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 12: formula_generator<12>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 13: formula_generator<13>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 14: formula_generator<14>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 15: formula_generator<15>::EXEC(p_result, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
	default:
		break;
	}
}

void calculation(double currentmass, vector<char*> compositions, int mincount[], int maxcount[], vector<char*> rules, float tolerance, int charge, char *agentformula, int agentcharge, char units[])
{
	double _mz = mz(currentmass, 0, charge, agentformula, agentcharge);
	int elcount = compositions.size();
	double finalmass=_mz;
	double loMass,hiMass;
	if (strcmp(units,"ppm") == 0 && charge == 0){
		loMass = finalmass - (finalmass/1000000) * tolerance;
		hiMass = finalmass + (finalmass/1000000) * tolerance;
	}
	else if (charge != 0 && strcmp(units,"ppm") == 0)
	{
		loMass = finalmass - abs(charge) * (finalmass/1000000) * tolerance;
		hiMass = finalmass + abs(charge) * (finalmass/1000000) * tolerance;
	}
	else if (charge!=0 && strcmp(units,"Da")==0)
	{
		loMass = finalmass - abs(charge) * tolerance;
		hiMass = finalmass + abs(charge) * tolerance;
	}
	else 
	{
		loMass = finalmass - tolerance;
		hiMass = finalmass + tolerance;
	}

	double *_elmass = elmass(compositions, elcount);

	for(int i = 0; i < elcount; i++)
	{
		if (strcmp(compositions[i], "C") == 0)
		{
			if (_mz < 500)
			{
				maxcount[i] = min(maxcount[i], 40);
			}
			else if (_mz < 1000)
			{
				maxcount[i] = min(maxcount[i], 80);
			}
			else if (_mz < 2000)
			{
				maxcount[i] = min(maxcount[i], 160);
			}
			else maxcount[i] = min(maxcount[i], 180);
		}
		else if (strcmp(compositions[i], "H") == 0)
		{
			if (_mz < 500)
			{
				maxcount[i] = min(maxcount[i], 80);
			}
			else if (_mz < 1000)
			{
				maxcount[i] = min(maxcount[i], 130);
			}
			else if (_mz < 2000)
			{
				maxcount[i] = min(maxcount[i], 250);
			}
			else maxcount[i] = min(maxcount[i], 300);
		}
		else if (strcmp(compositions[i], "O") == 0)
		{
			if (_mz < 500)
			{
				maxcount[i] = min(maxcount[i], 20);
			}
			else if (_mz < 1000)
			{
				maxcount[i] = min(maxcount[i], 30);
			}
			else if (_mz < 2000)
			{
				maxcount[i] = min(maxcount[i], 70);
			}
			else maxcount[i] = min(maxcount[i], 90);
		}
		else if (strcmp(compositions[i], "N") == 0)
		{
			if (_mz < 500)
			{
				maxcount[i] = min(maxcount[i], 20);
			}
			else if (_mz < 1000)
			{
				maxcount[i] = min(maxcount[i], 30);
			}
			else if (_mz < 2000)
			{
				maxcount[i] = min(maxcount[i], 40);
			}
			else maxcount[i] = min(maxcount[i], 60);
		}
		maxcount[i] = min(maxcount[i], int(hiMass / _elmass[i]));

	}

	result *p_result = new result [1000000];
	p_result->data = new int [1000000 * elcount];
	p_result->len = 0;
	p_result->mass = new double [1000000];
	int *current = new int [elcount];

	for (int i = 0; i < elcount; i++)
	{
		current[i] = mincount[i];
	}
	double pre_mass = 0.0;
	PFG(p_result, elcount, mincount, maxcount, _elmass, current, pre_mass, loMass, hiMass);

	FILE *fpt; 
	fpt=fopen("result.txt","w");
	fprintf(fpt,"formula\t mass\t mz\t error\t rdbe\n");
	char *temp = new char[100];
	char *s = new char[100];

	//start = clock();
	for(int i = 0; i < p_result->len; i++)
	{
		temp[0]='\0';
		s[0]='\0';
		int countC = 0;
		int countH = 0;
		int countO = 0;
		int countN = 0;
		int countP = 0;
		int countS = 0;
		double rdbevalue;
		vector<int> count;
		for(int j = 0; j < elcount; j++)
		{
			if (strcmp(compositions[j], "C") == 0)
			{
				countC = p_result->data[i*elcount + j];
			}
			else if (strcmp(compositions[j], "H") == 0)
			{
				countH = p_result->data[i*elcount + j];
			}
			else if (strcmp(compositions[j], "O") == 0)
			{
				countO = p_result->data[i*elcount + j];
			}
			else if (strcmp(compositions[j], "N") == 0)
			{
				countN = p_result->data[i*elcount + j];
			}
			else if (strcmp(compositions[j], "P") == 0)
			{
				countP = p_result->data[i*elcount + j];
			}
			else if (strcmp(compositions[j], "S") == 0)
			{
				countS = p_result->data[i*elcount + j];
			}			

			if(p_result->data[i*elcount + j] !=0)
			{
				sprintf(s, "%s%d",compositions[j], p_result->data[i*elcount + j]);
				strncat(temp, s, strlen(s));
			}

			count.push_back(p_result->data[i*elcount + j]);
		}
		rdbevalue = rdbe(compositions, count);
	
		int t=frules(temp, countC, countH, countO, countN, countP, countS, rdbevalue, rules);
		if(t == 1)
		{
			double mass = p_result->mass[i];
			double _mass = mass;
			double fmz = mz(_mass, charge, 0, agentformula, agentcharge);
			double error = delta(currentmass, fmz, units);
			fprintf(fpt, "%s\t %f\t %f\t %f\t %f\n", temp, mass, fmz, error, rdbevalue);
		}
	}
	fclose(fpt);
	delete []temp;
	delete []s;
	delete []_elmass;
	delete []p_result->mass;
	delete []p_result->data;
	delete []p_result;
	delete []current;
}