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

double delta(double measuredMass, double countedMass, string unit="ppm"){
	// Calculate error between measured Mass and counted Mass in specified units.
	if(strcmp(unit.c_str(),"ppm")==0)
		return (measuredMass - countedMass) / countedMass*1000000;
	else if(strcmp(unit.c_str(),"Da")==0)
		return (measuredMass - countedMass);
	else if(strcmp(unit.c_str(),"%")==0)
		return (measuredMass - countedMass) / countedMass*100;
	else printf("Unknown units for delta!");
	    return false;
}



void PFG(ResultsWriter & writer, int elcount, int minimum[], int maximum[], double masses[], int *current, double pre_mass,double loMass, double hiMass)
	//call the template metaprogramming
{
	int k = 0;
	switch (elcount)
	{
		case 1:  formula_generator<1>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 2:  formula_generator<2>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 3:  formula_generator<3>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 4:  formula_generator<4>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 5:  formula_generator<5>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 6:  formula_generator<6>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 7:  formula_generator<7>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 8:  formula_generator<8>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 9:  formula_generator<9>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 10: formula_generator<10>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 11: formula_generator<11>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 12: formula_generator<12>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 13: formula_generator<13>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 14: formula_generator<14>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
		case 15: formula_generator<15>::EXEC(writer, elcount, minimum, maximum, masses, current, pre_mass, loMass, hiMass, k); break;
	default:
		break;
	}
}

void calculation(double currentmass, vector<char*> compositions, int mincount[], int maxcount[], vector<char*> rules, float tolerance, int charge, char *agentformula, int agentcharge, string unit, string outfile)
{
	double _mz = mz(currentmass, 0, charge, agentformula, agentcharge);
	int elcount = compositions.size();
	double finalmass=_mz;
	double loMass,hiMass;
	if (strcmp(unit.c_str(),"ppm") == 0 && charge == 0){
		loMass = finalmass - (finalmass/1000000) * tolerance;
		hiMass = finalmass + (finalmass/1000000) * tolerance;
	}
	else if (charge != 0 && strcmp(unit.c_str(),"ppm") == 0)
	{
		loMass = finalmass - abs(charge) * (finalmass/1000000) * tolerance;
		hiMass = finalmass + abs(charge) * (finalmass/1000000) * tolerance;
	}
	else if (charge!=0 && strcmp(unit.c_str(),"Da")==0)
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
		/*
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
		}*/
		maxcount[i] = min(maxcount[i], int(hiMass / _elmass[i]));
	}

	int *current = new int [elcount];
	for (int i = 0; i < elcount; i++)
	{
		current[i] = mincount[i];
	}
	double pre_mass = 0.0;


	ResultsWriter writer(outfile, compositions, rules, charge, agentformula, agentcharge, currentmass, unit);
	PFG(writer, elcount, mincount, maxcount, _elmass, current, pre_mass, loMass, hiMass);
	writer.writeResults();

	delete []_elmass;
	delete []current;
}

ResultsWriter::ResultsWriter(string outfile, vector<char*> compositions, vector<char*> rules, int charge, char * agentformula, int agentcharge, double currentmass, string unit)
{
	m_outfile = outfile;
	m_compositions = compositions;
	m_rules = rules;
	m_charge = charge;
	m_agentformula = agentformula;
	m_agentcharge = agentcharge;
	m_currentmass = currentmass;
	m_unit = unit;

	m_elcount = compositions.size();

	m_poolSize = 1000000;

	p_result = new result[m_poolSize];
	p_result->data = new int[m_poolSize * m_elcount];
	p_result->len = 0;
	p_result->mass = new double[m_poolSize];

}

ResultsWriter::~ResultsWriter()
{
	delete[]p_result->mass;
	delete[]p_result->data;
	delete[]p_result;
}


void ResultsWriter::writeResults()
{
	FILE *fpt;

	if (m_outfile.size())
	{
		fpt = fopen(m_outfile.c_str(), "w");
		if (fpt == NULL)
		{
			fpt = stdout;
		}
	}
	else
	{
		fpt = stdout;
	}

	fprintf(fpt, "formula\t mass\t mz\t error\t rdbe\n");
	char *temp = new char[100];
	char *s = new char[100];

	//start = clock();
	for (int i = 0; i < p_result->len; i++)
	{
		temp[0] = '\0';
		s[0] = '\0';
		int countC = 0;
		int countH = 0;
		int countO = 0;
		int countN = 0;
		int countP = 0;
		int countS = 0;
		double rdbevalue;
		vector<int> count;
		for (int j = 0; j < m_elcount; j++)
		{
			if (strcmp(m_compositions[j], "C") == 0)
			{
				countC = p_result->data[i*m_elcount + j];
			}
			else if (strcmp(m_compositions[j], "H") == 0)
			{
				countH = p_result->data[i*m_elcount + j];
			}
			else if (strcmp(m_compositions[j], "O") == 0)
			{
				countO = p_result->data[i*m_elcount + j];
			}
			else if (strcmp(m_compositions[j], "N") == 0)
			{
				countN = p_result->data[i*m_elcount + j];
			}
			else if (strcmp(m_compositions[j], "P") == 0)
			{
				countP = p_result->data[i*m_elcount + j];
			}
			else if (strcmp(m_compositions[j], "S") == 0)
			{
				countS = p_result->data[i*m_elcount + j];
			}

			const unsigned int elnum = p_result->data[i*m_elcount + j];
			if (elnum != 0)
			{
				strncat(temp, m_compositions[j], strlen(m_compositions[j]));
				if (elnum != 1)
				{
					sprintf(s, "%d", p_result->data[i*m_elcount + j]);
					strncat(temp, s, strlen(s));
				}
			}

			count.push_back(p_result->data[i*m_elcount + j]);
		}
		rdbevalue = rdbe(m_compositions, count);

		bool t = frules(temp, countC, countH, countO, countN, countP, countS, rdbevalue, m_rules);
		if (t == 1)
		{
			double mass = p_result->mass[i];
			double _mass = mass;
			double fmz = mz(_mass, m_charge, 0, m_agentformula, m_agentcharge);
			double error = delta(m_currentmass, fmz, m_unit);
			fprintf(fpt, "%s\t %f\t %f\t %f\t %f\n", temp, mass, fmz, error, rdbevalue);
		}
	}
	fclose(fpt);
	delete[]temp;
	delete[]s;
}