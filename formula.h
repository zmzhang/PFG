#ifndef __FORMULA_H   
#define __FORMULA_H  
#define SUBSLEN 10             
#define EBUFLEN 128         
#define BUFLEN 1024  
#define ELECTRON_MASS  0.00054857990924
#define min(a,b) a<b?a:b

#include <string>
#include <map>
#include <iostream>
#include <vector>
using namespace std;

typedef struct{
	char *name;
	char *symbol;
	int atomicNumber;
	int isotopecount;
	int massnumber[10];
	double isotope[10];
	double abundance[10];
	int valence;
}ellist;

typedef struct{
	char **comp;
	int *atomcount;
	int elcount;
}composition;

typedef struct{
	int *data;
	int len;
	double *mass;
}result;

typedef struct{
	char *_formula;
	double *_mass;
	double *_mz;
	double error;
	double errorDa;
	double rdbe;
}output;

//elements
double elementmass(char *string);
double rdbe(vector<char*> compositions, vector<int> count);
double formulaMass(vector<char*> compositions, vector<int> count);

//compound
double mz(double currentmass,int finalcharge, int currentcharge, char *agentformula, int agentcharge);
static double GLHC[2]={0.2, 3};
static double GLNOPSC[4]={2, 1.2, 0.32, 0.65};
static double GLRDBE[2]={0, 40};
bool frules(char *string, int countC, int countH, int countO, int countN, int countP, int countS, double rdbevalue, 
		   vector<char*> rules);

//formula
double *elmass(vector<char*> compositions, int size);
double delta(double measuredMass, double countedMass, char units[]);
void PFG(result *p_result, int elcount, int minimum[], int maximum[], double masses[], int *current, double pre_mass,double loMass, double hiMass);
void calculation(double currentmass, vector<char*> compositions, int mincount[], int maxcount[], vector<char*> rules, float tolerance=1,  
		  int charge = 0, char *agentformula = "H", int agentcharge=1, string unit="ppm", string outfile="");
#endif