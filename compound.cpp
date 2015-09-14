#include <iostream>
#include <string>
#include <map>
#include <string.h>
#include <stdlib.h>
#include <ctime>
#include "formula.h"
#include "temp.h"
#include <sstream>
#include <fstream>
#include <iomanip>
#include "getoptpp/getopt_pp.h"
#include "getoptpp/getopt_pp_standalone.h"

using namespace GetOpt;
#pragma GCC diagnostic ignored "-Wwrite-strings"

double mz(double currentmass,int finalcharge, int currentcharge, char *agentformula="H",int agentcharge=1)
	//calculate the m/z value for given mass and charge
{
	double agentMass;
	double massMo = 0.0;
	double massAv = 0.0;

	if(strcmp(agentformula,"e") == 0)
	{
		agentMass = ELECTRON_MASS;
	}
	else{
		agentMass = elementmass(agentformula);
		agentMass = agentMass - agentcharge*ELECTRON_MASS;
	}
	int agentCount = currentcharge/agentcharge;
	if (currentcharge != 0)
	{
		massMo = currentmass * abs(currentcharge) - agentMass * agentCount;
		currentmass = massMo;
	}
	else
	{
		massMo = currentmass;
	}
	if (finalcharge == 0)
	{
		return currentmass;
	}
	agentCount = finalcharge / agentcharge;

	currentmass = (currentmass + agentMass * agentCount) / abs(finalcharge);
	return currentmass;
}

int frules(char *string, int countC, int countH, int countO, int countN, int countP, int countS, double rdbevalue, vector<char*> rules)
	//Check formula rules for a given compound.
{
	int size = rules.size();

	if(countC){
		double ratioHC = countH / 1.000 / countC;
		double ratioNC = countN / 1.000 / countC;
		double ratioOC = countO / 1.000 / countC;
		double ratioPC = countP / 1.000 / countC;
		double ratioSC = countS / 1.000 / countC;

		for (int i = 0; i < size; i++)
		{ 
			if(strcmp( rules[i], "HC" ) == 0 ){
				if ((ratioHC < GLHC[0] ) || ( ratioHC > GLHC[1])){
					return false;
				}
			}

			if (strcmp( rules[i], "NOPSC" ) == 0 ){
				if (ratioNC > GLNOPSC[0] || ratioOC > GLNOPSC[1] || ratioPC > GLNOPSC[2] || ratioSC > GLNOPSC[3]){
					return false;
				}
			}

			if ((strcmp(rules[i], "NOPS") == 0) && (countN > 1 && countO > 1 && countP > 1 && countS > 1)){
				if((countN >= 10) || (countO >= 20) || (countP >= 4) || (countS >= 3)){
					return false;
				}
			}

			if ((strcmp(rules[i], "NOPS") == 0) && (countN > 3 && countO > 3 && countP > 3)){
				if ((countN >= 11) || (countO >= 22) || (countS >= 6)){
					return false;
				}
			}

			if ((strcmp(rules[i], "NOPS") == 0) && (countO > 1 && countP > 1 && countS > 1)){
				if ((countO >= 14) || (countP >= 3) || (countS >= 3)){
					return false;
				}
			}

			if ((strcmp(rules[i], "NOPS") == 0) && (countP > 1 && countS > 1 && countN > 1)){
				if ((countP >= 3) || (countS >= 3) || (countN >= 4)){
					return false;
				}
			}

			if ((strcmp(rules[i], "NOPS") == 0) && ( countN > 6 && countO > 6 && countS > 6 )){
				if ((countN >= 19) || (countO >= 14) || (countS >= 8)){
					return false;
				}
			}
			if (strcmp(rules[i], "RDBE") ==0 ){
				if (rdbevalue < GLRDBE[0] || rdbevalue > GLRDBE[1]){
					return false;
				}
			}
			if (strcmp(rules[i], "lewis") ==0 ){
				if (rdbevalue - floor(rdbevalue) != 0){
					return false;
				}
			}
		}
	}
	return true;
}

vector<string> split(string str, char delimiter) {
	vector<string> internal;
	stringstream ss(str); // Turn the string into a stream.
	string tok;

	while(getline(ss, tok, delimiter)) {
		internal.push_back(tok);
	}

	return internal;
}

vector<char*> split_char(string str, char delimiter) {
	vector<char*> internal;
	stringstream ss(str); // Turn the string into a stream.
	string tok;

	while(getline(ss, tok, delimiter)) {
		stringstream s1;
		char* x = new char[100];
		s1<<tok;
		s1>>x;
		internal.push_back(x);
		delete []x;
	}

	return internal;
}
void help()
{
	cout<<"PFG is a program to calculate the possible elemental compositions for a given mass." << endl
		<<"Authors: Mingjing Zhang & Zhimin Zhang" << endl
		<<"usage: PFG valid command line options are: " << endl
		<<"-h or --help       The help screen." <<endl
		<<"-m mass            Set mass." <<endl
		<<"-t tol             Set tolerance to tol 'ppm' ( default 5 )." << endl
		<<"-c charge          Set charge to be calculated." << endl
		<<"-r rules           Set rules to constrain formulas"<<endl
		<<"--X a-b            Set atom range a (min) to b (max) of element X. "<<endl
		<<"                   some of the valid elements:" << endl
		<<"          X           key      mass(6 decimals shown)" << endl
		<<"     ----------------------------------------------------" << endl
		<<"          C           --C           12.000000" << endl
		<<"          13C         --13C         13.003355"<<endl
		<<"          H           --H           1.007825" << endl
		<<"          D           --D           2.014102" <<endl
		<<"          N           --N           14.003074" << endl
		<<"          15N         --15N         15.000109"<<endl
		<<"          O           --O           15.994915" << endl
		<<"          P           --P           30.973762" << endl
		<<"          S           --S           31.972071" << endl
		<<"          F           --F           18.998403" << endl
		<<"          Cl          --Cl          34.968853" << endl
		<<"          Br          --Br          78.918338" << endl
		<<"          I           --I           126.904468" << endl
		<<"          Si          --Si          27.976927" << endl
		<<"          Na          --Na          22.989770" << endl
		<<"          K           --K           38.963707" << endl
		<<"--agentformula af  Set agent formula." << endl
		<<"--agentcharge  ac  Set the charge of the agent formula." << endl
		<< "for example:"<<endl
		<<"PFG  'Mesuximide' -m 203.0946 -t 5 --C 0-20 --H 0-40 --O 0-5 --N 0-5"<<endl
		<<"the result were stored in the result.txt"<<endl;
}

int main(int argc, char* argv[])
{
	stringstream ss;
	string _help;
	string tol;
	string rule;
	string _charge;
	string af;
	string ac;
	string _cmass;
	string C, H, N, O, P, S, Na, Si, F, Cl, Br, I, K, CX, D, NX;

	vector<char*> compositions;
	vector<int> mini;
	vector<int> maxi;

	GetOpt_pp ops(argc, argv);
	ops >> Option('h', "help", _help, "help");
	for (GetOpt_pp::short_iterator it = ops.begin(); it != ops.end(); ++it)
	{
		char short_opt[2] = {0,0};
		short_opt[0] = *it;
		if (short_opt[0] == 'h')
		{
			
			help();
		}
		
		//short_opt[0] = *it;
		//it >> vec;

		//std::cout << "\t" << short_opt << " has " << vec.size() << " integer arguments." << std::endl;
		//vec.clear();
	}

	ops >> Option('m', "mass", _cmass, "500");
	double cmass;
	ss.clear();
	ss << _cmass;
	ss >> cmass;

	ops >> Option('t', "tolerance", tol, "5");
	int tolerance;
	ss.clear();
	ss << tol;
	ss >> tolerance;

	ops >> Option('r', "rules", rule/*, "HC,NOPSC,NOPS,RDBE,lewis"*/);
	vector<char*> rules= split_char(rule, ',');

	ops >> Option('c', "charge", _charge, "0");
	ss.clear();
	int charge;
	ss << _charge;
	ss >> charge;

	ops >> Option("agentformula", af, "H");
	char *agentformula = new char [10];
	ss.clear();
	ss << af;
	ss >> agentformula;

	ops >> Option("agentcharge", ac, "1");
	int agentcharge;
	ss.clear();
	ss << ac;
	ss >> agentcharge;

	for (GetOpt_pp::long_iterator it = ops.begin(); it != ops.end(); ++it)
	{
		if (*it == "help")
		{
			help();
		}
		
		if (*it == "C")
		{
			ops >> Option("C", C, "0-40");
			compositions.push_back("C");
			vector<string> sep = split(C, '-');
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);	
		}
		if (*it == "H")
		{
			ops >> Option("H", H, "0-160");
			compositions.push_back("H");
			vector<string> sep = split(H, '-');
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "O")
		{
			ops >> Option("O", O, "0-30");
			vector<string> sep = split(O, '-');
			compositions.push_back("O");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "N")
		{
			ops >> Option("N", N, "0-20");
			vector<string> sep = split(N, '-');
			compositions.push_back("N");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "P")
		{
			ops >> Option("P", P, "0-0");
			vector<string> sep = split(P, '-');
			compositions.push_back("P");
			
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "S")
		{
			ops >> Option("S", S, "0-0");
			vector<string> sep = split(S, '-');
			compositions.push_back("S");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "F")
		{
			ops >> Option("F", F, "0-0");
			vector<string> sep = split(F, '-');
			compositions.push_back("F");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "Cl")
		{
			ops >> Option("Cl", Cl, "0-0");
			vector<string> sep = split(Cl, '-');
			compositions.push_back("Cl");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "Br")
		{
			ops >> Option("Br", Br, "0-0");
			vector<string> sep = split(Br, '-');
			compositions.push_back("Br");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "I")
		{
			ops >> Option("I", I, "0-0");
			vector<string> sep = split(I, '-');
			compositions.push_back("I");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "Si")
		{
			ops >> Option("Si", Si, "0-0");
			vector<string> sep = split(Si, '-');
			compositions.push_back("Si");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "Na")
		{
			ops >> Option("Na", Na, "0-0");
			vector<string> sep = split(Na, '-');
			compositions.push_back("Na");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "K")
		{
			ops >> Option("K", K, "0-0");
			vector<string> sep = split(K, '-');
			compositions.push_back("K");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "13C")
		{
			ops >> Option("13C", CX, "0-0");
			vector<string> sep = split(CX, '-');
			compositions.push_back("(13)C");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "D")
		{
			ops >> Option("D", D, "0-0");
			vector<string> sep = split(D, '-');
			compositions.push_back("D");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		if (*it == "15N")
		{
			ops >> Option("15N", NX, "0-0");
			vector<string> sep = split(NX, '-');
			compositions.push_back("(15)N");
			int tmp;
			ss.clear();
			ss << sep[0];
			ss >> tmp;
			mini.push_back(tmp);
			ss.clear();
			ss << sep[1];
			ss >> tmp;
			maxi.push_back(tmp);
		}
		/*std::cout << "\t" << *it << " has ";*/
	}

	int mincount[30];
	int maxcount[30];
	for (int i = 0; i < mini.size(); i++)
	{
		mincount[i] = mini[i];
	}
	for (int i = 0; i < maxi.size(); i++)
	{
		maxcount[i] = maxi[i];
	}
	calculation(cmass, compositions, mincount, maxcount, rules, tolerance, charge, agentformula, agentcharge);

	delete []agentformula;

	//ops >> Option('af', "agentformula", af, "H");
	//char *agentformula = new char [10];
	//ss.clear();
	//ss << af;
	//ss >> agentformula;

	//ops >> Option('ac', "agentcharge", ac, "1");
	//int agentcharge;
	//ss.clear();
	//ss << ac;
	//ss >> agentcharge;

	//vector<char*> compositions;
	//vector<int> mini;
	//vector<int> maxi;

	//ops >> Option('C', C, "0-40");
	//vector<string> sep1 = split(C, '-');
	//int a;
	//ss.clear();
	//ss << sep1[0];
	//ss >> a;
	//mini.push_back(a);
	//ss.clear();
	//ss << sep1[1];
	//ss >> a;
	//maxi.push_back(a);	

	//ops >> Option('H', H, "0-160");
	//vector<string> sep2 = split(H, '-');
	//ss.clear();
	//ss << sep2[0];
	//ss >> a;
	//mini.push_back(a);
	//ss.clear();
	//ss << sep2[1];
	//ss >> a;
	//maxi.push_back(a);

	//ops >> Option('N', N, "0-20");
	//vector<string> sep3 = split(N, '-');
	//ss.clear();
	//ss << sep3[0];
	//ss >> a;
	//mini.push_back(a);
	//ss.clear();
	//ss << sep3[1];
	//ss >> a;
	//maxi.push_back(a);

	//ops >> Option('O', O, "0-30");
	//vector<string> sep4 = split(O, '-');
	//ss.clear();
	//ss << sep4[0];
	//ss >> a;
	//mini.push_back(a);
	//ss.clear();
	//ss << sep4[1];
	//ss >> a;
	//maxi.push_back(a);




	//double cmass[1156];
	//ifstream file ( "masses.csv" );
	//string value, line; 
	//int i=0;
	//while (getline(file, line))  

	//{  
	//	istringstream is(line);
	//	while(getline(is, value, ','))
	//	{
	//		stringstream ss;
	//		ss<<value;
	//		ss>>cmass[i];
	//		i++;
	//	}
	//}
	//vector<char*> compositions;
	//vector<char*> rules;
	//compositions.push_back("C");
	//compositions.push_back("H");
	//compositions.push_back("O");
	//compositions.push_back("N");
 //	compositions.push_back("P");
 //	compositions.push_back("S");
	//rules.push_back("HC");
	//rules.push_back("NOPS");
 //	rules.push_back("NOPSC");
	//rules.push_back("RDBE");
	//rules.push_back("lewis");

	//FILE *fpts; 
	//fopen_s(&fpts,"caltime_P5S5_withoutparallel.xls","w");
	//for(int k = 0; k < 1155; k++){
	//	cout<<cmass[k]<<endl;
	//	fprintf (fpts, "%f\t" ,cmass[k]);
	//	double ttime = 0;
	//	double avtime = 0;
	//	for (int j = 0; j < 30; j++)
	//	{
	//		int mincount[]={0, 0, 0, 0, 0, 0};
	//		int maxcount[]={180, 300, 90, 60, 5, 5};
	//		double  start, finish;
	//		start = clock();
	//		test(cmass[k], compositions, mincount, maxcount, rules,GLHC,GLNOPSC,GLRDBE,3,"ppm",0,"H",1);
	//		finish = clock();
	//		double time = (finish - start) / CLOCKS_PER_SEC;
	//		fprintf (fpts, "%f\t", time);
	//		ttime = ttime + time;
	//	}
	//	avtime = ttime / 30;
	//	cout<<avtime<<endl;
	//	fprintf(fpts, "%s\t %f\n", "Average Time", avtime);
	//}
	//fclose(fpts);


	return 0;
}