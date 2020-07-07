/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <cmath>
#include "random.h"
#include "func.h"

using namespace std;
 
int main (int argc, char *argv[]) {

// Starting parameters
	vector<double> pos_iniz = {1., 1., 1.};
	vector<double> pos_in; 
	int thr = 1E6;
	ifstream p100, p210;
	int blocks = 100;
	int init = 1000;

// sampling..
	pos_in = pos_iniz; 
	cout << "psi_100" << endl;
	Metropolis(*psi_100, pos_in, thr, init, "Unif", 1.5, "Risultati/psi_100.uniform");
	cout << endl;
	pos_in = pos_iniz;
	cout << "psi_210" << endl;
	Metropolis(*psi_210, pos_in, thr, init, "Unif", 3.3, "Risultati/psi_210.uniform");
	
	cout << endl;
	
	pos_in = pos_iniz;
cout << "psi_100" << endl;
	Metropolis(*psi_100, pos_in, thr, init, "Gauss", 0.8, "Risultati/psi_100.norm");
	cout << endl;
	pos_in = pos_iniz;
	cout << "psi_210" << endl;
	Metropolis(*psi_210, pos_in, thr, init, "Gauss", 2.2, "Risultati/psi_210.norm");
	
	cout << endl; 


//distances and blocking method for uncertainties
	vector<double> r100u = Radii("Risultati/psi_100.uniform");
	vector<double> r210u = Radii("Risultati/psi_210.uniform");
	vector<double> r100g = Radii("Risultati/psi_100.norm");
	vector<double> r210g = Radii("Risultati/psi_210.norm");

	block_met(r100u, blocks, "Risultati/r100.uniform");
	block_met(r210u, blocks, "Risultati/r210.uniform");
	block_met(r100g, blocks, "Risultati/r100.norm");
	block_met(r210g, blocks, "Risultati/r210.norm");

	return 0;
}

/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
