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

// --- Part 1: a test aplication with casual parameters --- 
	int thr = 11000;
	int init = 1000;
	vector<double> param;

// Initial position and mu, sigma
	param = {1., 0.5};



// sampling... 
	cout << "Sampling trial |psi_T|^2 and <H>" << endl;
	cout << "Using mu = "<<param[0]<< " and sigma = "<<param[1] << endl << endl;
	
	vector<double> psi = Metropolis(psi2_T, 0., param, thr, init, 2.1, "Risultati/psi2_cam.trial"); 

	vector<double> ene;

	for(double elem : psi)
		ene.push_back( En(elem, param) );

	block_met(ene, 100, "Risultati/en.trial");

//-----------------------------

// ----- Part 2: optimization of the parameters

	double sigma_min = 0.3 , sigma_max= 1., sigma_step = 0.01;
	double mu_min = 0.5 , mu_max = 1.2, mu_step = 0.01;
        double E_min = 10000., sigma_opt = sigma_min, mu_opt = mu_min, sum;

	ofstream of;
	of.open("Risultati/var_param.dat");
	ofstream ofopt;
	ofopt.open("Risultati/opt_param.dat");
	int n=1;

	for(double mu = mu_min; mu < mu_max; mu += mu_step) {
		param[0] = mu;
		for(double sigma = sigma_min; sigma < sigma_max; sigma += sigma_step) {
			sum = 0.;
			param[1] = sigma;
			cout << "Simulation "<<n<<"/"<< (mu_max-mu_min)/mu_step * (sigma_max-sigma_min)/sigma_step << endl; 
			psi = Metropolis(psi2_T, 0., param, thr, init, 2.1, "Risultati/Psi.min");
		
			for(double elem : psi)
				sum += En(elem, param);
			
			sum = sum/(double)psi.size();
			
			if(sum < E_min) {
				E_min = sum;
				sigma_opt = sigma;
				mu_opt = mu;
			}
			n++;
			of << mu << "    " << sigma << "    " << sum << endl;
		}
	}
	ofopt << "mu_opt" <<"   "<< "sigma_opt" << "   "<< "E_min" << endl;
	ofopt << mu_opt <<"   "<< sigma_opt << "   "<< E_min << endl;
	of.close();
	ofopt.close();

	cout << "Minimized parameters:" << endl;
	cout << "E = " << E_min << " mu = " << mu_opt << " sigma = " << sigma_opt << endl << endl;


// perform with the best couple of parameters the simulation and evalue the uncertainties
	param = {mu_opt, sigma_opt};

	psi = Metropolis(psi2_T, 0., param, thr, init, 2.1, "Risultati/Psi.min");
	
	vector<double> ene_opt;

	for(double elem : psi)
		ene_opt.push_back( En(elem, param) );

	block_met(ene_opt, 100, "Risultati/En.min");
	

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
