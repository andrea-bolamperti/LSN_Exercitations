#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "random.h"

using namespace std;

//Computing statistical uncertainties: the blocking method
//Substantially I translated in c++ the given Python script

void block_met(vector<double>& data, int blocks, string file) {
	int throws = data.size();
	double sum, sum2;
	int ratio = throws/blocks;	//Number of throws in each block
	vector<double> mean, mean2, sum_prog, sum2_prog, err_prog;

	ofstream outfile;
	outfile.open(file);   //file di output	
	
	for(int i=0; i< blocks; i++) {
		sum=0.;
		for(int j=0; j < ratio; j++) {
			double k = j+i*ratio;
			sum += data[k];
		}
		mean.push_back(sum/(double)ratio);
		mean2.push_back( pow(mean[i],2) );
	}
	for(int i=0; i< blocks; i++) {
		sum=0.;
		sum2=0.;
		//err_prog[i]=0.;
		for(int j=0; j < i+1; j++) {
			sum += mean[j];
			sum2 += mean2[j];
		}
		sum_prog.push_back(sum/(double)(i+1));
		sum2_prog.push_back(sum2/(double)(i+1));
		err_prog.push_back( sqrt(sum2_prog[i] - pow(sum_prog[i], 2))/sqrt( (double)(i+1)) );
	
		outfile << (i+1)*ratio << "    " << sum_prog[i] << "    " << err_prog[i] << endl;
		
		//outfile << sum_prog[i] << endl;

		//outfile << (i+1)*ratio << "    " << sum/(double)(i+1) << "    "
		//<< sqrt(sum2/(double)(i+1) - pow(sum/(double)(i+1), 2) )/sqrt( (double)(i+1) ) 
		//<< endl;
	}

	outfile.close();
}

//Potential
double V(double x) {
	return pow(x,4) - (2.5) * pow(x,2);
}

//Probability distribution
//Trial wave function
double psi_T(double x, vector<double>& param) {
 	double mu = param[0];
	double sigma = param[1] * param[1];
	
	double g1 = exp(-pow(x-mu,2)/(2.*sigma));
	double g2 = exp(-pow(x+mu,2)/(2.*sigma));
	
	return g1+g2;
}

double psi2_T(double x, vector<double>& param) {
	return pow(psi_T(x, param), 2);
}


// Energy
double En(double x, vector<double>& param) {
	double mu = param[0];
	double sigma = param[1] * param[1];

	double e1 = pow(x-mu, 2)/sigma;
	double e2 = pow(x+mu, 2)/sigma;
	
	double psi_sec = ( -psi_T(x, param) + e1 * exp(-e1/2.) + e2 * exp(-e2/2.) )/sigma;
	
	return -0.5 * psi_sec/psi_T(x, param) + V(x);
}

// Metropolis algorithm
vector<double> Metropolis( double (*distr)(double,vector<double>&), double pos_in, vector<double>& param, int thr, int init, double step, string file) {
	Random rnd;
	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");

	if (Primes.is_open()) {
		Primes >> p1 >> p2 ;
	}
	else cerr << "PROBLEM: Unable to open Primes" << endl;
	
	Primes.close();

	ifstream input("seed.in");
	string property;

	if (input.is_open()) {
		while ( !input.eof() ) {
			input >> property;
			if( property == "RANDOMSEED" ) {
				input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
				rnd.SetRandom(seed,p1,p2);
			}
		}
		input.close();
	}
	else cerr << "PROBLEM: Unable to open seed.in" << endl;
	
	ofstream of;
	of.open(file);

	
// make the steps
	vector<double> ris;
	double passo, new_pos;
	int n = 0;
	double alpha;

	for (int k=0; k< thr; k++) {
		passo = rnd.Rannyu(-step,step);   	
		new_pos = pos_in + passo;
		
		alpha = distr(new_pos,param)/distr(pos_in, param);

		if(alpha>=1.) {	//ACCEPT THE POINT
			pos_in = new_pos;	//subscribe position
			if (k>=init) {	//not consider initial rejected point
				of << pos_in << endl;
				ris.push_back(pos_in);
			}
			n++;
		}
		else {
			double ran = rnd.Rannyu();
			if (ran <= alpha ) { // ACCEPT IT
				pos_in = new_pos;//subscribe position
				if (k>=init){ 	//not consider initial rejected point
					of << pos_in << endl;
					ris.push_back(pos_in);
				}
				n++;
			}
			
			else { //DON'T MOVE IN THE NEXT POINT AND COUNT THE ACTUAL AGAIN
				if (k>=init) {
					of << pos_in << endl;
					ris.push_back(pos_in);
				}
			}
		}
	}

	cout << "Sampling points accepted = " << (double)n/(double)thr * 100. << " %" << endl;
	

	of.close();
	rnd.SaveSeed();

	return ris;		
}



