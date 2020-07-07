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

//Probability distributions
//Wave function in a0 units
double psi_100(vector<double>& coord) {
	double r = 0.;
	for(double elem : coord)
		r += elem * elem;
	r = sqrt(r);
	
	return pow(exp(-r)/sqrt(M_PI) ,2);
}

double psi_210(vector<double>& coord) {
	double r = 0.;
	for(double elem : coord)
		r += elem * elem;
	r = sqrt(r);
	double cos_th = coord[2]/r;	//cos_th = z / r
	return pow(r* cos_th * exp(-r/2.)/(8.*sqrt(M_PI/2.)) ,2);
}


// Metropolis algorithm
void Metropolis( double (*distr)(vector<double>&), vector<double>& pos_in, int thr, int init, string prob, double step, string file) {
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

	if(prob == "Unif")
		cout << "Metropolis algorithm with uniform probability transition" << endl;

	else if(prob == "Gauss")
		cout << "Metropolis algorithm with normal probability transition" << endl;

	else {
		cerr << "Error: you must choose 'Unif' or 'Gauss'" << endl;
		return;
	}

	of << "x " << "y " << "z" << endl;	//inestation

// make the steps
	int dim = pos_in.size();
	vector<double> passo(dim), new_pos(dim);
	int n = 0;
	double alpha;

	for (int k=0; k< thr; k++) {
		for (int i = 0; i<dim; i++) {
			if(prob == "Unif")
				passo[i] = rnd.Rannyu(-step, step);
			else
				passo[i] = rnd.Gauss(0., step);

			new_pos[i]= pos_in[i] + passo[i];
		}
		
		alpha = distr(new_pos)/distr(pos_in);

		if(alpha>=1.) {	//ACCEPT THE POINT
			for(int i=0; i<dim; i++) {
				pos_in[i] = new_pos[i];	//subscribe position
				if (k>=init) 	//not consider initial rejected point
					of << pos_in[i] << " ";
			}
			if (k>=init)
				of << endl;
			n++;
		}
		else {
			double ran = rnd.Rannyu();
			if (ran <= alpha ) { // ACCEPT IT
				for(int i=0; i<dim; i++) {
					pos_in[i] = new_pos[i];	//subscribe position
					if (k>=init) 	//not consider initial rejected point
						of << pos_in[i] << " ";
				}
				if (k>=init)
					of << endl;
				n++;
			}
			
			else { //DON'T MOVE IN THE NEXT POINT AND COUNT THE ACTUAL AGAIN!
				for(int i=0; i<dim; i++) {
					if (k>=init)
						of << pos_in[i] << " ";
				}
				if (k>=init)
					of << endl;
			}
		}
	}

cout << "Sampling points accepted = " << (double)n/(double)thr * 100. << " %" << endl;
	cout << "The first " << init << " have been rejected" << endl;

	of.close();
	rnd.SaveSeed();
}

// Take a file with coordinates and compute the distances from the origin
vector<double> Radii(string infile) {
	ifstream inf;
	inf.open(infile);

	double x,y,z;
	vector<double> r;

	string Intestation;	//ignore the first line
	getline(inf, Intestation);
	
	while(inf>> x >> y >> z)
		r.push_back ( sqrt(x*x + y*y + z*z) );

	inf.close();

	return r;
}			
