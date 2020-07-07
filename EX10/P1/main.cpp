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
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <fstream>
#include "genetic.h"
#include "random.h"

using namespace std;
 
int main (int argc, char *argv[]){

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


//The salesman problem
	cout << "---THE TRAVELING SALESMAN PROBLEM---" << endl;
	cout << "Performing a simulated annealing algorithm" << endl << endl;

	int ncity = 32;		// Number of cities
	int nmut = 14;	// Number of mutations
	cout<< "Traveling around "<<ncity<<" cities"<< endl;
	vector<double> x(ncity), y(ncity);
	double r = 1.;	//radius of the circle

	cout<<"Cities placed ON a circle of radius "<<r<< endl <<endl;

	x[0] = rnd.Rannyu(-1., 1); //the first city

	for(int i=1; i<ncity; i++)  //generate x coordinates of the cities
		x[i] = rnd.Rannyu(-1.,1.);
	
	Circle(x, y, r, rnd); // generate y coord to be ON a circle

	Print(x, y, "Risultati/circ.0");
	cout << "Initial configuration printed" << endl << endl;
	

//look for the minimum
	double bi =0., bf = 50., db = 0.0005;
	int thr, n = 1000;

	Individual Ind_circle(x,y);

	vector<double> depx(ncity), depy(ncity);
	ofstream bp;

	bp.open("Risultati/circ.shortest");		// File with best path
	cout << "Beta in range = [ " << bi << ";" << bf << " ] with step = " << db << endl;
	cout << "Number of MC steps = " << n << endl << endl;

	for(double beta=bi; beta<bf; beta=beta+db) {
		thr = 0;
		for(int i=0; i<n; i++) {
			depx = x;
			depy = y;
			// Find a new configuration
			for(int j=0; j<nmut; j++) {
				double prob = rnd.Rannyu();
				Permutation(depx, depy, rnd, prob);
			}

			double prob = rnd.Rannyu();
			// Substitute if better

			SA(x, y, depx, depy, prob, beta, thr);
			
			Ind_circle.setx(x);
			Ind_circle.sety(y);
		}

		bp << beta << "    " << Ind_circle.L2() << endl;
		
		if((int)(beta+1) % 10 == 0) {
			cout << "Beta = " << beta+1 << "/" << bf << endl;
			cout << "Acceptance ratio = " << (double)thr/(double)n * 100. << endl;
		}
	}
		
	Print(x, y, "Risultati/circ.fin");
	cout << "Optimized configuration printed" << endl << endl;
	bp.close();
	cout << "-------------------------------------" << endl << endl;

//SQUARE
	rnd.SetRandom(seed,p1,p2);
	double l = 1.;

	cout<<"Cities placed IN a square of edge "<< 2*l << endl <<endl;

	for(int i=0; i<ncity; i++){  //generate x coordinates of the cities
		x[i] = rnd.Rannyu(-1.,1.);
		y[i] = rnd.Rannyu(-1.,1.);
	}
	
	Print(x, y, "Risultati/sq.0");
	cout << "Initial configuration printed" << endl << endl;
	

//look for the minimum
	Individual Ind_sq(x,y);

	vector<double> depxs(ncity), depys(ncity);
	ofstream bps;

	bps.open("Risultati/sq.shortest");		// File with best path
	cout << "Beta in range = [ " << bi << ";" << bf << " ] with step = " << db << endl;
	cout << "Number of MC steps = " << n << endl << endl;

	for(double beta=bi; beta<bf; beta=beta+db) {
		thr = 0;
		for(int i=0; i<n; i++) {
			depxs = x;
			depys = y;
			// Find a new configuration
			for(int j=0; j<nmut; j++) {
				double prob = rnd.Rannyu();
				Permutation(depxs, depys, rnd, prob);
			}

			double prob = rnd.Rannyu();
			// Substitute if better

			SA(x, y, depxs, depys, prob, beta, thr);
			
			Ind_sq.setx(x);
			Ind_sq.sety(y);	
		}

		bps << beta << "    " << Ind_sq.L2() << endl;
		
		if((int)(beta+1) % 10 == 0) {
			cout << "Beta = " << beta+1 << "/" << bf << endl;
			cout << "Acceptance ratio = " << (double)thr/(double)n * 100. << endl;
		}
	}
		
	Print(x, y, "Risultati/sq.fin");
	cout << "Optimized configuration printed" << endl << endl;
	bps.close();
	cout << "-------------------------------------" << endl << endl;




rnd.SaveSeed();
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
