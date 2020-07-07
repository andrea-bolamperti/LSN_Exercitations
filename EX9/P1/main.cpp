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
	cout << "Performing a genetic algorithm" << endl << endl;

	int ncity = 32;		// Number of cities
	cout<< "Traveling around "<<ncity<<" cities"<< endl;
	int nmut = 14;	// Number of mutations
	vector<double> x(ncity), y(ncity);
	double r = 1.;	//radius of the circle
	//int i = 1;

	cout<<"Cities placed ON a circle of radius "<<r<< endl <<endl;

	x[0] = rnd.Rannyu(-1., 1); //the first city

	for(int i=1; i<ncity; i++)  //generate x coordinates of the cities
		x[i] = rnd.Rannyu(-1.,1.);
	
	Circle(x, y, r, rnd); // generate y coord to be ON a circle

	Print(x, y, "Risultati/circ.0");
	cout << "Initial configuration printed" << endl << endl;
	

	int npop = 1000; 
	int ngen =  500;

	Individual CircleInd(x,y); //first individual
	CircleInd.Check();

	//create the first population (vector of individual)
	vector<Individual> popc = CircleInd.Create_Population(npop);
	int dim_popss= popc.size();
	vector< vector<double> > xpopc;
	vector< vector<double> > ypopc;
	for(int i=0; i<dim_popss; i++) {
		xpopc.push_back(popc[i].getx());
		ypopc.push_back(popc[i].gety());
	}	
	
	//make new generations
	ofstream bp, al;
	bp.open("Risultati/circ.shortest");		// File with best path for generation
	al.open("Risultati/circ.ave");		// File with average length
			
	for(int i=0; i<ngen; i++) {
		
		if ((i+1)%50 == 0 )
			cout << "Generation " << i+1 << "/" << ngen << endl;

	
		New_Gen(xpopc, ypopc, rnd, nmut);

		for(int i=0; i<npop; i++) {
			popc[i].setx(xpopc[i]);
			popc[i].sety(ypopc[i]);
		}
			
		for(int i=0; i<nmut; i++) {//probability of mutation
			double probm = rnd.Rannyu(); 
			int pos = rnd.Rannyu(0, npop);
			if(probm < 0.025) 
				popc[pos].Permutation(probm);
			else if(probm > 0.025 && probm <= 0.05)
				popc[pos].Shift(probm);
			else if(probm > 0.05 && probm <= 0.075)
				popc[pos].Block_Permutation(probm);
			else if(probm > 0.075 && probm <= 0.1)
				popc[pos].Inversion(probm);
		}		
		


		bp << i+1 << "    " << Shortest(popc) << endl;
		//cout << "shortest fatto" << endl;
		al << i+1 << "    " << (Stat(popc))[0] << "    " << (Stat(popc))[1] << endl;
		//cout << "stat fatto" << endl;
	}

	bp.close();
	al.close();
	cout << endl;

	Print_best_path(popc, "Risultati/circ.fin");

	cout << "Optimized configuration printed" << endl << endl;
	cout << "-------------------------------" << endl << endl; 	

//SQUARE
	rnd.SetRandom(seed,p1,p2); //or I will have different cities when I repeat	

	cout << "Performing a genetic algorithm" << endl << endl;
	cout<< "Traveling around "<<ncity<<" cities"<< endl;
	double l = 1.;
	cout << "Cities IN a square with edge = " << 2.*l << endl << endl;

	// Generate random cities in a square and check
	for(int i=0; i<ncity; i++) {
		x[i] = rnd.Rannyu(-l, l);
		y[i] = rnd.Rannyu(-l, l);
	}

	Print(x, y, "Risultati/sq.0");
	cout << "Initial configuration printed" << endl << endl;

	Individual SquareInd(x,y); //first individual
	SquareInd.Check();


	//create the first population (vector of individual)
	vector<Individual> pops = SquareInd.Create_Population(npop);
	vector< vector<double> > xpops;
	vector< vector<double> > ypops;	
	
	for(int i=0; i<dim_popss; i++) {
		xpops.push_back(pops[i].getx());
		ypops.push_back(pops[i].gety());
	}	
	
	//make new generations
	ofstream bps, als;
	bps.open("Risultati/sq.shortest");		// File with best path for generation
	als.open("Risultati/sq.ave");		// File with average length
			
	for(int i=0; i<ngen; i++) {
		
		if ((i+1)%50 == 0 )
			cout << "Generation " << i+1 << "/" << ngen << endl;

	
		New_Gen(xpops, ypops, rnd, nmut);

		for(int i=0; i<npop; i++) {
			pops[i].setx(xpops[i]);
			pops[i].sety(ypops[i]);
		}
			
		for(int i=0; i<nmut; i++) {//probability of mutation
			double probm = rnd.Rannyu(); 
			int pos = rnd.Rannyu(0, npop);
			if(probm < 0.025) 
				pops[pos].Permutation(probm);
			else if(probm > 0.025 && probm <= 0.05)
				pops[pos].Shift(probm);
			else if(probm > 0.05 && probm <= 0.075)
				pops[pos].Block_Permutation(probm);
			else if(probm > 0.075 && probm <= 0.1)
				pops[pos].Inversion(probm);
		}		
		


		bps << i+1 << "    " << Shortest(pops) << endl;
		//cout << "shortest fatto" << endl;
		als << i+1 << "    " << (Stat(pops))[0] << "    " << (Stat(pops))[1] << endl;
		//cout << "stat fatto" << endl;
	}

	bps.close();
	als.close();
	cout << endl;

	Print_best_path(pops, "Risultati/sq.fin");

	cout << "Optimized configuration printed" << endl << endl << endl;


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
