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
#include "block.h"

using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   int seed[4];
   int p1, p2;
   ifstream Primes("Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } 
   else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } 
   else cerr << "PROBLEM: Unable to open seed.in" << endl;




   int nrw= 1E4;	//number of random-walk simulations
   int steps = 1E2;	//number of steps
   int blocks = 100;	//number of blocks
   double x, y, z, num;
   vector<double> r_d(nrw);	
   ofstream risd, risc;

	risd.open("Risultati/RW.discr");
	risc.open("Risultati/RW.cont");

// Discrete random walk on a lattice with lenght "a"
   int a_d = 1;	//lenght of the lattice

   for(int k=1; k<=steps; k++) {	//cycle to print result as a function of steps
   	for(int i=0; i< nrw; i++) {
   		x = 0.;
		y = 0.;		//initialize the position 
		z = 0.;
		for (int j=0; j< k; j++) {
			num = rnd.Rannyu(0.,6.);	//random number [0,6]
			if ((int)num % 6 == 0)
				x -= a_d; 
			else if ((int)num % 6 == 1)
				y -= a_d; 
			else if ((int)num % 6 == 2)
				z -= a_d; 		
			else if ((int)num % 6 == 3)
				x += a_d; 
			else if ((int)num % 6 == 4)
				y += a_d; 
			else	//((int)num % 6 = 5)
				z += a_d;
		}
		r_d[i] = pow(x, 2) + pow(y, 2) + pow(z, 2);
	}
	
	block_met(r_d, blocks, "Risultati/RWb.discr" );
	risd << k <<" "<< sqrt(Mean(r_d)) <<" "<< sqrt(devstd(r_d)) << endl;
   }
	
   risd.close();
// Continuum Random Walk
   double th, phi;	//random angles of direction
   vector<double> r_c(nrw);	
   double a_c = 1;	//lenght of the "jump"  
	
   for(int k=1; k<=steps; k++) {	//cycle to print result as a function of steps
   	for(int i=0; i< nrw; i++) {
   		x = 0.;
		y = 0.;		//initialize the position 
		z = 0.;
		for (int j=0; j< k; j++) {
			th = rnd.Rannyu(0., M_PI);
			phi = rnd.Rannyu(0., 2.*M_PI);

			x += a_c * sin(th) * cos(phi);
			y += a_c * sin(th) * sin(phi);
			z += a_c * cos(th);
		}
	r_c[i] = pow(x, 2) + pow(y, 2) + pow(z, 2);
	}
	block_met(r_c, blocks, "Risultati/RWb.cont" );
	risc << k <<" "<< sqrt(Mean(r_c)) <<" "<< sqrt(devstd(r_c)) << endl;
   }
	
   risd.close();
			
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
