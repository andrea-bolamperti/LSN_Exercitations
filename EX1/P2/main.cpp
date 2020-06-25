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


// I modified the code from here

// implementation of the second exercise of the first exercitation
	int n = 1E4;
	int N[4]={1,2,10,100};
	double Sd, Se, Sl;
	ofstream dic, ex, lor;
			
 	for(int c=0; c<4; c++) {	//cycle with the cases 1,2,10,100
		dic.open("Risultati/dic_" + to_string(N[c]) + ".dice");
		ex.open("Risultati/exp_" + to_string(N[c]) + ".dice");
		lor.open("Risultati/lorentz_" + to_string(N[c]) + ".dice"); 
		for(int i=0; i<n; i++) {
			Sd=0.;
			Se=0.;
			Sl=0.;
			for(int j=0; j< N[c]; j++) {
				Sd += (int)rnd.Rannyu(1., 7.);
				Se += rnd.Exp(1.);
				Sl += rnd.Lor(1.,0.);
			}
			dic << (double)Sd/(double)N[c] << endl;
			ex << (double)Se/(double)N[c] << endl;
			lor << (double)Sl/(double)N[c] << endl;
		}
		dic.close();
		ex.close();
		lor.close();
	}


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
