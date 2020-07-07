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


// Parameters

   double S_0 = 100.;
   double s;
   double T = 1.;
   double K = 100.;
   double r = 0.1;
   double sigma = 0.25;  

   int M = 1E5;
   int bl = 100;
   vector<double> C, P;

//Sampling directly the final asset price
   for (int i=0; i<M; i++) {
   	s = S_0 * exp((r-0.5*pow(sigma,2))*T + sigma*rnd.Gauss(0., T));
	C.push_back( exp(-r*T) * max(s-K, 0.) );	//vector of 'calls(T)'
	P.push_back( exp(-r*T) * max(K-s, 0.) );	//vector of 'puts(T)'
   }
   block_met(C, bl, "Risultati/Call.dir");
   block_met(P, bl, "Risultati/Put.dir");


//Sampling the discretized path of the asset price
   int steps=100;	//number of discretized times
   double t0=0.;	//initial time	
 
   double it =( (T-t0)/(double)steps );	//lenght of all the intervals of time

   for (int i=0; i<M; i++) {
	s=S_0;
   	for (double t=t0; t<=T; t += it) {	//don't define a vector S[i] because I'm only interested in the final value S[T].
		s *= exp( (r - 0.5 * pow(sigma, 2)) * it + sigma * rnd.Gauss(0., 1.) * sqrt(it) );
	}
	C[i]=( exp(-r*T) * max(s-K, 0.) );	//vector of 'calls(T)'
	P[i]=( exp(-r*T) * max(K-s, 0.) );	//vector of 'puts(T)'
   }	
   block_met(C, bl, "Risultati/Call.disc");
   block_met(P, bl, "Risultati/Put.disc");


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
