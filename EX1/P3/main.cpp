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


 // I modified the code from here 


   //for(int i=0; i<20; i++){
   //   cout << rnd.Rannyu() << endl;

	double L=1;	// lenght of the needle
	double d=1.5;	// distance of the horizontal lines
	double x,y,D;	// D is the distance from the center of the needle to the lower horizontal line
	int count,n;
	int thr = 10000;	//number of throws
	int M = 100;	//number of blocks
	vector<double> pi;

	for( int i=0; i<thr; i++) {
		count = 0;
		n=0;
		while(n<=thr){		//cycle until the desired number of generated points
			D = rnd.Rannyu(0., d);	//random number [0,d]
			x = rnd.Rannyu(-1., 1.);	//random numbers x,y to sample an angle [0,pi] (not using pi) 
			y = rnd.Rannyu();
			if(x*x + y*y < 1.) {		
				if(L*y/sqrt(x*x + y*y)/2. + D >= d || L*y/sqrt(x*x + y*y)/2. >= D )	// If cross a line (upper or below)
					count++;		
				
				n++;		
			}
		}
	pi.push_back( 2.*L*(double)thr/( d*(double)count) );
	}
	

	//evaluation of the uncertainties with the blocking method
	block_met(pi,M,"Risultati/Pi.ext");		

	
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
