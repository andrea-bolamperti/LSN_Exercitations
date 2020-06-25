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

	int throws = 100000;	//number of throws  
	int blocks = 100;	//number of blocks
	//int L = throws/blocks;	//throws in each block
	vector<double> nrs;
	vector<double> nrs_s;

//1.  Code to evaluate the first integral
	for(int i=0; i<throws; i++) 
		nrs.push_back( rnd.Rannyu() );	//generation of a throws number of random numbers
	block_met(nrs, blocks, "Risultati/r.mean");


//2.  Code to evaluate the second integral
	for(int i=0; i<throws; i++){ 
		nrs.push_back( rnd.Rannyu() );	
		nrs_s.push_back(pow(nrs[i]-0.5,2)); //vector of the quantity we want to integrate
	}
	block_met(nrs_s, blocks, "Risultati/r.sigma");

//3.  Chi square
	int M = 100;	// Subintervals
	int n = 10000;
	double exp = (double)n/(double)M;	// Expected
	int dim = 1E6;	// Total throws

	double step = 1./(double)M;
	vector<double> num;
	int count;
	double chi;
	ofstream chi2;

	chi2.open("Risultati/Chi2.test");
   	for(int i=0; i<dim; i++)	//generate dim randoms
		num.push_back( rnd.Rannyu() );
	for(int k=0; k<M; k++) {
		chi=0.;
		for(int i=0; i<M; i++) {	//cycle for counting data in each subinterval
			count = 0;
			for(int j=0; j<n; j++)	//count how many fall in each interval
				if( num[j + k*n] >= step*(double)i && num[j + k*n] < step*(double)(i+1) )
				count++;
			chi += pow( (double)count - exp, 2)/exp;
			
		}
		//chi.push_back(pow( (double)count - exp, 2)/exp);
		chi2 << (k+1) << "    " << chi << endl;	
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
