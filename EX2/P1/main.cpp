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
   
   int M = 1E5; //throws
   int blocks = 100; //blocks
   vector<double> f_u(M);
   vector<double> f_i(M);
   ofstream of;
   

	of.open("Risultati/sampled.points");
        of << "sampled_points" << endl; 

	for(int i=0; i<M; i++) 	//evalue the integral in M throws
		f_u[i] = M_PI/2. * cos(M_PI * rnd.Rannyu()/2.);
	
	block_met(f_u, blocks, "Risultati/int.uniform");	
   


// Importance sampling
//I choose p= (4/3)(1-x/2)
   	for(int i=0; i<M; i++) {	//evalue the integral in n points
		double x = rnd.pcos();
		of << x << endl;
		f_i[i] = 3.*M_PI/8. * cos(M_PI * x/2.) / (1.-x/2.);
	}
	block_met(f_i, blocks, "Risultati/int.impsam");	


	of.close();
   		

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
