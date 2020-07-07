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
#include "mpi.h"

using namespace std;
 
int main (int argc, char *argv[]){

	MPI::Init(argc,argv);	//initialization
	int size = MPI::COMM_WORLD.Get_size();	// how many processes
	int rank = MPI::COMM_WORLD.Get_rank();	// who am I
	MPI_Status stat1, stat2;
	MPI_Request req;	

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
	int ncity = 32;		// Number of cities
	int nmut = 2;	// Number of mutations
	vector<double> x(ncity), y(ncity);
	int npop = 1500; 
	int ngen = 500;
	int nmigr = 10; 	

	double l = 1.;

	if(rank ==0) {
		cout << "---THE TRAVELING SALESMAN PROBLEM---" << endl;
		cout << "Performing a genetic algorithm" << endl << endl;

		cout<< "Traveling around "<<ncity<<" cities"<< endl;
		cout << "Cities IN a square with edge = " << 2.*l << endl << endl;
	}

	// Generate random cities in a square and check
	for(int i=0; i<ncity; i++) {
		x[i] = rnd.Rannyu(-l, l);
		y[i] = rnd.Rannyu(-l, l);
	}

	Print(x, y, "Risultati/sq.0");
	cout << "Initial configuration printed" << endl;

	Individual SquareInd(x,y); //first individual
	//SquareInd.Check();

	// Prepare the vectors for the results
	double x_irecv[size][ncity];
	double y_irecv[size][ncity];
	vector<double> best_irecv(size);


	// Set different random for each node
	int pa[size], pb[size];
	for(int i=0; i<size; i++)
		Primes >> pa[i] >> pb[i];

	rnd.SetRandom(seed, pa[rank], pb[rank]);


	//create the first population (vector of individual)
	vector<Individual> pops = SquareInd.Create_Population(npop);
	int dim_popss= pops.size();
	vector< vector<double> > xpops;
	vector< vector<double> > ypops;	
	vector<double> xb(ncity), yb(ncity);

	for(int i=0; i<dim_popss; i++) {
		xpops.push_back(pops[i].getx());
		ypops.push_back(pops[i].gety());
	}		
	
	//make new generations	
	int n =0; 		
	for(int i=0; i<ngen; i++) {
		
		if ((i+1)%50 == 0 )
			cout << "Generation " << i+1 << "/" << ngen << endl;

		
		if( (i+1)%nmigr == 0 ) {  //migrations every n_migr generations
			n++;
			double* imesgx = new double[ncity]; 
			double* imesgy = new double[ncity];
			int itag=1; int itag2=2; int itag3=3; int itag4=4;
			for(int i=0; i<ncity; i++) {
				imesgx[i] = xb[i];
				imesgy[i] = yb[i];
			}
			if(rank==1){ //swap x-coordinates of the best between rank 0 and 1
				MPI_Isend(&imesgx, ncity, MPI_DOUBLE, 0, itag, MPI_COMM_WORLD, &req);
				MPI_Recv(&imesgx, ncity, MPI_DOUBLE, 0, itag2, MPI_COMM_WORLD, &stat2);
			} 
			else if(rank==0){
				MPI_Send(&imesgx, ncity, MPI_DOUBLE, 1, itag2, MPI_COMM_WORLD);
				MPI_Recv(&imesgx, ncity, MPI_DOUBLE, 1, itag, MPI_COMM_WORLD, &stat1);
			}
			if(rank==1){ //swap y-coordinates of the best between rank 0 and 1
				MPI_Isend(&imesgy, ncity, MPI_DOUBLE, 0, itag, MPI_COMM_WORLD, &req);
				MPI_Recv(&imesgy, ncity, MPI_DOUBLE, 0, itag2, MPI_COMM_WORLD, &stat2);
			} 
			else if(rank==0){
				MPI_Send(&imesgy, ncity, MPI_DOUBLE, 1, itag2, MPI_COMM_WORLD);
				MPI_Recv(&imesgy, ncity, MPI_DOUBLE, 1, itag, MPI_COMM_WORLD, &stat1);
			}
			if(rank==3){ //swap x-coordinates of the best between rank 2 and 3
				MPI_Isend(&imesgx, ncity, MPI_DOUBLE, 2, itag3, MPI_COMM_WORLD, &req);
				MPI_Recv(&imesgx, ncity, MPI_DOUBLE, 2, itag4, MPI_COMM_WORLD, &stat2);
			} 
			else if(rank==2){
				MPI_Send(&imesgx, ncity, MPI_DOUBLE, 3, itag4, MPI_COMM_WORLD);
				MPI_Recv(&imesgx, ncity, MPI_DOUBLE, 3, itag3, MPI_COMM_WORLD, &stat1);
			}
			if(rank==3){ //swap y-coordinates of the best between rank 2 and 3
				MPI_Isend(&imesgy, ncity, MPI_DOUBLE, 2, itag3, MPI_COMM_WORLD, &req);
				MPI_Recv(&imesgy, ncity, MPI_DOUBLE, 2, itag4, MPI_COMM_WORLD, &stat2);
			} 
			else if(rank==2){
				MPI_Send(&imesgy, ncity, MPI_DOUBLE, 3, itag4, MPI_COMM_WORLD);
				MPI_Recv(&imesgy, ncity, MPI_DOUBLE, 3, itag3, MPI_COMM_WORLD, &stat1);
			}
		}

		New_Gen(xpops, ypops, rnd, nmut);

		xb = xpops[0]; //x e y best 
		yb = ypops[0];
	}	
	Print(xb, yb, to_string(rank) + "cities.best"); //print the best of each node
	if (rank ==0 ) cout << n << " migrations done!" << endl;
	// Prepare the results for each node
	double l_isend = L2(xb, yb);
	ofstream lrank; // L2 of each rank
	lrank.open(to_string(rank) + "L2.best");
	lrank << l_isend << endl;
	lrank.close();
	double * x_isend = new double[ncity];
	double * y_isend = new double[ncity];

	for(int i=0; i<ncity; i++) {
		x_isend[i] = xb[i];
		y_isend[i] = yb[i];
	}
		
// Back to node 0
	MPI_Gather(x_isend, ncity, MPI_DOUBLE, x_irecv[rank], ncity, MPI_DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(y_isend, ncity, MPI_DOUBLE, y_irecv[rank], ncity, MPI_DOUBLE, 0, MPI::COMM_WORLD);
	MPI_Gather(&l_isend, 1, MPI_DOUBLE, &best_irecv[rank], 1, MPI_DOUBLE, 0, MPI::COMM_WORLD);


// Print results: best path and correspondent configuration
	if(rank==0) {
		//for(int i=0; i<4; i++)
		//	cout << best_irecv[i] << endl;
		double min = 10000.;
		int ib = 0;
		ofstream ris;

		ris.open("Risultati/Square.ris");
		
		for(int i=0; i<size; i++)
			if( best_irecv[i] < min ) {
				min = best_irecv[i];
				ib = i;
			}

		//ris << min << endl;
		for(int i=0; i<ncity; i++)
			ris << x_irecv[ib][i] << "    " << y_irecv[ib][i] << endl;
		ris << x_irecv[ib][0] << "    " << y_irecv[ib][0] << endl;

		ris.close();
		cout << "Square best path = " << min << endl;
		cout << "Square: results are printed" << endl;
	}

 rnd.SaveSeed();
 MPI::Finalize();// finish
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
