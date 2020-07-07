#include "genetic.h"
#include <vector>
#include <algorithm>

// Print
void Print(vector<double> vett) {
	int n = vett.size();

	for(int i=0; i<n; i++)
		cout << vett[i] << endl;
	cout << "Size= " << vett.size() << endl;
}
//function to generate cities on the circle.
//uses the vector x and costruct y
void Circle(vector<double>& x, vector<double>& y, double r, Random rnd) {
	if(x.size() != y.size() ) {
		cerr << "Error: you must use two vectors of the same dimension" << endl;
		return;
	}

	int n= x.size();
	for(int i=0; i<n; i++) {
		double a = rnd.Rannyu(0., 1.); //half +, half - 
		if(a < 0.5)
			y[i] = sqrt(pow(r, 2) - pow(x[i], 2) );

		else
			y[i] = -sqrt(pow(r, 2) - pow(x[i], 2) );
	}
}

void Print(vector<double>& x, string file) {
	
	int n = x.size();
	ofstream of;
	of.open(file);
	for(int i=0; i<n; i++)
	of << x[i] <<endl;	
}


//function to print an output file with two vectors (x,y coord)
void Print(vector<double>& x, vector<double>& y, string file) {
	if(x.size() != y.size() ) {
		cerr << "Error: you must use two vectors of the same dimension" << endl;
		return;
	}
	
	int n = x.size();
	ofstream of;
	of.open(file);
	of << x[0] << "    " << y[0] <<"   "<< endl; //hometown
	for(int i=1; i<n; i++)
	of << x[i] << "    " << y[i] <<endl;	
	//of << x[i] << "    " << y[i] <<"   "<< x[i]*x[i]+y[i]*y[i]<<endl; //verify that third column MUST be r^2, then remove it
	of << x[BC(n,n)] << "    " << y[BC(n,n)] << endl; //hometown
}


//periodic boundary conditions
int BC(int i, int dim) {
	if(i >= dim)
		i = i - dim;

	else if(i < 0)
		i = i + dim;

	return i;
}




// --- INDIVIDUAL CLASS --- 
// Costructor and destructor
Individual::Individual(vector<double>& cromox, vector<double>& cromoy) :
Random(){

	int seed[4];
	int p1, p2;
	ifstream Primes("Primes");

	if (Primes.is_open()){
	        Primes >> p1 >> p2 ;
	  } else cerr << "PROBLEM: Unable to open Primes" << endl;
   	Primes.close();
	
   	ifstream input("seed.in");
   	string property;
   	if (input.is_open()){
     		while ( !input.eof() ){
        	input >> property;
         	if( property == "RANDOMSEED" ){
            	input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            	_rnd.SetRandom(seed,p1,p2);
        	}
     	}
     	input.close();
  	} else cerr << "PROBLEM: Unable to open seed.in" << endl;

_cromox = cromox;
_cromoy = cromoy;

}

Individual::~Individual(){
  _rnd.SaveSeed();
}	

//Functions to call and modify the private methods
vector<double> Individual::getx() {
	return _cromox;
}
vector<double> Individual::gety() {
	return _cromoy;
}
void Individual::setx(vector<double>& cromox) {
	_cromox = cromox;
}
void Individual::sety(vector<double>& cromoy) {
	_cromoy = cromoy;
}


void Individual::Print(string file) {
	if( getx().size() != gety().size() ) {
		cerr << "Error: you must use two vectors of the same dimension in print function!" << endl;
		return;
	}
	
	int n = getx().size();
	ofstream of;
	of.open(file);
	of << getx()[0] << "    " << gety()[0] <<"   "<< endl; //hometown
	for(int i=1; i<n; i++)
	of << getx()[i] << "    " << gety()[i] <<endl;	
	//of << getx()[i] << "    " << getx()[i] <<"   "<< getx()[i]*getx()[i]+gety()[i]*gety()[i]<<endl; //verify that third column MUST be r^2, then remove it
	of << getx()[BC(n,n)] << "    " << gety()[BC(n,n)] << endl; //hometown
return;
}

// --- Check function ---
// It controls that every individual has different cities and the first and the last one are the same
//Useful to check the right working of mutation algorithms
void Individual::Check() {
	
	int dim = getx().size();

	for(int i=0; i<dim; i++)
		for(int j=i+1; j<dim; j++) {
			if( getx()[i] == getx()[j] ) {
				cerr << "Error: test Check False! Same cities in positions " << i << " , " << j << endl;
				return;
			}
			if( gety()[i] == gety()[j] ) {
				cerr << "Error: test Check False! Same cities in positions " << i << " , " << j << endl;
				return;
			}
		}
	//if ( getx()[0]==  getx()[BC(dim,dim)] && gety()[0]== gety()[BC(dim,dim)]) //last city is the first city
	if ( getx()[0]==  getx()[BC(dim,dim)] && gety()[0]==  gety()[BC(dim,dim)]) //last city is the first city
		cout << "Test check: True" << endl;
	else cerr << "Error: test Check False! First and last cities are different and the salesman doesn't return home! " << endl;
	return;
}
 




//--- Permutation function ---

void Individual::Permutation(double prob){
   int dim = ( getx() ).size();
   if (dim != signed(( gety() ).size()) ) {
   	cerr << "Error: vectors of different size in Permutation mutation" << endl;
        return;
   }

   if (prob < 0.025) { //probability of permutation
	int in = _rnd.Rannyu(1., dim);	//
	int fin = _rnd.Rannyu(1., dim); //cromosomes to switch

	vector<double> depx = getx() ;
	swap( depx[in], depx[fin] );
	setx(depx);
	vector<double> depy = gety() ;
	swap( depy[in], depy[fin] );
	sety(depy);
   }
return;
} 
  
//--- Shift function ---
//Shift of lenght "shift" of a block of m cromosomes 
void Individual::Shift(double prob){
   int dim = ( getx() ).size();
   if (dim != signed(( gety() ).size()) ) {
   	cerr << "Error: vectors of different size in Shift mutation" << endl;
        return;
   }

   if (prob < 0.1) { //probability of mutation of 10%
	int in = _rnd.Rannyu(1., dim-1); //first city to shift
	int m = _rnd.Rannyu(1. , dim-in+1); //shift of a block m
	int dim_sh = _rnd.Rannyu(1. , dim-in-m);   //shift lenght

	vector<double> depx = getx() ;
	for(int i =in+m-1 ; i>= in; i--)
		depx.emplace( depx.begin()+ in + m+ dim_sh , depx[i]); //move the train of components of a shift of dim_sh
	depx.erase(depx.begin()+in, depx.begin()+in+m ); //delete the moved components
	setx(depx);

	vector<double> depy = gety() ;
	for(int i =in+m-1 ; i>= in; i--)
		depy.emplace( depy.begin()+ in + m+ dim_sh , depy[i]); //move the train of components of a shift of dim_sh
	depy.erase(depy.begin()+in, depy.begin()+in+m ); //delete the moved components
	sety(depy);
   }
return;
} 


void Individual::Block_Permutation(double prob){
   int dim = ( getx() ).size();
   if (dim != signed(( gety() ).size()) ) {
   	cerr << "Error: vectors of different size in Block_Permutation mutation" << endl;
        return;
   }

   if (prob < 0.1) { //probability of mutation of 10%
	int in = _rnd.Rannyu(1., dim-2);	//first city to permutate
	int bl = _rnd.Rannyu(0., 0.5*(dim-in)-1); //lenght of the block

	vector<double> depx = getx() ;
	for(int i=in; i<=in+bl; i++)	
		swap( depx[BC(i,dim)], depx[BC(i+bl,dim)] );
	setx(depx);

	vector<double> depy = gety() ;
	for(int i=in; i<=in+bl; i++)	
		swap( depy[BC(i,dim)], depy[BC(i+bl,dim)] );
	sety(depy);
   }
return;
} 

void Individual::Inversion(double prob){
   int dim = ( getx() ).size();
   if (dim != signed(( gety() ).size()) ) {
   	cerr << "Error: vectors of different size in Block_Permutation mutation" << endl;
        return;
   }

   if (prob < 0.1) { //probability of mutation of 10%
	int in = _rnd.Rannyu(2., dim-1);	//first city to permutate
	int len = _rnd.Rannyu(0., min(in/2, (dim-in)/2)); //lenght of the block to mirror
	vector<double> depx = getx();
	for(int i=1; i<= len; i++)	
		swap( depx[BC(in+i,dim)], depx[BC(in-i,dim)] );
	setx(depx);
	
	vector<double> depy = gety();
	for(int i=1; i<= len; i++)	
		swap( depy[BC(in+i,dim)], depy[BC(in-i,dim)] );
	sety(depy);
   }
return;
} 


void Individual::Make_Mutations(double prob) {
	if(prob < 0.025) 
		Permutation(prob);
	else if(prob > 0.025 && prob <= 0.050)
		Shift(prob);
	else if(prob > 0.050 && prob <= 0.075)
		Block_Permutation(prob);
	else if(prob > 0.075 && prob <= 0.1)
		Inversion(prob);
}

vector<Individual> Individual::Create_Population(int npop){
   
   int dim = ( getx() ).size();

   int in, fin;

   vector<Individual> depi;
   for(int i=0; i<npop; i++) {
	in = _rnd.Rannyu(1., dim);	//
	fin = _rnd.Rannyu(1., dim); //cromosomes to switch
	
	vector<double> depx = getx() ;
	swap( depx[in], depx[fin] );
	vector<double> depy = gety() ;
	swap( depy[in], depy[fin] );
	setx(depx);
        sety(depy);
	Individual dep(depx,depy);

	depi.push_back(dep);

   } 
cout<< "Creata popolazione di " << depi.size() << " individui" << endl;
return depi; 
} 

// --- L2 calculation
double Individual::L2() {
   int dim = ( getx() ).size();
   double L_2 = 0.;

   for(int i=0; i<dim; i++)
   	L_2 += pow( getx()[ BC(i+1, dim) ]-getx()[i], 2) + pow( gety()[ BC(i+1, dim) ]- gety()[i], 2);
   //cout << L_2<< endl;;
   return L_2;
}

void Crossover(Individual p, Individual m, double prob, Random rnd) { //m , p scelti da roulette truccata
	int dim = (m.getx()).size();
	if (prob > 0.5) {
		int cut = rnd.Rannyu(1., dim-1); //coordinata del taglio
			
		vector<double> depx = p.getx();
		depx.erase(depx.begin()+cut, depx.end());
		vector<double> depy = p.gety();
		depy.erase(depy.begin()+cut, depy.end());
		for(int i=0; i<dim; i++) {
			for(int j=cut; j< dim; j++) { //cicla sulla madre
				if (m.getx()[i] == p.getx()[j]) //se un elemento della madre è uguale a quello che ho tolto al padre
					depx.push_back(m.getx()[i]); //aggiungilo al figlio
				if (m.gety()[i] == p.gety()[j]) //se un elemento della madre è uguale a quello che ho tolto al padre
					depy.push_back(m.gety()[i]); //aggiungilo al figlio
			}	
		}			 
	
	p.setx(depx);
	p.sety(depy);
	}	
}

double L2(vector<double> x, vector<double> y) {
	if(x.size() != y.size()) {
		cerr << "Error: the vectors have different size" << endl;
		return -1.;
	}
	
	int dim = x.size();
	double sum = 0.;

	for(int i=0; i<dim; i++)
		sum += pow( x[ BC(i+1, dim) ]-x[i], 2) + pow( y[ BC(i+1, dim) ]-y[i], 2);

	return sum;
}

void Permutation(vector<double>& x, vector<double>& y, Random rnd, double prob) {
	double dim = x.size();
	if( prob <= 0.1) {
		int in = rnd.Rannyu(1., dim);
		int fin = rnd.Rannyu(1., dim);

		// Don't switch the same city
		while(in == fin)
			fin = rnd.Rannyu(1., dim);

		swap( x[in], x[fin] );
		swap( y[in], y[fin] );
	}
}

void Crossover(vector<double>& x, vector<double>& y, vector<double>& xx, vector<double>& yy, double prob, int pos) {
	if(x.size() != y.size()) {
		cerr << "Error: the vectors have different size" << endl;
		return;
	}

	if(prob <= 0.5) {
		vector<double> depx = x, depy = y, dx = xx, dy = yy;
		
		Write(x, dx, pos);
		Write(xx, depx, pos);
		Write(y, dy, pos);
		Write(yy, depy, pos);
	}
}

// Overwrite (for crossover) from nth position included
void Write(vector<double>& x, vector<double>& y, int n) {
	int dim = x.size();
	vector<double> dep = x;
	int c = n;
	
	for(int i=0; i<dim; i++)
			for(int j=n; j<dim; j++)
				if( dep[j] == y[i] ) {
					x[c] = y[i];
					c++;
				}
}


void New_Gen( vector< vector<double> >& xpop, vector< vector<double> >& ypop, Random rnd, int nmut) { //choose 
		
	int dim_pop = xpop.size();
	int ncity = (xpop[0]).size();
	//cout << dim_pops << "   " << dim_pop << "   " << ncity << endl;
	vector< vector<double> > depx, depy;
	vector<double> d(dim_pop);
	
	//calculate the path lenght and sort the population
	map<int, double> p;

	for(int i=0; i<dim_pop; i++) {
		d[i] = L2(xpop[i], ypop[i]);
		p[i] = d[i];
	}
	
	sort( d.begin(), d.end() );

	double acc = (double)dim_pop/2.;
	int c = 0;
	while(c < dim_pop) {
		int k = acc * pow(rnd.Rannyu(), 2);
		for(int i=0; i<dim_pop; i++)
			if(d[k] == p[i]) {
				depx.push_back(xpop[i]);
				depy.push_back(ypop[i]);
				c++;
				break;
			}
	}
	//cout << "Size after selection  " << depx.size() << endl;
	
	for(int k=0; k<dim_pop-1; k++) {
		// Extract two parents
		int im = rnd.Rannyu(0., (double)dim_pop);
		int id = rnd.Rannyu(0., (double)dim_pop);
		
		while(im == id)
			id = rnd.Rannyu(0., (double)dim_pop);
		
		//crossover and then mutate
		double prob = rnd.Rannyu();
		int pos = rnd.Rannyu(1., (double)(ncity-1));
		Crossover(depx[im], depy[im], depx[id], depy[id], prob, pos);

		for(int i=0; i<nmut; i++) {
			double prob1 = rnd.Rannyu();
			double prob2 = rnd.Rannyu();
			Permutation(depx[im], depy[im], rnd, prob1);
			Permutation(depx[id], depy[id], rnd, prob2);
		}


	// Overwrite the sons in the old population
		xpop[k] = depx[im];
		ypop[k] = depy[im];
		xpop[k+1] = depx[id];
		ypop[k+1] = depy[id];
	}
}

vector<double> L2pop(vector<Individual> pop) {
	int dim_pop= pop.size();
	vector<double> L2;

	for(int i=0; i<dim_pop; i++) {
		L2.push_back(pop[i].L2());
	}
	
	return L2;
}

double Shortest(vector<Individual> pop) {
	int npop = pop.size();
	double min = 10000.;
	vector<double> l = L2pop(pop);
	
	for(int i=0; i<npop; i++)
		if( l[i] < min )
			min = l[i];


	return min;
}
	
// Print average path length (on the best half of the population)
vector<double> Stat(vector<Individual> pop) {
	int npop = pop.size();
	double sum = 0., sum2 = 0.;
	vector<double> ris(2);
	vector<double> d(npop);
	vector<double> l = L2pop(pop);
		
	for(int i=0; i<npop; i++)
		d[i] = l[i];
	
	sort( d.begin(), d.end() );
	
	npop = npop/2;
	for(int i=0; i<npop; i++) {
		sum += d[i];
		sum2 += pow(d[i], 2);
	}
	
	ris[0] = sum/(double)npop;
	ris[1] = sqrt( sum2/(double)npop - ris[0]*ris[0] );

	return ris;
}

void Print_best_path(vector< Individual> pop, string file) {
	int npop = pop.size();
	int b = 0;
	double min = 10000.;
	vector<double> l = L2pop(pop);
	
	for(int i=0; i<npop; i++)
		if( l[i] < min ) {
			min = l[i];
			b = i;
		}
	
	pop[b].Print(file);
}
