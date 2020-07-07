#ifndef _genetic__h_
#define _genetic__h_

#include <iostream>
#include <fstream>
#include <map>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>
#include "random.h"

using namespace std;

void Print(vector<double> vett);
void Circle(vector<double>&, vector<double>&, double, Random rnd);  //generates ciries randomly placed on a circumference
void Print(vector<double>& x, string file);
void Print(vector<double>&, vector<double>&, string);
int BC(int, int); //boundary conditions: the saleman wants toreturn home at the end

class Individual : public Random {
  private:
    Random _rnd;
    vector<double> _cromox;	//cromosomes component x, a vector with 32 cities x-coord
    vector<double> _cromoy;	//cromosomes component y, a vector with 32 cities y-coord   
    
  public:
    
    Individual(vector<double>& cromox, vector<double>& cromoy);
    ~Individual();

    vector<double> getx();
    vector<double> gety();
    void setx(vector<double>& cromox);
    void sety(vector<double>& cromoy);
    void Print(string);
    void Check();
    void Permutation(double);
    void Shift(double);
    void Block_Permutation(double);
    void Inversion(double);
    void Make_Mutations(double prob);
    vector<Individual> Create_Population(int);
    double L2();
};

    double L2(vector<double>, vector<double> );
    void Permutation(vector<double>& x, vector<double>& y, Random rnd, double);
    void Crossover(vector<double>& x, vector<double>& y, vector<double>& xx, vector<double>& yy, double prob, int pos);
    void Write(vector<double>& x, vector<double>& y, int n);
    vector<double> L2pop(vector<Individual> pop);
    double Shortest(vector<Individual> pop);	
    vector<double> Stat(vector<Individual> xpop);
    void New_Gen( vector< vector<double> >& xpop, vector< vector<double> >& ypop, Random rnd, int nmut);
    void Print_best_path(vector<Individual> pop, string);


#endif //_genetic__h_
