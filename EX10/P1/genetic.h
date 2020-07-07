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
void Circle(vector<double>&, vector<double>&, double, Random rnd);
void Print(vector<double>& x, string file);
void Print(vector<double>&, vector<double>&, string);
int BC(int, int);

class Individual : public Random {
  private:
    Random _rnd;
    vector<double> _cromox;	//cromosomes component x, a vector with 32 cities x-coord
    vector<double> _cromoy;	//cromosomes component y, a vector with 32 cities y-coord   
    double _L2;

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
    double L2();
};


void Permutation(vector<double>& x, vector<double>& y, Random rnd, double);
void SA(vector<double>& x, vector<double>& y, vector<double> depx, vector<double> depy, double prob, double beta, int& thr);
double L2(vector<double> x, vector<double> y);

#endif //_genetic__h_
