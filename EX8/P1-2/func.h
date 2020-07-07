#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "random.h"

using namespace std;

void block_met(vector<double>&, int, string);

double psi_T(double, vector<double>&);
double psi2_T(double, vector<double>&);

double V(double);
double En(double, vector<double>&);

vector<double> Metropolis(double (*distr)(double, vector<double>&), double pos_in, vector<double>& param, int thr, int init, double step, string file);

