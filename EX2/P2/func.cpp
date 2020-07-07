#include "func.h"

using namespace std;

double Mean(vector<double> data) {
	double N = data.size();
	int sum = 0;
	for (int i=0; i<N; i++) 
		sum += data[i];

	return sum/N;
}

double devstd(vector<double> data) {
	double N = data.size();
	int sum = 0;
	double mean = Mean(data);

	for (int i=0; i<N; i++) 
		sum += pow(data[i]-mean,2);

	return sqrt(sum/(N-1));
}
	
				
