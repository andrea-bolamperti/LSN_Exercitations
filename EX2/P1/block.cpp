#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>
#include <string>

using namespace std;

//Computing statistical uncertainties: the blocking method
//Substantially I translated in c++ the given Python script

void block_met(vector<double>& data, int blocks, string file) {
	int throws = data.size();
	double sum, sum2;
	int ratio = throws/blocks;	//Number of throws in each block
	vector<double> mean, mean2, sum_prog, sum2_prog, err_prog;

	ofstream outfile;
	outfile.open(file);   //file di output
	
	//outfile << "n " << "x " << "err_x" << endl; //if an intestation is desired in the output

	for(int i=0; i< blocks; i++) {
		sum=0.;
		for(int j=0; j < ratio; j++) {
			double k = j+i*ratio;
			sum += data[k];
		}
		mean.push_back(sum/(double)ratio);
		mean2.push_back( pow(mean[i],2) );
	}
	for(int i=0; i< blocks; i++) {
		sum=0.;
		sum2=0.;
		//err_prog[i]=0.;
		for(int j=0; j < i+1; j++) {
			sum += mean[j];
			sum2 += mean2[j];
		}
		sum_prog.push_back(sum/(double)(i+1));
		sum2_prog.push_back(sum2/(double)(i+1));
		err_prog.push_back( sqrt(sum2_prog[i] - pow(sum_prog[i], 2))/sqrt( (double)(i+1)) );
	
		outfile << (i+1)*ratio << "    " << sum_prog[i] << "    " << err_prog[i] << endl;
		
		//outfile << sum_prog[i] << endl;

		//outfile << (i+1)*ratio << "    " << sum/(double)(i+1) << "    "
		//<< sqrt(sum2/(double)(i+1) - pow(sum/(double)(i+1), 2) )/sqrt( (double)(i+1) ) 
		//<< endl;
	}

	outfile.close();
}

		
				
