/****************************************************************
*****************************************************************
    _/    _/  _/_/_/  _/       Numerical Simulation Laboratory
   _/_/  _/ _/       _/       Physics Department
  _/  _/_/    _/    _/       Universita' degli Studi di Milano
 _/    _/       _/ _/       Prof. D.E. Galli
_/    _/  _/_/_/  _/_/_/_/ email: Davide.Galli@unimi.it
*****************************************************************
*****************************************************************/
#include <stdlib.h>     // srand, rand: to generate random number
#include <iostream>     // cin, cout: Standard Input/Output Streams Library
#include <fstream>      // Stream class to both read and write from/to files.
#include <cmath>        // rint, pow
#include <string>
#include <vector>
#include "MolDyn_NVE.h"
#include "block.h"

using namespace std;

int main(){ 
  if(mode == "actual") {			//Inizialization
  	cout << "Reading the old configuration from config.0" << endl;
  	Input();
  }
  else if(mode == "old") {
  	cout << "Reading the old configuration from old.0" << endl;             
  	Input_old();
  }
  int nconf = 1;
  
  Eq_Temp();	//equilibration of the temperature
  
  for(int istep=1; istep <= nstep; ++istep){
     Move();           //Move particles with Verlet algorithm
     if(istep%iprint == 0) cout << "Number of time-steps: " << istep << endl;
     if(istep%10 == 0){
        Measure();     //Properties measurement
//        ConfXYZ(nconf);//Write actual configuration in XYZ format //Commented to avoid "filesystem full"! 
        nconf += 1;
     }
  }
  Measure_ave();	//block method
  ConfFinal();         //Write final configuration to restart

  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl << endl;
  ReadConf.open("config.0");
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//Prepare initial velocities
   cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   double sumv[3] = {0.0, 0.0, 0.0};
   for (int i=0; i<npart; ++i){
     vx[i] = rand()/double(RAND_MAX) - 0.5;
     vy[i] = rand()/double(RAND_MAX) - 0.5;
     vz[i] = rand()/double(RAND_MAX) - 0.5;

     sumv[0] += vx[i];
     sumv[1] += vy[i];
     sumv[2] += vz[i];
   }
   for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
   double sumv2 = 0.0, fs;
   for (int i=0; i<npart; ++i){
     vx[i] = vx[i] - sumv[0];
     vy[i] = vy[i] - sumv[1];
     vz[i] = vz[i] - sumv[2];

     sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   }
   sumv2 /= (double)npart;

   fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
   for (int i=0; i<npart; ++i){
     vx[i] *= fs;
     vy[i] *= fs;
     vz[i] *= fs;

     xold[i] = Pbc(x[i] - vx[i] * delta);
     yold[i] = Pbc(y[i] - vy[i] * delta);
     zold[i] = Pbc(z[i] - vz[i] * delta);
   }
   return;
}




void Input_old(void){ //Prepare all stuff for the simulation
  ifstream ReadInput,ReadConf;
  //double ep, ek, pr, et, vir;

  cout << "Classic Lennard-Jones fluid        " << endl;
  cout << "Molecular dynamics simulation in NVE ensemble  " << endl << endl;
  cout << "Interatomic potential v(r) = 4 * [(1/r)^12 - (1/r)^6]" << endl << endl;
  cout << "The program uses Lennard-Jones units " << endl;

  seed = 1;    //Set seed for random numbers
  srand(seed); //Initialize random number generator
  
  ReadInput.open("input.dat"); //Read input

  ReadInput >> temp;

  ReadInput >> npart;
  cout << "Number of particles = " << npart << endl;

  ReadInput >> rho;
  cout << "Density of particles = " << rho << endl;
  vol = (double)npart/rho;
  cout << "Volume of the simulation box = " << vol << endl;
  box = pow(vol,1.0/3.0);
  cout << "Edge of the simulation box = " << box << endl;

  ReadInput >> rcut;
  ReadInput >> delta;
  ReadInput >> nstep;
  ReadInput >> iprint;
  ReadInput >> mode;
  ReadInput >> Tst;
  
  if(mode != "actual" && mode != "old") {
		cerr << "Error: you should choose 'old' or 'actual' in the reading file mode..." << endl << endl;
		return;
	}
   

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  cout << "The mode of the input choosen is the '" << mode << "' mode" << endl << endl;
  ReadInput.close();

//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

//Read initial configuration
  cout << "Read initial configuration from file config.0 " << endl;
  
  ReadConf.open("config.0");
  
  for (int i=0; i<npart; ++i){
    ReadConf >> x[i] >> y[i] >> z[i];
    x[i] = x[i] * box;
    y[i] = y[i] * box;
    z[i] = z[i] * box;
  }
  ReadConf.close();

//if you want to use the old configuration
  if(mode == "old") {
  	cout << "Reading the old configuration from old.0" << endl;
  	
  	ReadConf.open("old.0");
  	for (int i=0; i<npart; ++i){
			ReadConf >> xold[i] >> yold[i] >> zold[i];
			xold[i] = xold[i] * box;
			yold[i] = yold[i] * box;
			zold[i] = zold[i] * box;
		}

		ReadConf.close();
  cout << endl;
  return;
  }
 

   else if(mode == "actual") {

   //Prepare initial velocities
   	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl << endl;
   	double sumv[3] = {0.0, 0.0, 0.0};
   	for (int i=0; i<npart; ++i){
     		vx[i] = rand()/double(RAND_MAX) - 0.5;
     		vy[i] = rand()/double(RAND_MAX) - 0.5;
     		vz[i] = rand()/double(RAND_MAX) - 0.5;

	     	sumv[0] += vx[i];
	     	sumv[1] += vy[i];
	     	sumv[2] += vz[i];
  	}
	for (int idim=0; idim<3; ++idim) sumv[idim] /= (double)npart;
	   
	double sumv2 = 0.0, fs;
	   
	for (int i=0; i<npart; ++i){
        	vx[i] = vx[i] - sumv[0];
	     	vy[i] = vy[i] - sumv[1];
	     	vz[i] = vz[i] - sumv[2];
	
	     	sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
	}
	sumv2 /= (double)npart;
	
	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor 
	
	for (int i=0; i<npart; ++i){
		vx[i] *= fs;
	     	vy[i] *= fs;
	    	vz[i] *= fs;
	
	     	xold[i] = Pbc(x[i] - vx[i] * delta);
	     	yold[i] = Pbc(y[i] - vy[i] * delta);
	     	zold[i] = Pbc(z[i] - vz[i] * delta);
	}
   cout << endl;
   return;
   }
}

//Equilibration of the temperature
void Eq_Temp(void){
	cout <<"Equilibration of temperature..." << endl << endl;
	
	double sumv[3];
	double dT = 0.0001;
	int n=0;	//accuracy of equilibration
	ofstream of;

        of.open("Risultati/Eq.temp",ios::app);
	do {
		for(int i=0; i<3; i++)
			sumv[i] = 0.;
		
		//Make a step to evaluate the velocity
			Move();
		
		for(int i=0; i<npart; i++) {
			v0x[i] = Pbc (x[i]-xold[i])/delta;
			v0y[i] = Pbc (y[i]-yold[i])/delta;// v(t + dt/2)
			v0z[i] = Pbc (z[i]-zold[i])/delta;
		}
		for (int i=0; i<npart; ++i){
     		
	     		sumv[0] += v0x[i];
	     		sumv[1] += v0y[i];
	     		sumv[2] += v0z[i];
  		}
		for (int idim=0; idim<3; ++idim) 
			sumv[idim] /= (double)npart;
	   
		double sumv2 = 0.0, eqf;
	   
		for (int i=0; i<npart; ++i){
        		v0x[i] = v0x[i] - sumv[0];
	     		v0y[i] = v0y[i] - sumv[1];
	     		v0z[i] = v0z[i] - sumv[2];
	
	     	sumv2 += v0x[i]*v0x[i] + v0y[i]*v0y[i] + v0z[i]*v0z[i];
		}
		sumv2 /= (double)npart;

		temp = sumv2/3. ;	//evaluation of temperature
		of << n << "    " << temp << endl;
		n++;
		if(n == 10000) {
			cout << "Warning! Can't reach the desired temperature" << endl;
			return;
		}
		eqf = sqrt(Tst/temp); 	//scaling factor
		for (int i=0; i<npart; ++i){
		v0x[i] *= eqf;
	     	v0y[i] *= eqf;
	    	v0z[i] *= eqf;
	
	     	xold[i] = Pbc(x[i] - v0x[i] * delta);
	     	yold[i] = Pbc(y[i] - v0y[i] * delta);
	     	zold[i] = Pbc(z[i] - v0z[i] * delta);
		}
	} while(abs(Tst - temp) >= dT);


	cout << "Equilbrium temperature: " << temp << "  --- Desired temperature: " << Tst << endl;
	if(abs(temp-Tst)<dT)
		cout << "The system reached the equilibrium temperature "<<temp<<" in " <<n<< " steps" << " and ready to perform" << endl << endl;
	cout << "This value will be used in the simulation" << endl;
	cout << "Velocities are scaled and positions saved" << endl << endl;
	
	

	of.close();
	return;

}

void Move(void){ //Move particles with Verlet algorithm
  double xnew, ynew, znew, fx[m_part], fy[m_part], fz[m_part];

  for(int i=0; i<npart; ++i){ //Force acting on particle i
    fx[i] = Force(i,0);
    fy[i] = Force(i,1);
    fz[i] = Force(i,2);
  }

  for(int i=0; i<npart; ++i){ //Verlet integration scheme
    
    dxold[i] = xold[i];		//to save old positions
    dyold[i] = yold[i];		//t-dt
    dzold[i] = zold[i];

    xnew = Pbc( 2.0 * x[i] - dxold[i] + fx[i] * pow(delta,2) );
    ynew = Pbc( 2.0 * y[i] - dyold[i] + fy[i] * pow(delta,2) );	//t + dt
    znew = Pbc( 2.0 * z[i] - dzold[i] + fz[i] * pow(delta,2) );

    vx[i] = Pbc(xnew - xold[i])/(2.0 * delta);
    vy[i] = Pbc(ynew - yold[i])/(2.0 * delta);	// v
    vz[i] = Pbc(znew - zold[i])/(2.0 * delta);

    xold[i] = x[i];
    yold[i] = y[i];	// t 
    zold[i] = z[i];

    x[i] = xnew;
    y[i] = ynew;	//t + dt
    z[i] = znew;
  }
  return;
}

double Force(int ip, int idir){ //Compute forces as -Grad_ip V(r)
  double f=0.0;
  double dvec[3], dr;

  for (int i=0; i<npart; ++i){
    if(i != ip){
      dvec[0] = Pbc( x[ip] - x[i] );  // distance ip-i in pbc
      dvec[1] = Pbc( y[ip] - y[i] );
      dvec[2] = Pbc( z[ip] - z[i] );

      dr = dvec[0]*dvec[0] + dvec[1]*dvec[1] + dvec[2]*dvec[2];
      dr = sqrt(dr);

      if(dr < rcut){
        f += dvec[idir] * (48.0/pow(dr,14) - 24.0/pow(dr,8)); // -Grad_ip V(r)
      }
    }
  }
  
  return f;
}

void Measure(){ //Properties measurement
  //int bin;
  double v, t, vij;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("Risultati/output_epot.dat",ios::app);
  Ekin.open("Risultati/output_ekin.dat",ios::app);
  Temp.open("Risultati/output_temp.dat",ios::app);
  Etot.open("Risultati/output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
//Potential energy
       v += vij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; 	//Potential energy per particle
    stima_kin = t/(double)npart; 	//Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; 	//Total energy per particle

    Epot << stima_pot  << endl;
    Ekin << stima_kin  << endl;
    Temp << stima_temp << endl;
    Etot << stima_etot << endl;

//save the current values in a vector needed to evaluate the means with blocks method in the following function
    utot.push_back(stima_pot);
    ktot.push_back(stima_kin);
    Ttot.push_back(stima_temp);
    etot.push_back(stima_etot);

    Epot.close();
    Ekin.close();
    Temp.close();
    Etot.close();

    return;
}

void Measure_ave() {
	block_met(utot, blocks, "Risultati/ave_epot.out");
	block_met(ktot, blocks, "Risultati/ave_ekin.out");
	block_met(Ttot, blocks, "Risultati/ave_temp.out");
	block_met(etot, blocks, "Risultati/ave_etot.out");
}


void ConfFinal(void){ //Write final configuration
  ofstream WriteConf;
  ofstream WriteConf_old;

  cout << endl;
  cout << "Print final configuration to file config.final " << endl;
  WriteConf.open("config.final");
  
  cout << "Print old final configuration to file old.final " << endl << endl;
  WriteConf_old.open("old.final");

  for (int i=0; i<npart; ++i){
    WriteConf << x[i]/box << "   " <<  y[i]/box << "   " << z[i]/box << endl;
  }
  WriteConf.close();

  for (int i=0; i<npart; ++i){
    WriteConf_old << xold[i]/box << "   " <<  yold[i]/box << "   " << zold[i]/box << endl;
  }	
  WriteConf.close();

  return;
}

void ConfXYZ(int nconf){ //Write configuration in .xyz format
  ofstream WriteXYZ;

  WriteXYZ.open("frames/config_" + to_string(nconf) + ".xyz");
  WriteXYZ << npart << endl;
  WriteXYZ << "This is only a comment!" << endl;
  for (int i=0; i<npart; ++i){
    WriteXYZ << "LJ  " << Pbc(x[i]) << "   " <<  Pbc(y[i]) << "   " << Pbc(z[i]) << endl;
  }
  WriteXYZ.close();
}

double Pbc(double r){  //Algorithm for periodic boundary conditions with side L=box
    return r - box * rint(r/box);
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
