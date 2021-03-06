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
  Input();             //Inizialization
  int nconf = 1;
  //Eq_Temp();

  ofstream Term;
  Term.open("Risultati/termal.dat",ios::app);

  cout<<endl<<endl<<"Thermalization of the system from T="<< temp << endl;
   //termalizzazione del sistema con metodo di partire da una temperatura più bassa
   for(int i = 0; i<n_term;i++){
      Move();

      if(i%10 == 0){
         Termalizzazione();
         Term << stima_temp<< endl;
      }
  }
  Term.close();
  cout << "Temperature after thermalization: " << stima_temp << endl;

  cout << endl;

  cout<<"Starting simulation"<<endl;


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

  Measure_gofr(); // FUNCTION ADDED TO CALCULATE g(r)

  ConfFinal();         //Write final configuration to restart
  return 0;
}


void Input(void){ //Prepare all stuff for the simulation
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
  ReadInput >> n_term;
  ReadInput >> phase;
  cout << "Performing -- " << phase << " -- phase" << endl;
  
  if(mode != "actual" && mode != "old") {
		cerr << "Error: you should choose 'old' or 'actual' in the reading file mode..." << endl << endl;
		return;
	}
   

  cout << "The program integrates Newton equations with the Verlet method " << endl;
  cout << "Time step = " << delta << endl;
  cout << "Number of steps = " << nstep << endl << endl;
  cout << "The mode of the input choosen is the '" << mode << "' mode" << endl << endl;
  ReadInput.close();

  T_target = temp;
//Prepare array for measurements
  iv = 0; //Potential energy
  ik = 1; //Kinetic energy
  ie = 2; //Total energy
  it = 3; //Temperature
  n_props = 4; //Number of observables

  igofr = 4;
  nbins = 100;
  n_props = n_props + nbins;
  bin_size = (box/2.0)/(double)nbins;

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
   	cout << "Prepare random velocities with center of mass velocity equal to zero " << endl;
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
	
	double sumv2, fs;

	for(int i=0; i<npart; ++i){ //VELOCITY
    		vx[i] = Pbc(x[i] - dxold[i])/(2.0 * delta);    //r(t+dt)-r(t-dt)----> guarda la funzione Move();
    		vy[i] = Pbc(y[i] - dyold[i])/(2.0 * delta);
    		vz[i] = Pbc(z[i] - dzold[i])/(2.0 * delta);
	}	
	
	double t=0.;
  	for (int i=0; i<npart; ++i) 
		t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
	
	sumv2 = (t/0.5)/(double)npart;
	stima_temp = (2.0 / 3.0) * t/(double)npart;

	fs = sqrt(3 * T_target / sumv2);   // fs = velocity scale factor
     	
	for (int i=0; i<npart; ++i){
	     vx[i] *= fs;
	     vy[i] *= fs;
	     vz[i] *= fs;

     	     xold[i] = x[i] - vx[i] * delta;   //ho messo io PBC e il 2
     	     yold[i] = y[i] - vy[i] * delta;
     	     zold[i] = z[i] - vz[i] * delta;
  	}

     return;


}

void Termalizzazione(){ 
	double sumv[3] = {0.0, 0.0, 0.0};

   	for(int i=0; i<npart; ++i){ //Verlet integration scheme

     		vx[i] = Pbc(x[i] - dxold[i])/(2.0 * delta);
     		vy[i] = Pbc(y[i] - dyold[i])/(2.0 * delta);
    		vz[i] = Pbc(z[i] - dzold[i])/(2.0 * delta);
     		sumv[0] += vx[i];
     		sumv[1] += vy[i];
     		sumv[2] += vz[i];
   	}

   	for (int idim=0; idim<3; ++idim) 
		sumv[idim] /= (double)npart;
   	double sumv2 = 0.0, fs;
    	for (int i=0; i<npart; ++i){
     		vx[i] = vx[i] - sumv[0];
     		vy[i] = vy[i] - sumv[1];
     		vz[i] = vz[i] - sumv[2];

     		sumv2 += vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i];
   	}
   	sumv2 /= (double)npart;

   	double t=0.;
  	for (int i=0; i<npart; ++i) 
		t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   	temp = (2.0 / 3.0) * t/(double)npart;
   	stima_temp=temp;

   	fs = sqrt(3 * temp / sumv2);   // fs = velocity scale factor
   	for (int i=0; i<npart; ++i){
     		vx[i] *= fs;
     		vy[i] *= fs;
     		vz[i] *= fs;
   	}
	
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
  int bin;
  double v, t, vij, wij, w;
  double dx, dy, dz, dr;
  ofstream Epot, Ekin, Etot, Temp;

  Epot.open("Risultati/output_epot.dat",ios::app);
  Ekin.open("Risultati/output_ekin.dat",ios::app);
  Temp.open("Risultati/output_temp.dat",ios::app);
  Etot.open("Risultati/output_etot.dat",ios::app);

  v = 0.0; //reset observables
  t = 0.0;
  w = 0.;

  //reset the hystogram of g(r)
  for (int k=igofr; k<igofr+nbins; ++k)
	walker[k]=0.0;

//cycle over pairs of particles
  for (int i=0; i<npart-1; ++i){
    for (int j=i+1; j<npart; ++j){

     dx = Pbc( xold[i] - xold[j] ); // here I use old configurations [old = r(t)]
     dy = Pbc( yold[i] - yold[j] ); // to be compatible with EKin which uses v(t)
     dz = Pbc( zold[i] - zold[j] ); // => EPot should be computed with r(t)

     dr = dx*dx + dy*dy + dz*dz;
     dr = sqrt(dr);

     //update of the histogram of g(r)
     bin = dr/bin_size;
     walker[igofr + bin] = walker[igofr + bin] + 2;

     if(dr < rcut){
       vij = 4.0/pow(dr,12) - 4.0/pow(dr,6);
//Potential energy
       v += vij;
       w += wij;
     }
    }          
  }

//Kinetic energy
  for (int i=0; i<npart; ++i) t += 0.5 * (vx[i]*vx[i] + vy[i]*vy[i] + vz[i]*vz[i]);
   
    stima_pot = v/(double)npart; 	//Potential energy per particle
    stima_kin = t/(double)npart; 	//Kinetic energy per particle
    stima_temp = (2.0 / 3.0) * t/(double)npart; //Temperature
    stima_etot = (t+v)/(double)npart; 	//Total energy per particle
 
    for(int i=igofr; i<igofr+nbins; i++) {	//g(r)
// evaluate Delta V(r)
		double r = (double)(i-igofr) * bin_size;
		double dV = 4.*M_PI/3. * ( pow(r+bin_size, 3) - pow(r, 3) );
		
		stima_gofr = walker[i] / (rho * dV * (double)npart);
		glob_av[i] += stima_gofr;
		glob_av2[i] += stima_gofr * stima_gofr;
	}


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

void Measure_gofr() {
	ofstream Gave;
	Gave.open("gofr.ave");
	int n = nstep/10;

	for(int i=igofr; i<igofr+nbins; i++) {
		double err_gdir = Error(glob_av[i], glob_av2[i], n);
		double r = (double)(i-igofr) * bin_size;

// print (r g(r) for each block and final from last block with error)
		Gave << r << "    " << glob_av[i]/(double)n << "    " << err_gdir << endl;
	}
	
	Gave.close();
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
double Error(double sum, double sum2, int iblk)
{
    if( iblk == 1 ) return 0.0;
    else return sqrt((sum2/(double)iblk - pow(sum/(double)iblk,2))/(double)(iblk-1));
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
