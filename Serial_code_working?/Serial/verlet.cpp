#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include<vector>
#include "force_calculation.h"

using namespace std;


const int N = 100; //Number of Particles
const int iters = 10000; //Number of Iterations
const double final_time  = 10;


const double sigma = 1;
const double epsilon= 1;
const double cut_off = 3;
const double delta = 0.1;
const double L = 1.0; //Length of Box
const double spacing = 1.3;


double x[iters][N];
double y[iters][N];
double vx[iters][N];
double vy[iters][N];
double x_positions[N];
double y_positions[N];

int main(){

//////////////////Initialization//////////////////////////
vector<vector<double>>F;
vector<vector<double>>F1;

// make it such that n**2 >N, can controll the spacing through this parameter
int n = int(ceil(pow(N,1.0/2)));

double dt  = final_time /iters;

// construct a trivial random generator engine from a time-based seed:
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
std::normal_distribution<double> distribution (0.0,1.0);

// Coefficients for smooth cut-off
double *polynomial_coeffs;
polynomial_coeffs = determine_polynomial_coeffs(sigma,epsilon,cut_off,delta);


/////////////////////////Output////////////////////////
ofstream out {"parameters.csv"};
out<<fixed<<setprecision(4);
out << "N" <<" " <<"iters" <<" "<< "sigma" <<" "<< "epsilon" <<" "<<"cut_off" <<" "<<"delta"<<" "<< "a" <<" "<< "b" <<" "<< "c" << " "<< "d" << endl;
out << N << " " << iters <<" " <<sigma << " " << epsilon << " " << cut_off << " " <<delta << " " << polynomial_coeffs[0]  << " " << polynomial_coeffs[1] << " " << polynomial_coeffs[2] << " " << polynomial_coeffs[3]  << endl;

out.close();


if (spacing <= pow(2,1/6)*sigma ){
   cout << "Warning: Particles are initialized within the repulsive range" <<
   endl;
}
cout << "Beginning initialisation..."<<endl;
//initialize positions
int k =0; //counter
for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
     if (k<N){
       //initialize positions
       x_positions[k] = x[0][k] = (i + 0.5)*spacing;
       y_positions[k] = y[0][k] =  (j + 0.5)*spacing;
       //initialize velocities
       vx[0][k] = 0;//distribution(generator);
       vy[0][k] = 0;//distribution(generator);
       //cout << x[0][k] << " " << y[0][k] << " " << vx[0][k] << " " << vy[0][k] << endl;
     }
     k++;

  }

}
/////////////////////////////////////////////////////////
cout<< "...done." << endl;

// Updating positions and Velocities Using Verlet method
F = force_calculation(x_positions,y_positions,polynomial_coeffs);
 for (int t=0; t<iters-1; t++) {

      for (int j=0; j<N; j++) {
           x_positions[j] = x[t+1][j] = x[t][j] + vx[t][j]*dt + 0.5*F[j][0]*dt*dt;
           y_positions[j] = y[t+1][j] = y[t][j] + vy[t][j]*dt + 0.5*F[j][1]*dt*dt;
          }

      F1 = force_calculation(x_positions,y_positions,polynomial_coeffs);

      for (int j=0; j<N; j++) {
          vx[t+1][j] = vx[t][j] + 0.5 * dt * (F1[j][0] + F[j][0]) ;
          vy[t+1][j] = vy[t][j] + 0.5 * dt * (F1[j][1] + F[j][1]) ;

      }
      F = F1;
   cout<< "Timestep "<< t+1 << " of " << iters << " done."<<endl;   
 }


/////////////////////////Output////////////////////////
cout<< "Beginning output..." << endl;
ofstream mout {"molecular_data.csv"};
mout<<fixed<<setprecision(4);
mout << "x" <<" "<< "y" <<" "<< "vx" <<" "<<"vy" << endl;
  for (int i = 0; i<iters; i++)
	for(int j=0; j<N; ++j){
	    mout << x[i][j] << " " << y[i][j] << " " << vx[i][j] << " " << vy[i][j] << endl;
		}

mout.close();
cout<< "...done.";


}
