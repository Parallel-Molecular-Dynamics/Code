#include <iostream>
#include<stdlib.h>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include<vector>
#include "force_calculation.h"

using namespace std;

const int N = 2; //Number of Particles
const double sigma = 1;
const double epsilon= 1;
const double cut_off = 3;
const double delta = 0.1;

int main(){
//////////////////Initialization//////////////////////////
double L = 1.0; //Length of Box
int iters = 1000; //Number of Iterations
double spacing = 1.3;



double final_time  = 1;
double x[iters][N];
double y[iters][N];
double vx[iters][N];
double vy[iters][N];
double x_positions[N];
double y_positions[N];
vector<vector<double>>F;
vector<vector<double>>F1;



//number of atoms in each direction
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


// Updating positions and Velocities Using Verlet method
 for (int t=0; t<iters; t++) {
      F = force_calculation(x_positions,y_positions,polynomial_coeffs);
      if(t == 1){      cout<<F[0][0]<<' ' << F[0][1]<<endl;cout<<F[1][0]<<' ' << F[1][1]<<endl;}
      for (int j=0; j<N; j++) {
           x_positions[j] = x[t+1][j] = x[t][j] + vx[t][j]*dt + 0.5*F[j][0]*dt*dt;
           y_positions[j] = y[t+1][j] = y[t][j] + vy[t][j]*dt + 0.5*F[j][1]*dt*dt;
          }
      if(t == 1){cout<<x_positions[0]<<' ' <<y_positions[0]<<endl;cout<<x_positions[1]<<' ' <<y_positions[1]<<endl;}

      F1 = force_calculation(x_positions,y_positions,polynomial_coeffs);
      if(t == 1){      cout<<F1[0][0]<<' ' << F1[0][1]<<endl;cout<<F1[1][0]<<' ' << F1[1][1]<<endl;}

      for (int j=0; j<N; j++) {
          vx[t+1][j] = vx[t][j] + 0.5 * dt * (F1[j][0] + F[j][0]) ;
          vy[t+1][j] = vy[t][j] + 0.5 * dt * (F1[j][1] + F[j][1]) ;

      }
 }


/////////////////////////Output////////////////////////
ofstream out {"molecular_data.csv"};
out<<fixed<<setprecision(4);
out << "x" <<" "<< "y" <<" "<< "vx" <<" "<<"vy" << endl;
  for (int i = 0; i<iters; i++)
	for(int j=0; j<N; ++j){
	    out << x[i][j] << " " << y[i][j] << " " << vx[i][j] << " " << vy[i][j] << endl;
		}

out.close();



}
