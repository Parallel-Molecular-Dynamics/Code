#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>

using namespace std;

int main(){

//////////////////Initialization//////////////////////////
int N = 16; //Number of Particles
double L = 1.0; //Length of Box
int iters = 10; //Number of Iterations
double final_time  = 1;
double x[iters][N];
double y[iters][N];
double vx[iters][N];
double vy[iters][N];

//number of atoms in each direction
// make it such that n**2 >N, can controll the spacing through this parameter
int n = int(ceil(pow(N,1.0/2)));

double spacing = L/n;
double dt  = final_time /iters;

// construct a trivial random generator engine from a time-based seed:
unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
std::default_random_engine generator (seed);
std::normal_distribution<double> distribution (0.0,1.0);


//initialize positions
int k =0; //counter
for (int i=0;i<n;i++){
  for (int j=0;j<n;j++){
     if (k<N){
       //initialize positions
       x[0][k] = (i + 0.5)*spacing;
       y[0][k] =  (j + 0.5)*spacing;
       //initialize velocities
       vx[0][k] = distribution(generator);
       vy[0][k] = distribution(generator);
       //cout << x[0][k] << " " << y[0][k] << " " << vx[0][k] << " " << vy[0][k] << endl;
     }
     k++;

  }

}
/////////////////////////////////////////////////////////


// Updating positions and Velocities Using Verlet method
 for (int t=1; t<iters; t++) {
      Fx, Fy = forces(x[t],y[t]);
      for (int j=0; j<N; j++) {
           x[t+1][j] = x[t][j] + vx[t][j]*dt + 0.5*Fx[j]*dt*dt;
           y[t+1][j] = y[t][j] + vy[t][j]*dt + 0.5*Fy[j]*dt*dt;
          }

      Fx1, Fy1 = forces(x[t+1],y[t+1]);

      for (int j=0; j<N; j++) {
          vx[t+1][j] = vx[t][j] + 0.5 * dt * (Fx1[j] + Fx[j]) ;
          vy[t+1][j] = vy[t][j] + 0.5 * dt * (Fy1[j] + Fy[j]) ;

      }
 } 
 

/////////////////////////Output////////////////////////
ofstream out {"molecular_data.csv"};
out<<fixed<<setprecision(4);
  for (int i = 0; i<iters; i++)
	for(int j=0; j<N; ++j){
	    out << x[i][j] << " " << y[i][j] << " " << vx[i][j] << " " << vy[i][j] << endl;
		}

out.close();



}
