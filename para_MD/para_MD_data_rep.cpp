#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <mpi.h>
#include "parallel_force_calculation.h"

using namespace std;

const int nproc          = 2;
const int    iters       = 30000; //Number of Iterations
const double final_time  = 300;
const double sigma       = 1;
const double epsilon     = 1;
const int    n           = 2;
const int    N           = n*n; //Number of Particles
const double spacing     = 1.3; //1.1;
const double L           = n*spacing;
const double cut_off     = L/2;//0.95 * L/2;
const double delta       = 0.1;
const double dt          = final_time /iters;
const int Np             = N/nproc;


struct Particle particles[N];
struct Force F[N];
struct Particle particle_subset[Np];

int main(){

    MPI_Init(NULL,NULL);

    double *polynomial_coeffs;
    int rank,idx;

    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    // Warnings for non-useful parameter choices.
    if(N%nproc != 0){
        cout<<"N/nproc must leave no remainder"<<endl;
        return 1;
    }

    if (spacing <= pow(2,1/6)*sigma ){
        cout << "Warning: Particles are initialized within the repulsive range" <<
        endl;
    }

    MPI_Datatype particle_type = make_MPI_struct_particle();

    // Coefficients for smooth cut-off in the LJ potential.
    polynomial_coeffs = determine_polynomial_coeffs(sigma,epsilon,cut_off,delta);

    /////////////////////////Output Parameters////////////////////////
    if(rank == 0){output_parameters(iters, polynomial_coeffs);}

    initialise_particles(rank,particles,particle_subset);
    
    force_calculation(rank,F,particles,polynomial_coeffs);

  ////Open an output file and read out initial conditions/////
    ofstream mout {"molecular_data.csv"};
    ofstream Eout {"energy_data.csv"};
    if(rank == 0){
        mout<< "x,y,vx,vy"<<endl;
        Eout << "K,W"<<endl;

        calculate_energy_and_output(particles,polynomial_coeffs,Eout,mout);
    }

     for (int t=0; t<iters-1; t++) {
        for (int i=0; i<Np; i++) {
            idx =(Np*rank)+i;

            particle_subset[i].vx += 0.5*dt*F[idx].x;
            particle_subset[i].vy += 0.5*dt*F[idx].y;

            particle_subset[i].x += dt*particle_subset[i].vx;
            particle_subset[i].x = modulo(particle_subset[i].x,L);

            particle_subset[i].y += dt*particle_subset[i].vy;
            particle_subset[i].y = modulo(particle_subset[i].y,L);
        }

        MPI_Allgather(&particle_subset,Np,particle_type, &particles,Np,particle_type,MPI_COMM_WORLD);

        force_calculation(rank,F,particles,polynomial_coeffs);

        for (int i=0; i<Np; i++) {
            idx =(Np*rank)+i;
            particle_subset[i].vx += 0.5*dt*F[idx].x;
            particle_subset[i].vy += 0.5*dt*F[idx].y;
        }

        MPI_Allgather(&particle_subset,Np,particle_type,&particles,Np,particle_type,MPI_COMM_WORLD); 

        if(rank == 0){
            cout<< "Timestep "<< t+1 << " of " << iters-1 << " done."<<endl;
            calculate_energy_and_output(particles,polynomial_coeffs,Eout,mout);
        }
        
    }
    mout.close();
    Eout.close();    

    MPI_Finalize(); 
    return 0;
}

