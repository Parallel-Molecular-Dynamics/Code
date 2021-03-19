#include <iostream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <random>
#include <mpi.h>
#include "force_calculation.h"

using namespace std;


const int N = 200;
const int nproc = 2;
const int iters = 10000; 
const double final_time  = 10;


const double sigma = 1;
const double epsilon= 1;
const double cut_off = 3;
const double delta = 0.1;
const double L = 1.0; 
const double spacing = 1.3;
const double dt  = final_time /iters;
const int Np = N/nproc;

struct Particle all_particles[N];
struct Force F[N];
struct Particle particle_subset[Np];

int main(){

    MPI_Init(NULL,NULL);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    if(N%nproc != 0){
        cout<<"N/nproc must leave no remainder"<<endl;
        return 1;
    }

    MPI_Datatype particle_type = make_MPI_struct_particle();
    MPI_Datatype force_type = make_MPI_struct_force();

    int idx;
    double K=0;
    double r;
    double W=0;
    
    if (spacing <= pow(2,1/6)*sigma ){
        cout << "Warning: Particles are initialized within the repulsive range" <<
        endl;
    }

    // construct a trivial random generator engine from a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);

    // Coefficients for smooth cut-off in the LJ potential.
    double *polynomial_coeffs;
    polynomial_coeffs = determine_polynomial_coeffs(sigma,epsilon,cut_off,delta);

    /////////////////////////Output Parameters////////////////////////
    if(rank == 0){
        ofstream out {"parameters.csv"};
        out<<fixed<<setprecision(4);

        out << "N" <<"," <<"iters" <<","<< "sigma" <<","<< "epsilon" <<","<<"cut_off" <<","<<"delta"<<","<< "a" <<","<< "b" <<","<< "c" << ","<< "d" << endl;
        out << N << "," << iters <<"," <<sigma << "," << epsilon << "," << cut_off << "," <<delta << "," << polynomial_coeffs[0]  << "," << polynomial_coeffs[1] << "," << polynomial_coeffs[2] << "," << polynomial_coeffs[3]  << endl;

        out.close();
    }


    // make it such that n**2 >N, can control the spacing through this parameter
    int n = int(ceil(pow(N,1.0/2)));

    if(rank==0){cout << "Beginning initialisation..."<<endl;}
    int k =0; 
    for (int i=0;i<n;i++){
        for (int j=0;j<n;j++){
            if (k<N){
            
                all_particles[k].x = (i + 0.5)*spacing;
                all_particles[k].y = (j + 0.5)*spacing;
              
                all_particles[k].vx = 0;//distribution(generator);
                all_particles[k].vy = 0;//distribution(generator);
                //if(rank == 0){cout << all_particles[k].x << " " << all_particles[k].y << " " << all_particles[k].vx << " " << all_particles[k].vy << endl;}
            }
            k++;
        }
    }
    if(rank == 0){cout<< "...done." << endl;}

    // Pick out the relevant subset of particles for the given processor.
    for(int i=0; i <Np;i++){
        idx =(Np*rank)+i;
        particle_subset[i].x = all_particles[idx].x;
        particle_subset[i].y = all_particles[idx].y;
        particle_subset[i].vx = all_particles[idx].vx;
        particle_subset[i].vy = all_particles[idx].vy;

    }

    force_calculation(rank,F,all_particles,polynomial_coeffs);

    ////Open an output file and read out initial conditions/////
    ofstream mout {"molecular_data.csv"};
    ofstream Eout {"energy_data.csv"};
    if(rank == 0){
        mout<< "x,y,vx,vy"<<endl;

        for(int i = 0; i<N; i++){
            mout << all_particles[i].x << "," << all_particles[i].y<<","<<all_particles[i].vx<<","<<all_particles[i].vy<<endl;

            K += 0.5*(pow(all_particles[i].vx,2)+pow(all_particles[i].vy,2));
            for(int j=i+1;j<N;j++){
                r = distance(all_particles[i],all_particles[j]);
                W +=LJ_potential(r,polynomial_coeffs);
            }
        }

        Eout << "K,W"<<endl;
        Eout << K << "," << W << endl; 
    }

     for (int t=0; t<iters-1; t++) {
        W = 0;
        for (int i=0; i<Np; i++) {
            idx =(Np*rank)+i;

            particle_subset[i].vx += 0.5*dt*F[idx].x;
            particle_subset[i].vy += 0.5*dt*F[idx].y;

            particle_subset[i].x += dt*particle_subset[i].vx;
            particle_subset[i].y += dt*particle_subset[i].vy;
        }

        MPI_Allgather(&particle_subset,Np,particle_type, &all_particles,Np,particle_type,MPI_COMM_WORLD);


        force_calculation(rank,F,all_particles,polynomial_coeffs);

        for (int i=0; i<Np; i++) {
            idx =(Np*rank)+i;
            particle_subset[i].vx += 0.5*dt*F[idx].x;
            particle_subset[i].vy += 0.5*dt*F[idx].y;
        }

        MPI_Allgather(&particle_subset,Np,particle_type,&all_particles,Np,particle_type,MPI_COMM_WORLD); 
        if(rank == 0){cout<< "Timestep "<< t+1 << " of " << iters-1 << " done."<<endl;}

        if(rank == 0){
            K = 0;
            W = 0;
            for(int i = 0; i<N; i++){
                mout << all_particles[i].x << "," << all_particles[i].y<<","<<all_particles[i].vx<<","<<all_particles[i].vy<<endl;

                K += 0.5*(pow(all_particles[i].vx,2)+pow(all_particles[i].vy,2));
                for(int j=i+1;j<N;j++){
                    r = distance(all_particles[i],all_particles[j]);
                    W +=LJ_potential(r,polynomial_coeffs);
                }
            }
            Eout << K << "," << W << endl;  
        }
        
    }
    mout.close();
    Eout.close(); 

    MPI_Finalize();
    return 0;
}

