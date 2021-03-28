#include <iostream>
#include <stdlib.h>
#include<cmath>
#include <chrono>
#include <random>
#include<mpi.h>


using namespace std;

extern const int N;
extern const int Np; 
extern const double sigma;
extern const double epsilon;
extern const double cut_off;
extern const double delta;
extern const int n;
extern const double spacing;
extern const double L;

struct Particle

{
    double x;
    double y;
    double vx;
    double vy;
};

struct Force

{
    double x;
    double y;
};

MPI_Datatype make_MPI_struct_particle(){
    MPI_Datatype particle_type;
    int lengths[4] = { 1,1,1,1};
    const MPI_Aint displacements[4] = { 0, sizeof(double), 2*sizeof(double),3*sizeof(double)};
    MPI_Datatype types[4] = { MPI_DOUBLE, MPI_DOUBLE, MPI_DOUBLE,MPI_DOUBLE};
    MPI_Type_create_struct(4, lengths, displacements, types, &particle_type);
    MPI_Type_commit(&particle_type);

    return particle_type;
}

double distance(struct Particle p1, struct Particle p2){

    double xdiff = p1.x - p2.x;
    double ydiff = p1.y - p2.y;

return sqrt(xdiff*xdiff + ydiff*ydiff);
}

double power(double x,int n){

    double xn = 1;
    for(int i = 0;i<n;i++){
        xn = xn* x;
    }
    
    return xn;
}

double * determine_polynomial_coeffs(double sigma,double epsilon,int cut_off,double delta){
    int n =4;

    double A[n][n+1];
    static double coeffs[4];
    double rc = cut_off*sigma;
    
    A[0][0] = 1;
    A[0][1] = rc;
    A[0][2] = power(rc,2);
    A[0][3] = power(rc,3);
    A[0][4] = 4*epsilon*(power(sigma/rc,12)-power(sigma/rc,6));

    A[1][0] = 0;
    A[1][1] = 1;
    A[1][2] = 2*rc;
    A[1][3] = 3*power(rc,2);
    A[1][4] = 24*epsilon*((power(sigma,6)/power(rc,7))-(2*power(sigma,12)/power(rc,13)));

    A[2][0] = 1;
    A[2][1] = rc + delta;
    A[2][2] = power(rc+delta,2);
    A[2][3] = power(rc+delta,3);
    A[2][4] = 0;

    A[3][0] = 0;
    A[3][1] = 1;
    A[3][2] = 2*(rc+delta);
    A[3][3] = 3*power(rc+delta,2);
    A[3][4] = 0;

    for(int i=0;i<n;i++){                   
        for(int j=i+1;j<n;j++){
            if(abs(A[i][i]) < abs(A[j][i])){
                for(int k=0;k<n+1;k++){

                    A[i][k]=A[i][k]+A[j][k];
                    A[j][k]=A[i][k]-A[j][k];
                    A[i][k]=A[i][k]-A[j][k];

                }
            }
        }
    }

     /* performing Gaussian elimination */
    for(int i=0;i<n-1;i++){
        for(int j=i+1;j<n;j++){

            double f=A[j][i]/A[i][i];

            for(int k=0;k<n+1;k++){

                A[j][k]=A[j][k]-f*A[i][k];
            }
        }
    }

    /* Backward substitution for discovering values of unknowns */
    for(int i=n-1;i>=0;i--){                     
        coeffs[i]=A[i][n];
        for(int j=i+1;j<n;j++){
            if(i!=j){
              coeffs[i]=coeffs[i]-A[i][j]*coeffs[j];
            }          
        }
        coeffs[i]=coeffs[i]/A[i][i];  
    }
    
    return coeffs;
}

double LJ_potential_derivative(double r,double coeffs[]){

    double rc = cut_off*sigma;
    double potential_derivative;

    if(r<rc){
        potential_derivative = 24*epsilon*((power(sigma,6)/power(r,7))-(2*power(sigma,12)/power(r,13)));
    }
    else if(r<rc+delta){
        potential_derivative = coeffs[1] + 2*coeffs[2]*r + 3 *coeffs[3]*power(r,2);
    }
    return potential_derivative;
}

double LJ_potential(double r,double coeffs[]){

    double rc = cut_off*sigma;
    double potential=0;

    if(r<rc){
        potential = 4*epsilon*(power(sigma/r,12) - power(sigma/r,6));
    }
    else if(r<rc+delta){
        potential = coeffs[0] + coeffs[1]*r + coeffs[2]*power(r,2)+ coeffs[3]*power(r,3);
    }
    return potential;
}

int force_calculation(int rank, struct Force F[], struct Particle particles[],double coeffs[]){

    double r;
    double potential_derivative;
    double dx;
    double dy;
    struct Force Delta_F;
    int idx;

    for(int i = 0; i<N;i++){
        F[i].x = F[i].y  = 0;
    }

    int break_loop = 0;
    
    for(int i=0;i<Np;++i){
        idx = (Np*rank)+i;
        for(int j=0;j<Np*rank;++j){ 

            dx  = particles[j].x - particles[idx].x;
            dy  = particles[j].y - particles[idx].y;
		    
		    if (dx > L/2){
                dx -= L;
		    }
		    else if (dx < -L/2){
                dx += L;    
		    }

		    if (dy > L/2){
                dy -= L;
		    }
		    else if (dy < -L/2){
                dy += L;    
	 	    }
		    
            r = sqrt(dx*dx + dy*dy);
       
            if(r<cut_off*sigma+delta){
                potential_derivative = LJ_potential_derivative(r,coeffs);
                
                Delta_F.x = (potential_derivative/r) * dx;
                Delta_F.y = (potential_derivative/r) * dy;

                F[idx].x +=  Delta_F.x;
                F[idx].y +=   Delta_F.y;

                F[j].x -=  Delta_F.x;
                F[j].y -=   Delta_F.y;
            }

        }

        for(int j=idx+1;j<N;++j){ 

            dx  = particles[j].x - particles[idx].x;
            dy  = particles[j].y - particles[i].y;
		    
		    if (dx > L/2){
                dx -= L;
		    }
		    else if (dx < -L/2){
                dx += L;    
		    }

		    if (dy > L/2){
                dy -= L;
		    }
		    else if (dy < -L/2){
                dy += L;    
	 	    }
		    
            r = sqrt(dx*dx + dy*dy);
       
            if(r<cut_off*sigma+delta){
                potential_derivative = LJ_potential_derivative(r,coeffs);

                Delta_F.x = (potential_derivative/r) * dx;
                Delta_F.y = (potential_derivative/r) * dy;

                F[idx].x +=  Delta_F.x;
                F[idx].y +=   Delta_F.y;

                F[j].x -=  Delta_F.x;
                F[j].y -=   Delta_F.y;
            }
        }
    }

    return 0;
}

void output_parameters(int iters,double polynomial_coeffs[]){
    ofstream out {"parameters.csv"};
    out<<fixed<<setprecision(4);

    out << "N" <<"," <<"iters" <<","<< "sigma" <<","<< "epsilon" <<","<<"cut_off" <<","<<"delta"<<","<< "a" <<","<< "b" <<","<< "c" << ","<< "d" << ',' <<"L" << endl;
    out << N << "," << iters <<"," <<sigma << "," << epsilon << "," << cut_off << "," <<delta << "," << polynomial_coeffs[0]  << "," << polynomial_coeffs[1] << "," << polynomial_coeffs[2] << "," << polynomial_coeffs[3] <<','<< L  << endl;

    out.close();
}

void calculate_energy_and_output(struct Particle particles[],double coeffs[],ofstream &Energy_out,ofstream &MD_out){

    double W = 0;
    double K = 0;
    double r,dx,dy;

    for(int i = 0; i<N; i++){
        MD_out << particles[i].x << "," << particles[i].y<<","<<particles[i].vx<<","<<particles[i].vy<<endl;

        K += 0.5*(pow(particles[i].vx,2)+pow(particles[i].vy,2));

        for(int j=i+1;j<N;j++){
            dx  = particles[j].x - particles[i].x;
            dy  = particles[j].y - particles[i].y;
		    
		    if (dx > L/2){
                dx -= L;
		    }
		    else if (dx < -L/2){
                dx += L;    
		    }

		    if (dy > L/2){
                dy -= L;
		    }
		    else if (dy < -L/2){
                dy += L;    
	 	    }
		    
            r = sqrt(dx*dx + dy*dy);
            W +=LJ_potential(r,coeffs);
        }
    }
    Energy_out << K << "," << W << endl;
}

int initialise_particles(int rank,struct Particle all_particles[], struct Particle particle_subset[]){
    // construct a trivial random generator engine from a time-based seed:
    unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
    std::default_random_engine generator (seed);
    std::normal_distribution<double> distribution (0.0,1.0);
    
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
    // Pick out the relevant subset of particles for the given processor.
    int idx;
    for(int i=0; i <Np;i++){
        idx =(Np*rank)+i;
        particle_subset[i].x = all_particles[idx].x;
        particle_subset[i].y = all_particles[idx].y;
        particle_subset[i].vx = all_particles[idx].vx;
        particle_subset[i].vy = all_particles[idx].vy;
    }

    return 0;
}

double modulo(double n, double d){

    double mod = fmod(n,d);

    if (mod < 0){
       mod = d+mod;
    }
return mod;
}
