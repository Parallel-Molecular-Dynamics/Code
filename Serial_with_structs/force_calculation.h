#include <iostream>
#include <stdlib.h>
#include<cmath>


using namespace std;

extern const int N;
extern const double sigma;
extern const double epsilon;
extern const double cut_off;
extern const double delta;

struct Particle

{
    double x;
    double y;
    double vx;
    double vy;

    double K;
    double W;
};

struct Force

{
    double x;
    double y;
};

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
    double potential;

    if(r<rc){
        potential = 4*epsilon*(power(sigma/r,12) - power(sigma/r,6));
    }
    else if(r<rc+delta){
        potential = coeffs[0] + coeffs[1]*r + coeffs[2]*power(r,2)+ coeffs[3]*power(r,3);
    }
    return potential;
}

int force_calculation(struct Force F[], struct Particle particles[],double coeffs[]){

    double r;
    double potential_derivative;
    struct Force Delta_F;

    for(int i = 0; i<N;i++){
        F[i].x = F[i].y = particles[i].W = 0;
    }

    for(int i=0;i<N;++i){
        for(int j=i+1;j<N;++j){

            r = distance(particles[i], particles[j]);

            if(r<cut_off*sigma+delta){
                potential_derivative = LJ_potential_derivative(r,coeffs);
            } else{
                potential_derivative = 0;
            }

            Delta_F.x = (potential_derivative/r)*(particles[j].x - particles[i].x);
            Delta_F.y = (potential_derivative/r)*(particles[j].y - particles[i].y);

            F[i].x +=  Delta_F.x;
            F[j].x -=  Delta_F.x; 

            F[i].y +=   Delta_F.y;
            F[j].y -=   Delta_F.y;

            particles[i].W += 0.5 * LJ_potential(r,coeffs);
            particles[j].W += 0.5 * LJ_potential(r,coeffs);
        }
    }
    return 0;
}