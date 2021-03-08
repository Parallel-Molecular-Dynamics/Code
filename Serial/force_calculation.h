#include <iostream>
#include <stdlib.h>
#include<vector>
#include<cmath>


using namespace std;

extern const int N;
extern const double sigma;
extern const double epsilon;
extern const double cut_off;
extern const double delta;

double distance(double x1,double y1, double x2, double y2){

    double xdiff = x1 - x2;
    double ydiff = y1 - y2;

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

double LJ_pot_deriv(double r,double coeffs[]){
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

vector<vector<double>> force_calculation(double x_positions[], double y_positions[],double coeffs[]){
    double r;
    double potential_derivative;
    double change_in_x_force,change_in_y_force;
    vector<vector<double>> F(N,vector<double>(2,0));
    
    for(int i=0;i<N;++i){
        for(int j=i+1;j<N;++j){
            r = distance(x_positions[i],y_positions[i], x_positions[j], y_positions[j]);
            if(r<cut_off*sigma+delta){
                potential_derivative = LJ_pot_deriv(r,coeffs);
            } else{
                potential_derivative = 0;
            }
            change_in_x_force = (potential_derivative/r)*(x_positions[j] - x_positions[i]);
            change_in_y_force = (potential_derivative/r)*(y_positions[j] - y_positions[i]);

            F[i][0] = F[i][0] + change_in_x_force;
            F[j][0] = F[j][0] - change_in_x_force; 

            F[i][1] = F[i][1] + change_in_y_force;
            F[j][1] = F[j][1] - change_in_y_force; 
        }
    }
        
    return F;
}
