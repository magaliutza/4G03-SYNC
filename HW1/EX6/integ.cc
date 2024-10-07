// Numerical Integration
#include <iostream>
#include <math.h>
#include <cmath>
using namespace std;

double piknown = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

double myfunc (double x)
{ double r;
  r=4./(1+x*x);
  return r;
}

double Trapezoidal(double (*f)(double x), double xl,double xh,int n){ 

    // Segment, find the del, initiate the sum for the integral.
    double len = xh - xl;
    //cout << "len: " << len << endl;
    double delx = len/n;
    //cout << "delx: " << delx << endl;
    double summa; 
    
    for (int i = 0; i<n; i++){
        double a = xl + delx*i;
        double b = a + delx;
        summa += (b-a)*0.5*(f(a)+f(b));
    }

    return summa;
}

double Romberg(double (*f)(double x), double xl,double xh,int N){ 
    
    // Start with intializing romburg R(0,0) = (f(xl) + f(xh))/2
    double Roo = (f(xh)+f(xl))/2;

    // If someone were to call upon R(0,0), let's not waste their time.
    if (N==0){
        return Roo;
    }

    // Initiate the matrix to store all possible values of R(n,m) for 0<= n,m <=N. We will use them inductively.
    
    double r_n_m[N+1][N+1];
    r_n_m[0][0] = Roo;

    // Start the loop to achieve R(n,0) first. We want to first find R(1,0) so we'll set n = 1 to begin and move up.
    for (int n = 1; n<=N; n++){
        
        // Define the h_n we'll used, based on the step size and the n we're on.
        double h = (xh - xl)/pow(2,n);

        // Initiate a sum and loop for the individual sum to get the next R(n,0).
        double summa = 0;
        for (int k = 1; k<= pow(2,(n-1)); k++){

            double val = xl + (2*k-1)*h;
            summa += f(val);

        }
        
        r_n_m[n][0] = r_n_m[n-1][0]/2 + h*summa;
        
    }
    // We now should have all R(n,0) available.
    // Then move on to induct the rest of R(n,m) recursively.

    for (int m = 1; m<=N; m++){

        for (int n = m; n<=N; n++){

            r_n_m[n][m] = r_n_m[n][m-1] + 1/(pow(4,m)-1) * (r_n_m[n][m-1] - r_n_m[n-1][m-1]);

        }


    }

    return r_n_m[N][N];
}

int main ()
{
    // cout.precision(24);
    // cout << "The Trapezoidal result is     " <<  Trapezoidal(&myfunc,0.0,1.0,64) << endl;
    // cout << "TRAPEZOID Theory minus calculated:     " << Trapezoidal(&myfunc,0.0,1.0,64) - piknown << endl;
    // cout << "The Romberg result is         " << Romberg(&myfunc,0.0,1.0,6) << endl;\
    // cout << "ROMBERG Theory minus calculated:     " << Romberg(&myfunc,0.0,1.0,6) - piknown << endl;


    // cout << "Romberg 6,6:   " << Romberg(&myfunc,0.0,1.0,6) << endl;
    // cout << "Romberg 8,8:   " << Romberg(&myfunc,0.0,1.0,8) << endl;
    // cout << "Romberg 10,10: " << Romberg(&myfunc,0.0,1.0,10) << endl;

    

}
