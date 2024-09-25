// Numerical Integration
#include <iostream>
#include <math.h>
#include <cmath>
#include <ctime>
#include <fstream>
#include <vector>

using namespace std;

double piknown = 3.141592653589793238462643383279502884197169399375105820974944592307816406286;

double myfunc (double x)
{ double r;
  r=4./(1+x*x);
  return r;
}
// 6.1 Creating the trapezoidal approx function.
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

// 6.2 Creating the romberg approx function.
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
    // cout << "ROMBERG Theory minus calculated:     " << Romberg(&myfunc,0.0,1.0,6) - piknown << endl

    // 6.3 TIMING THE METHODS

    // Initiate vectors to place the data into so I can export it later.
    vector<double> trap_time;
    vector<double> trap_intervals;
    vector<double> trap_precision;

    vector<double> romb_time;
    vector<double> romb_intervals;
    vector<double> romb_precision;

    // First loop will work with trapezoid. 15 different 'n's, each of them calculated 50,000 times. Using <clock> to measure time, start before the loop and end after. 

    for (int i=0;i<15;i++){
        // Set the start time.
        clock_t t1 = clock();

        for (int k=1;k<50000;k++){
            Trapezoidal(&myfunc,0.0,1.0,pow(2,i));
        }
        // End the clock
        clock_t t2 = clock();
        
        // Give the list the information on how many intervals, the time it took, and the precision of the result.
        trap_intervals.push_back(pow(2,i));
        // Below line in necessary to convert cpu time to seconds, while preserving it as a double. (originally came out as a clock type? Threw errors.)
        double timeinsec = static_cast<double>(t2-t1)/CLOCKS_PER_SEC;
        trap_time.push_back(timeinsec);
        trap_precision.push_back(piknown/Trapezoidal(&myfunc,0.0,1.0,pow(2,i))*100);
    }

    // Now Romberg, same idea.

    for (int i=0;i<15;i++){
        // Set the start time.
        clock_t t1 = clock();

        for (int k=1;k<50000;k++){
            Romberg(&myfunc,0.0,1.0,i);
        }
        // End the clock
        clock_t t2 = clock();
        
        // Give the list the information on how many intervals, the time it took, and the precision of the result.
        romb_intervals.push_back(i);
        double timeinsec = static_cast<double>(t2-t1) /CLOCKS_PER_SEC;
        romb_time.push_back(timeinsec);
        romb_precision.push_back(piknown/Romberg(&myfunc,0.0,1.0,i)*100);
    }


    // CODE TO OUTPUT THIS DATA INTO A CSV.
    ofstream output("dat.csv");

    output << "type,intervals,time,precision \n";

    for (int i=0; i<trap_intervals.size();i++){
        
        double inter = trap_intervals[i];
        double timer = trap_time[i];
        double precicer = trap_precision[i];

        output << "trap," << inter << "," << timer << "," << precicer << "\n";
    }

    for (int i=0; i<romb_intervals.size();i++){

        double inter = romb_intervals[i];
        double timer = romb_time[i];
        double precicer = romb_precision[i];

        output << "romb," << inter << "," << timer << "," << precicer << "\n";
    }

    output.close();
    cout << "written." << endl;
}