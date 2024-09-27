#include <iostream>
#include <cmath>
#include <ctime>
using namespace std;

double piknown =  3.141592653589793238462643383279502884197169399375105820974944592307816406286;
double sqrtpiknown = 1.7724538509055160272981674833411451827975494561223871282138077898;

// Making a factorial using the gamma function +1.
double fact(double a){
	return tgamma(a+1);
}

double sqrt_taylor(double x, int n = 50){
	// initiate the parameters.
	double summa = 0;
	double coeff = 1;
	double exponent = 0.5;

	for (int i = 0; i < n; i++){
		
		// Add each term of the taylor series, starting with 1. Since we set x-1, this is centered around 1.
		summa += coeff * pow(x-1,i) / fact(i);
		
		// Update the coeff and exponent.
		coeff = coeff*exponent;
		exponent-=1;

	}
	return summa;
}

double sqrt_newt(double x, double tolerance = 1e-24){
	
	// Initial guess, apparently x/2 is a good start.
	double step = x/2;
	double next;
	// Set a max number of steps to prevent runaway.
	for (int i = 0; i<1e15; i++){
		
		// Newtons method...
		next = 0.5 * (step + x / step);

		// Check the tolerance, end if its too small.
		if (next - step < tolerance){
			//cout << "Broken." << endl;
			break;
		}
		// Iterate if necessary.
		step = next;
	}
	return step;
}

double sqrt_frac(double x, int n){

	// We'll hard code the first approximation to be 2, since were looking for 1.7~

	double step = 2;

	// Return when at the bottom of the fraction.
	if (n==0)
		return step + (x - step*step) / (step);
	
	// Otherwise, find the next fraction.	
	else {
		//cout << "n: " << n << endl;
		return step + (x - step*step) / (step + sqrt_frac(x,n-1)); 
	}

}


int main(){
	int n = 13;	

	// cout << "Taylor approx: " << sqrt_taylor(piknown-1) << endl;
	// cout <<  "Taylor precision: " << sqrt_taylor(piknown-1) - sqrtpiknown << endl;

	// cout << "Newton approx: " << sqrt_newt(piknown) << endl;
	// cout <<  "Newton precision: " << sqrt_newt(piknown) - sqrtpiknown << endl;

	// cout << "Fraction approx: " << sqrt_frac(piknown,n) << endl;
	// cout <<  "Fraction precision: " << sqrt_frac(piknown,n) - sqrtpiknown << endl;


	int iterations = 500000;
	// TAYLOR TIMING.

    clock_t t1 = clock();

    for (int k=1;k<iterations;k++){
        sqrt_taylor(piknown-1);
    }
    // End the clock
    clock_t t2 = clock();
    
    // Below line in necessary to convert cpu time to seconds, while preserving it as a double. (originally came out as a clock type? Threw errors.)
    double timeinsec = static_cast<double>(t2-t1)/CLOCKS_PER_SEC;
    cout << "taylor time: " << timeinsec << endl;
    cout << "taylor precision: " << sqrtpiknown-sqrt_taylor(piknown-1) << endl;
    

    // NEWT TIMING.

    clock_t t3 = clock();

    for (int k=1;k<iterations;k++){
        sqrt_newt(piknown);
    }
    // End the clock
    clock_t t4 = clock();
    
    // Below line in necessary to convert cpu time to seconds, while preserving it as a double. (originally came out as a clock type? Threw errors.)
    timeinsec = static_cast<double>(t4-t3)/CLOCKS_PER_SEC;
    cout << "newt time: " << timeinsec << endl;
    cout << "newt precision: " << sqrtpiknown-sqrt_newt(piknown)<< endl;

    // FRAC TIMING.

    clock_t t5 = clock();

    for (int k=1;k<iterations;k++){
        sqrt_frac(piknown,n);
    }
    // End the clock
    clock_t t6 = clock();
    
    // Below line in necessary to convert cpu time to seconds, while preserving it as a double. (originally came out as a clock type? Threw errors.)
    timeinsec = static_cast<double>(t6-t5)/CLOCKS_PER_SEC;
    cout << "frac time: " << timeinsec << endl;
    cout << "frac precision: " << sqrtpiknown-sqrt_frac(piknown,n)<< endl;

}