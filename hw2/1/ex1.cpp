#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

// Define our f(x) which we'd like to integrate. For 1.1, we use exp(-x). For 1.2, we use exp(-x*x)
double f_x(double x){
	return exp(-x*x);
}

// Write the definite integral f(x) from 0-1 to create the normalizing weighting func. This is for 1.1. 
double weight = (exp(1) - 1) / exp(1);

// The weighting func itself, g(x). We use g(x) = exp(-x) since it is 'similar' to the integral we actually want to estimate. 
double g_x(double x){
	return exp(-x) / weight;
}

// The integral of the weigh func, inversed, G^-1(x), to apply a normal distribution on [0,1]. 
double G_inv(double u){
	return -log(1 - (u * (weight)));
}


// Begin the actual calculations.
int main(){


	// We define the known integral results to compare to later. Int1 is exp(-x) on [0,1], Int2 is exp(-x*x) on [0,10].
	double int1known = 0.6321205588285576784044762298385391325541888689682321654;
	double int2known = 0.9999546000702375151484644;
	
	// This is for 1.2 when I want to replace the weighting.
	weight = int2known;
	
	// RNG related code. Initiating a random device to get a random seed, which is used for the twister engine.. 
	random_device rand;
	seed_seq seed{rand(), rand(), rand(), rand(), rand(), rand(), rand()};
	mt19937 eng(seed);
	uniform_real_distribution<> dist(0,1);
	

	// Initialize number of steps, and the sum at 0.
	int n = 5;
	double summa = 0;

	// Looping number of steps...
	for (int i = 0; i < n; i++){

		// u is randomly generated on [0,1]
		double u = dist(eng);
		cout << u << endl;
		// We get an value of x from G^-1(u) = x
		double xval = G_inv(u);
		// Add to the sum the estimated point.
		summa += f_x(xval)/g_x(xval);
		cout << "summa:" << summa << endl;
	}

	// Find the result by dividing the sum by the number of steps, basically averaging. 
	double integ = summa/n;
	
	cout << setprecision(10) << weight << endl;
	cout << "Integral estiamte: " << integ << endl;
	cout << "Integral precision: " << abs(integ - int1known) << endl;

	ofstream file;
	file.open("dat.csv");
	file << "n, error\n";

	for (int i = 0; i < 6; i++){

		n = pow(10,i);
		summa = 0;

		for (int k = 0; k < n; k++){
			double u = dist(eng);
			double xval = G_inv(u);
			summa += f_x(xval)/g_x(xval);
		}

		integ = summa/n;
		double integ_calc = 0.88622692545275801364908374167057259139877;
		

		cout << "result: " << integ << endl;
		cout << "steps (n): " << n << endl;
		cout << "error: " << abs(integ-integ_calc) << endl;

		file << n << ", " << abs(integ-integ_calc) << endl;
	}

	file.close();


}
