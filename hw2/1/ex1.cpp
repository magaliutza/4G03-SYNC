#include <iostream>
#include <random>
#include <cmath>
#include <iomanip>
#include <fstream>

using namespace std;

// define our f(x) which we think is similar to the integral to be estimated.
double f_x(double x){
	return exp(-x*x);
}

// Define the integral f(x) from 0-1 to create the normalizing weighting func.
double weight = (exp(1) - 1) / exp(1);

// The weighting func itself
double g_x(double x){
	return exp(-x) / weight;
}

// The integral of the weigh func, inverse, to use a normal dist [0,1] on.
double G_inv(double u){
	return -log(1 - (u * (weight)));
}

int main(){

	double int1known = 0.6321205588285576784044762298385391325541888689682321654;
	double int2known = 0.9999546000702375151484644;
	weight = int2known;
	random_device rand;
	seed_seq seed{rand(), rand(), rand(), rand(), rand(), rand(), rand()};
	mt19937 eng(seed);
	uniform_real_distribution<> dist(0,1);
	

	// Initialize number of steps, and the sum at 0.
	int n = 5;
	double summa = 0;
	cout << setprecision(10) << weight << endl;
	for (int i = 0; i < n; i++){

		double u = dist(eng);
		cout << u << endl;
		double xval = G_inv(u);
		cout << "xvals: " << xval << endl;
		summa += f_x(xval)/g_x(xval);
		cout << "summas: " << summa << endl << endl;

	}

	double integ = summa/n;
	cout << "Integral estiamte: " << integ << endl;
	//cout << "Integral precision :" << abs(integ - 0.886227) << endl;

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
