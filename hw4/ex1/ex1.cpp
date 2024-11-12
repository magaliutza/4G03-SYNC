#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;

// MDM3 !!!! Welcome to my program I am NOT going insane. 


// define our potential so we can summon it for ksqr (more readable this way, also can use it externally if needed)
// of course we preserve units of 2m/hbar^2, because why would we work with silly units other than the natural ones!
// double V(double x){
// 	if (x >= -5 and x <= 5)
// 		return (x * x) / 5.0  - 14.0;
// 	else
// 		return 0.0;
// }

// test potential which worked very well thank you.
// double V(double x){
// 	int lambda = 4;
// 	double pot = lambda * (lambda - 1) * (1.0/2.0 - 1/(cosh(x)*cosh(x)));
// 	return pot;
// }

// second potential for ex 1.6

double V(double x){
	return (0.05/(sin(x)*sin(x))) + (5.0/(cos(x)*cos(x)));
}

// lets define k^2 for numerov's method!
// again preserving units and ignoring the hbar and 2m.
// DOUBLE CHECK THIS BEFORE DEBUGGING
double ksqr(double x, double E){
	return E - V(x);
}


// possibly going to be used to automate the guesses for E based on the infinite square well of an identical size. DEFUNCT
// use E_n = n^2hbar^2/8mL^2, but pull out a factor of hbar^2/2m to become n^2/4L^2
// double infE(double L){
// 	return 1.0/(4*L*L);
// }

// not used by the program, but leave in for updating. 
double turning(double E = -9){
	return sqrt(5*(E+14));
}

double numerov(vector<double>& phi_left, vector<double>& phi_right, vector<double>& x, double E, int n, double h){

	// define the turning point. Might nomrally use turning(), but it sucks a bit. 5 works best.
	//double xm = turning();
	double xm = 1.5;
	

	//cout << "BEGINNING NUMEROV AT E = " << E << endl;

	
	// start with left to right integration. Yes I baked it into once function.

	// set up ksqrs ahead of time and calculate only once...
	vector<double> ksqrlist(n);
	for (int i = 0; i < n; i++){
		ksqrlist[i] = ksqr(x[i], E);
	}

	// set the primary two values to zero and 'some small number' which we know doesnt REALLY matter
	phi_left[0] = 0;
	phi_left[1] = 1e-5;

	// start the numerov integration, for left to right integration
	for (int i = 1; i < n - 1; i++){
		
		// stop when we get to the turning point.
		if (x[i] >= xm){
			break;
		}
		
		// find the k^2's we need for i-1, i, i+1...
		double ksqrmin1 = ksqrlist[i-1];
		double ksqri = ksqrlist[i];
		double ksqrplus1 = ksqrlist[i+1];

		// use numerov expression to directly calculate the next value.
		phi_left[i+1] = (2.0*(1.0 - (5.0/12.0)*h*h*ksqri)*phi_left[i] - (1.0 + h*h/12.0*ksqrmin1)*phi_left[i-1])/(1.0 + h*h/12.0*ksqrplus1);
	}

	// lets do the same for right to left!
	phi_right[n-1] = 0;
	phi_right[n-2] = 1e-5;

	for (int i = n - 2; i > 0; i--){
		
		// stopping condition.
		if (x[i] <= xm)
			break;	
		
		// already calculated, just pull from the vector. 
		double ksqrmin1 = ksqrlist[i-1];
		double ksqri = ksqrlist[i];
		double ksqrplus1 = ksqrlist[i+1];

		// numerov it, but in reverse!
		phi_right[i-1] = (2.0*(1.0-5.0/12.0*h*h*ksqri)*phi_right[i] - (1.0+h*h/12.0*ksqrplus1)*phi_right[i+1])/(1.0+h*h/12.0*ksqrmin1);

	}

	// lets rescale the left to the right based on the matching point

	// xm is cast as a double and wont refer to the index where the turning point is.
	// xmi finds the index where xm is using the evenly spaced array.
	int xmi = static_cast<int>((xm-x[0])/h+0.5);

	// find how much we need to scale by!
	double scaler = phi_right[xmi] / phi_left[xmi];
	for (int j = 0; j <= xmi; j++){
		phi_left[j] *= scaler;
	}

	// caculate G(E) and RETURN it out of numerov, this is what we show for all our work!
	double G = (phi_left[xmi + 1] - phi_left[xmi - 1]) - (phi_right[xmi + 1] - phi_right[xmi - 1]);
	return G;

}

double bisection(double min, double max, double tol, vector<double>* G_E, vector<double>* x, int n, double E_step, double E_min, double E_max, double xm = 5){
	// Initialize some variables, the starting a and b, the tolerance, and the centerpoint.
		// since I copied this RIGHT out of the previous assignment, I just recast the variable to what they are in this assignment.
		double a = min;
		double b = max;

		double tolerance = tol;
		double centerpt = 0;
		int index = 0;
		cout << "prior to trigger" << endl;
		// stolen STRAIGHT from assignment 1 ex 4, modified to this program.
		while (abs(b-a) > tolerance){
			cout << index << endl;
			index +=1;
			// Calculate the centerpoint.
			centerpt = (b+a)*0.5;
			
			// lots of work to convert numbers into indicies, but works out well since we know the spacing!
			// int index_a = static_cast<int>((a - (*x)[0]) / h + 0.5);
			// int index_b = static_cast<int>((b - (*x)[0]) / h + 0.5);
			// int index_center = static_cast<int>((centerpt - (*x)[0]) / h + 0.5);

			int index_a = static_cast<int>((a - E_min) / E_step + 0.5);
			int index_b = static_cast<int>((b - E_min) / E_step + 0.5);
			int index_center = static_cast<int>((centerpt - E_min) / E_step + 0.5);

			// some brutal looking code to make sure we remain INSIDE the bounds of things, otherwise we throw some segmentation errors.
			index_a = std::max(0, std::min(index_a, static_cast<int>(G_E->size() - 1)));
			index_b = std::max(0, std::min(index_b, static_cast<int>(G_E->size() - 1)));
			index_center = std::max(0, std::min(index_center, static_cast<int>(G_E->size() - 1)));


			// Check if the centerpoint is near enough...
			if (abs((*G_E)[index_center]) < tol)
				return (*x)[index_center];

			// Check if the sign change is on between the center and a (the only way this result would be negative is a + * - , thus indicating a crossing of the 0 line).
			else if ((*G_E)[index_a] * (*G_E)[index_center] < 0)
				// Shrink the bisection towards a if true.
				b = centerpt;
			// If its not anything before, it must be somewhere between center and b!
			else
				a = centerpt;
			
		}
	// make sure something is returned, worst case scenario.
	cout << "worst case" << endl;
	return (b + a) * 0.5;
}


int main(){

	// initialize certain vars. Tolerance, steps n, size L, stepsize, and initial energy guess.
	double tol = 0.001;
	int n = 1000;
	double L = 3.14159265/2.0;
	double step = (L - 2 * 1e-5) / n; 
	//double E = -9.0;

	// initialize the x array since it is unchanging between energy values. 
	vector<double> x(n);
	for (int i = 0; i < n; i++){
		x[i] = 1e-5 + i*step;
	}

	double E_min = 0;
	double E_max = 50;
	double E_step = 0.0005;
	int Gsize = (E_max - E_min)/E_step;
	
	vector<double> G_E(Gsize);
	vector<double> E_arr(Gsize);
	

	// open a file for G(E) and E, write a header. 
	ofstream file("ge.txt");
	file << "E, G" << "\n";

	// index so i can place G(E)'s in an array while looping over E's instead of i's
	int index = 0;
	
	// start a loop to numerov different energy values. These return their G(E) value at the end, so well save and plot that.
	for (double E = E_min; E <=E_max; E+=E_step){
		vector<double> phi_left(n);
		vector<double> phi_right(n);
		E_arr[index] = E;

		
		//cout << "E: " << E << endl;


		// numerov at energy E
		double G = numerov(phi_left, phi_right, x, E, n, step);
		G_E[index] = G;

		//cout << "G: " << G << endl;
		file << E << "," << G << "\n";

		
		index+=1;
		


	}

	file.close();

	double eig_min, eig_max;
	// ask user for bracketing to find an eigenvalue more precisely with BISECTION.
	cout << "bracket an eigenvalue to find it more precisely..." << endl;
	cout << "Min val?" << endl;
	cin >> eig_min;
	cout << "Max val?" << endl;
	cin >> eig_max;

	double eigenval;
	eigenval = bisection(eig_min, eig_max, tol, &G_E, &x, n, E_step, E_min, E_max);

	cout << "Eigenvalue found at: " << eigenval << endl;
	// second file for the eigenfuncs.
	ofstream eig_file("eigenfunc.txt");
	eig_file << "x, phil, phir" << "\n";

	//eigenval = 47.8691;

	vector<double> phi_left(n);
	vector<double> phi_right(n);


	// numerov at energy E for the eigenavls.
	numerov(phi_left, phi_right, x, eigenval, n, step);

	for (int i = 0; i < n; i++){
		eig_file << x[i] << "," << phi_left[i] << "," << phi_right[i] << "\n";
	}
	eig_file.close();




	cout << "I DID THE THING." << endl;
	// initiate arrays 

	


	//cout << infE(L) << endl;
	//cout << turning() << endl;


	//cout << numerov(phi_left, phi_right, x, E, n, step) << endl;


	//cout << V(10) << ", " << V(1) << endl;
	return 0;
}