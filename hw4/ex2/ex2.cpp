#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
using namespace std;

typedef vector<double> Row;
typedef vector<Row> Matrix;

double pi = 3.14159265358979323846264338327950288419716939937510;

// calcualte the analytical solution with that very odd presnetaion of how to get it...
double theoretical(double x, double y, double L){
	double sum = 0;
	// no real reason to loop 100 times, just waht i thought was close enough to infinity to pass.
	for (int n = 1; n <= L; n+=2){
		double add = (400.0/n/pi)*sin(n*pi*x)*sinh(n*pi*y)/sinh(n*pi);
		sum += add;
	}

	return sum;
}

double jacobi(Matrix& U_old, Matrix& U_new, int n){

	double topdiff = 0;
	// iterate overall x and y in the interior grid (the outside should be set up to the figure) also I know x and y are TECHNICALLY swapped but it REALLY doesn't matter.
	for (int x = 1; x <= n - 2; x++){
		for (int y = 1; y <= n -2; y++){
			U_new[x][y] = 1.0/4.0 * (U_old[x+1][y] + U_old[x-1][y] + U_old[x][y+1] + U_old[x][y-1] );

			// find the difference on this point, then update the maximum difference as needed.
			double herediff = fabs(U_new[x][y] - U_old[x][y]);
			if ( herediff > topdiff){
				topdiff = herediff;
			}
		}
	}

	return topdiff;
}


int main(){

	// variables.
	double tol = 1e-6;
	int L = 100;
	double dd = 1.0 / (L - 1);

	// initialize two matricies, old and new with 0's everywhere.
	Matrix U_old(L, Row(L, 0.0)), U_new(L, Row(L, 0.0));

	// condition to set up the upper potential of 100, as in the diagram.
	for (int x = 0; x < L; x++){
		U_old[L-1][x] = 100;
		U_new[L-1][x] = 100;
	}
	// if i dont initialize this as greater than tolerance, the program just wont work.
	double topdiff = 1;
	// prevent infinite loops.
	int iter = 0;

	while (topdiff > tol){

		// perform the jacobi sweep.
		double max_diff = jacobi(U_old, U_new, L);
		// swap the old for the new
		swap(U_old, U_new);

		topdiff = max_diff;
		//cout << "topdiff = " << topdiff << endl;
		
		// lets make sure we dont go on forever, worst case scenario.
		iter += 1;
		//cout << "ITER: " << iter << endl;
		if (iter >= 100000){
			cout << "broken at max iterations!" << endl;
			break;
		}

	}
	// put this into a file for plotting.
	ofstream output("Umid.txt");

	for (int i = 0; i < L; i++){
		for (int j = 0 ; j < L; j++){
			output << U_new[i][j];
			if (j < L - 1){
				output << ","; 
			}
		}
		output << endl;
	}


	output.close();

	cout << "DID THE THING." << endl;


	// theoretically plotting. make an empty matrix, populate it, shovel it out to another file. BEGONE.
	Matrix U_theo(L, Row(L,0.0));

	for (int i = 0; i < L; i++){
		for (int j = 0; j < L; j++){
			// gotta convert indicies to values for the function, so we use the spacing
			double x = i * dd;
			double y = j * dd;
			U_theo[i][j] = theoretical(y,x,L);
		}
	}
	
	// identical to the above output, but to another file. 
	ofstream output2("Utheo.txt");

	for (int i = 0; i < L; i++){
		for (int j = 0 ; j < L; j++){
			output2 << U_theo[i][j];
			if (j < L - 1){
				output2 << ","; 
			}
		}
		output2 << endl;
	}

	output2.close();

	cout << "I DID THE OTHER THING." << endl;
}