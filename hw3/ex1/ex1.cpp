#include <cmath>
#include <vector>
#include <iostream>
#include <iomanip>
#include <algorithm>
using namespace std;
// MOO DENG MOMENT (MDM)
// setting the types for a matrix
typedef vector<double> Row;
typedef vector<Row> Matrix;



// basic definitions because im absurdly lazy (And traumatized from manually setting equations.)
double alphcalc(double aqq, double app, double apq){
		return (aqq-app)/(2*apq);
	}

double tcalc(double alph){
	return (alph <0 ? -1 : 1)/(abs(alph)+sqrt(alph*alph+1));
}

double ccalc(double t){
	return 1/(sqrt(t*t+1));
}

double scalc(double t, double c){
	return t*c;
}

double taucalc(double s, double c){
	return s / (c+1);
}

void jacdiag(Matrix& A, vector<double>& d){

	double thresh = 2.2e-16;
	int len = A.size();

	for (int sweeps = 0; sweeps <=100; sweeps++){

		// sum checker so we don't sweep unnecessarily (if we've converged to 0 already...)
		double diagsum = 0;
		for (int p = 0; p < len; p++){
			for (int q = p + 1; q < len; q++){
				diagsum += A[p][q];
			}
		}	
		cout << "Diagonal element sum at sweep " << sweeps << " = " << diagsum << endl;		
		if (diagsum ==0){
			break;
		}


		// set the threshold values depending on what sweep number we're on , based on handout
		if (sweeps>4){
			thresh = 2e-25;
		}
		// If an entire sweep occurs with no rotations, we crunch the loop.
		bool significant = false;

		for (int p = 0; p < len; p++){
			for (int q = p + 1; q < len; q++){
				
				// break check for ignoring low value rotations.
				// if (abs(A[p][q])<thresh){
				// 	//cout << "BROKEN!" << endl;
				// 	continue;
				// }

				significant = true;

				// pull the values of App, Aqq, and Apq for calcs
				double app = A[p][p];
				double aqq = A[q][q];
				double apq = A[p][q];

			
				// calculate the necessary variables. some might just be to calculate others 
				// and im just lazy
				double alph = alphcalc(aqq, app, apq);
				double t = tcalc(alph);
				double c = ccalc(t);
				double s = scalc(t,c);
				double tau = taucalc(s,c);
				//cout << "tau at " << p << "," << q << " is " << tau << endl;

				// the golden moment; adjust the diagonals and set the chosen point to 0.
				A[p][p] = app - t*apq;
				A[q][q] = aqq + t*apq;
				A[p][q] = 0;
				A[q][p] = 0;
				// deal with the rotation of the row and column.
				for (int i = 0; i < len; ++i) {
                    if (i != p && i != q) {
                        double aip = A[i][p];
                        double aiq = A[i][q];
                        A[i][p] = aip - s * (aiq + tau * aip);
                        A[i][q] = aiq + s * (aip - tau * aiq);
                        A[p][i] = A[i][p];
                        A[q][i] = A[i][q];
                    }
                }
			}
		}
		
		// if (!significant){
		// 	break;
		// }

		// square sum printout as demanded by 1.2 
		double squaresum = 0;
		for (int p = 0; p < len; p++){
			for (int q = p + 1; q < len; q++){
				squaresum += A[p][q]*A[p][q];
			}
		}
		//cout << "square sum at sweep " << sweeps << " = " << squaresum << endl;	
	}

	for (int k = 0; k < len; k++) {
        	d[k] = A[k][k];
	}
}

void Hnm(Matrix& H){
	int len = H.size();
	for (int n=1; n <= len; n++){
		for (int m=1; m <= len; m++){
			double val = (n == m ? 1 : 0)*(n*n) + ( 0.05*min(n,m) + 5*pow(-1,abs(m-n))*min(n,m));
			H[n-1][m-1] = val;
		}
	}
}

int main(){

	// Matrix A = {{1.0, -2.0, 3.0},{-2.0,4.0,6.0},{3.0,6.0,9.0}};
	// Row d(3);

	// jacdiag(A,d);

	// cout << "Eigenvalues: " << endl;
	// for (auto val : d){
	// 	cout << setprecision(18) << val << endl;
	// }
	
	// return 1;

	for (int size = 10; size <=40; size+=10){
		cout << "M = " << size << endl;
		cout << "Jacobi for M = " << size << endl;
		Matrix H(size, Row(size,0.0));
		Hnm(H);

		Row d(size);
		jacdiag(H,d);

		sort(d.begin(), d.end());
		for (int i=0; i<3; i++){
			cout << "Eigenvalue " << i+1 << ": " << setprecision(15) << d[i] << endl;
		}
		

		cout << endl << endl << endl << endl;
	}

	
	
	

}
