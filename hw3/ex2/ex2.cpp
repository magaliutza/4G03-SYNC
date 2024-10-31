#include <random>
#include <chrono>
#include <iostream>
#include <cmath>
#include <vector>
#include "hv.cc"
#include <fstream>
using namespace std;

typedef vector<double> Row; // One row of the matrix
typedef vector<Row> Matrix; // Matrix: a vector of rows

void hv(vector<double>& y, const vector<double>& x,int L);
void jacdiag(Matrix& A, vector<double>& d);

inline double dotprod(const vector<double>& y, const vector<double>& x){ 
	double sum = 0.0;
	for (int i=0;i<y.size();i++){
		sum += x[i] * y[i];
	}
	return (sum);
}




int main(){

	// define the geneator without seeding it.
	mt19937 generator;
	uniform_real_distribution<double> distribution(0,1.0);



	//int L,m;

	// cout << "Size of system, L: ";
	// cin >> L;
	// cout << "Size of Lanczos matrix, m: ";
	// cin >> m;
	

	int m = 50;
	
	ofstream output("energies_2.csv");
	output << "L,E0,E1" << endl;

	for (int L = 2; L <=24; L+=2){
		cout << "L = " << L << endl;
		int N=pow(2,L);

		vector<double> v1(N),v2(N),f(N),omega(N);
		Matrix Lan(m,Row(m)),v(m,Row(m));
		
		// my best recreation of the algorithm in the assignment.

		// |v0⟩ ← random vector
		for (int i=0;i<N;i++){
			v1[i] = 1.0-2.0*distribution(generator);
			v2[i]=0.0;
		}

		// normalizing v1 , |v0⟩ ← |v0⟩/||v0||
		double normalize = sqrt(dotprod(v1,v1));
		for (int j =0; j < N; j++){
			v1[j] = v1[j]/normalize;
		}

		// courteously premade lanczos generator, for |w⟩ ← A|v0⟩
		hv(omega, v1, L);

		// α0 ← ⟨v0|w⟩
		double alph = dotprod(v1,omega);

		// we finally calculated alpha, so lets set it as such!
		Lan[0][0] = alph;


		// |f0⟩ ← |w⟩ − α0|v0⟩
		for (int k = 0; k < N; k++){
			f[k] = omega[k] - alph * v1[k];
		}



		// for j = 0, . . . m − 2
		for (int n = 0; n <= m-2; n++){
			
			// βj ← ||fj||
			double beta = sqrt(dotprod(f,f));

			// |vj+1⟩ ← |fj ⟩/βj
			//v2[n] = if (beta = 0 ? )f[n]/beta;
			for (int md = 0; md <= N; md++){
				v2[md] = f[md] / beta;
			}

			// |w⟩ ← A|vj+1⟩ − βj|vj ⟩ (new omega time), create using hv from the next vector
			hv(omega,v2,L);
			// now actually set omega, we already have omega filled with v2 values from hv() call
			for (int dm = 0; dm <= N; dm++){
				omega[dm] -= beta * v1[dm];
			}

			// αj+1 ← ⟨vj+1|w⟩
			alph = dotprod(v2,omega);

			// |fj+1⟩ ← |w⟩ − αj+1|vj+1⟩
			for (int k = 0; k < N; k++){
				f[k] = omega[k] - alph * v2[k];
			}
			// Fill in the LAN tridiagonals with the calculated and normalized and orthogonalized values.
			Lan[n+1][n+1] = alph; 
			if (n < m - 1) {
	    		Lan[n][n + 1] = beta; 
	    		Lan[n + 1][n] = beta; 
			}	

			// update for the next cycle, when we move down the diagonal.
			v1 = v2;
		}

		Row d(m);

		jacdiag(Lan,d);
		sort(d.begin(), d.end());

		cout << setprecision(18) << "0th Energy Level: " << d[0] << endl;
		cout << setprecision(18) << "1st Energy Level: " << d[1] << endl << endl;
		output << L << "," << d[0] << "," << d[1] << endl;
	}
	output.close();
	return 0;
}