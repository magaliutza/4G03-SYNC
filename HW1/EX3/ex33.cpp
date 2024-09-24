//Bessel question
//3.3 

#include <iostream>
#include <iomanip>
#include <vector>
using namespace std;

double bessel(double x, int o, int N){
	
	// Function to calculate bessel func of order 'o' and 'x', with steps 'N'.

	// 'N' is a variable so we can change how many order above what we want to calculate we start at without replacing stuff in the code manually.
	// 'N' should be 11 because in testing, this is the upper limit for maximum precision and minimum computation power. Any more repetitions are a waste.
	//int N = 11;
	//change
	int order = N + o;

	double js[order+1];
	js[order] = 0;
	js[order-1] = 1;
	// We want the iterator in the loop to conform to n in the equation, so we start at n = N - 2 since we want n to be the first descended order,
	// n + 1 is the highest order which is 0, and n - 1 is the next decsended order.
	for (int n = order - 1; n > 0; n--){

		// J_(n-1)(x) = (2n/x)J_n(x) - J_(n+1)(x) --> TRANSLATE TO CODE!
		js[n-1] = (2*(n)/x)*js[n] - js[n+1]; 
		//cout << js[n-1] << endl;
	}
	// Let's do the normalization.
	
	int j = 2;
	double sum;
	
	sum += js[0];
	while(j<=order){
		sum+= 2*(js[j]);
		//cout << "J: " << j << " Sum: " << sum << endl;
		j+=2;
	}
	double besselval = js[o]/sum;
	return besselval;
}

int main(){

	// Initialize two variables for x and order of Bessel, and have a dialogue with the user to input which they'd like to calculate.

	double x;
	int n; 
	
	cout << "Bessel function estimator, using a descending recursion relation." << endl;
	cout << "First specify the x value to determine the bessel funciton at: ";
	cin >> x; 
	cout << "Bessel defined as J_n(" << x << ")." << endl;
	cout << "Specify the order desired to calculate: ";
	cin >> n; 
	cout << "You are searching for the Bessel function of order " << n << " at x-value " << x << " or, " << endl << "J_" << n << "(" << x << ")" << endl; 
	cout << "Result: " << setprecision(16) << bessel(0.5,0,11) <<endl;

	int k = 1;
	while(k<20){

		cout << "k: " << k << " TO: " << setprecision(16) << bessel(0.5,0,k) << endl;
		k += 1;

	}
}


