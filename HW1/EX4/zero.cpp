#include <iostream>
#include <math.h>
#include <iomanip>
#include <cmath>
#include "bessel.h"
using namespace std;


int main(){

	// Reasonable starting values, by inspection, are (2,3), (5,6), and (8,9).
	double starting[6] = {2,3,5,6,8,9};

	int j = 1;
	// Looping through all 3 starting values.
	for (int i=0;i<6;i+=2){
		// Initialize some variables, the starting a and b, the tolerance, and the centerpoint.
		double a = starting[i];
		double b = starting[i+1];

		double tolerance = 0.00000000001;
		double centerpt = 0;

		// Make sure to avoid runaways, we keep it under 15 decimals (the limit of double).
		while (b-a > tolerance){

			// Calculate the centerpoint.
			centerpt = (b+a)*0.5;
			
			// Check if the centerpoint is magically 0!
			if (bessel(centerpt,0,11)==0){
				break;
			}
			// Check if the sign change is on between the center and a (the only way this result would be negative is a + * - , thus indicating a crossing of the 0 line).
			else if (bessel(a,0,11)*bessel(centerpt,0,11) < 0){
				// Shrink the bisection towards a if true.
				b = centerpt;
			}
			// If its not anything before, it must be somewhere between center and b!
			else{
				a = centerpt;
			}
		}

		cout << setprecision(15) << "Zero of 0th order Bessel function number " << j << " is calculated to be: " << centerpt << endl;
		j++;
	}
}