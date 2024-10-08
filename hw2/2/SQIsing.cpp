#include "MCvar.h"
#include "MCSweeps.h"

// A type for the int array so that I don't call that whole line every time. 
typedef vector<vector<int>> spinvec;  

// Making a random engine from rand seed.
// random_device ran; 
// seed_seq seed{ran(), ran(), ran(), ran()};
// mt19937 engine(seed);

void printSpinConfiguration(const spinvec& SpinConf) {
    cout << "Spin Configuration:" << endl;
    
    // Print the dimensions
    cout << "Array Dimensions: " << SpinConf.size() << "x" << SpinConf[0].size() << endl;

    // Print the spin configuration in a grid format
    for (const auto& row : SpinConf) {
        cout << "| "; // Add a border for better visibility
        for (int spin : row) {
            cout << spin << " | "; // Print each spin with a border
        }
        cout << endl; // Move to the next line after each row
        cout << "---------------------" << endl; // Separator line
    }
}

   

void SQIsing(int L, double T, int NWarmup, int Nmeas, int NStep){



	double k = 1.380649e-23;
	//T = 1;

	// INITIALIZATION OF SpinConf...

	// Start by initializing the deltaE exponenet values for use in MCSweeps...
	unordered_map<int, double> deltaE_exp;
	deltaE_exp[-8] = exp(-8/T);
	deltaE_exp[-4] = exp(-4/T);
	deltaE_exp[0]  = 1.0;  
	deltaE_exp[4]  = exp(4/T);
	deltaE_exp[8]  = exp(8/T);

	cout << deltaE_exp[-8] << endl;
	cout << deltaE_exp[-4] << endl;
	cout << deltaE_exp[-0] << endl;
	cout << deltaE_exp[4] << endl;
	cout << deltaE_exp[8] << endl;


	// SpinConf Initialized array to L size. 
	spinvec SpinConf(L, vector<int>(L));

	// Randomize the spins either -1 or 1 first pass.
	uniform_int_distribution<> lindist(0,1);	

	for (int n = 0; n < L; n++){
		for (int m = 0; m < L; m++){
			SpinConf[n][m] = lindist(engine) * 2 - 1;
			//cout << lindist(engine);
			//cout << SpinConf[n][m] << endl;
		}
	}
	printSpinConfiguration(SpinConf);
	//cout << "INITIALIZATION SpinConf" << endl;


	// WARM-UP SWEEPS...

	MCSweeps(NWarmup, SpinConf, L, deltaE_exp);

	printSpinConfiguration(SpinConf);
	//cout << "WARMUP SpinConf" << endl;

	// SWEEP, THEN DO_MEASUREMENT...

	MCvar<double> E, E_2, E_4, M, M_2, M_4;

	for (int meas = 0; meas < Nmeas; meas++){

		MCSweeps(NStep, SpinConf, L, deltaE_exp);

		do_measurement(L, SpinConf, E, E_2, E_4, M, M_2, M_4);
	}

	// BEGIN THE PRINTOUT...

	cout << "BEGINNING OUTPUT............" << endl;
	cout << "Average Energy: " << E.Avrg() << " +/- " << E.Err_Avrg() << endl;
	cout << "Average Magnetisation: " << M.Avrg() << " +/- " << M.Err_Avrg() << endl;


	




	

}

int main(){

	int L = 10;
	double T = 2.5;
	int NWarmup = 1000;
	int Nmeas = 100;
	int NStep = 100;

	SQIsing(L,T,NWarmup,Nmeas,NStep);

}