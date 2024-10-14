#include "MCvar.h"
#include "MCSweeps.h"

// A type for the int array so that I don't call that whole line every time. 
typedef vector<vector<int>> spinvec;  

void SQIsing(int L, double T, int NWarmup, int Nmeas, int NStep, ofstream& out2){


	// INITIALIZATION OF SpinConf...

	// Start by initializing the deltaE exponenet values for use in MCSweeps...
	unordered_map<int, double> deltaE_exp;
	deltaE_exp[-8] = exp(8/T);
	deltaE_exp[-4] = exp(4/T);
	deltaE_exp[0]  = 1.0;  
	deltaE_exp[4]  = exp(-4/T);
	deltaE_exp[8]  = exp(-8/T);

	cout << deltaE_exp[-8] << endl;
	cout << deltaE_exp[-4] << endl;
	cout << deltaE_exp[-0] << endl;
	cout << deltaE_exp[4] << endl;
	cout << deltaE_exp[8] << endl;


	// SpinConf Initialized array to L size. 
	spinvec SpinConf(L, vector<int>(L));

	// Randomize the spins either -1 or 1 first pass.
	uniform_int_distribution<int> lindist(0,1);	

	for (int n = 0; n < L; n++){
		for (int m = 0; m < L; m++){
			SpinConf[n][m] = lindist(engine) * 2 - 1;
		}
	}


	// WARM-UP SWEEPS...

	MCSweeps(NWarmup, SpinConf, L, deltaE_exp);	
	cout << "WARMUP DONE..." << endl;
	// SWEEP, THEN DO_MEASUREMENT...

	MCvar<double> E, E_2, E_4, M, M_2, M_4;

	ofstream output("dat.txt");

	output << "Magnetisation" << endl;

	for (int meas = 0; meas < 1000*Nmeas; meas++){

		MCSweeps(NStep, SpinConf, L, deltaE_exp);

		do_measurement(output, L, SpinConf, E, E_2, E_4, M, M_2, M_4);
		cout << "MEASURING AT: " << meas << endl;
	}

	output.close();

	// BEGIN THE PRINTOUT...

	cout << "BEGINNING OUTPUT............" << endl;
	cout << "Average Energy: " 	<< E.Avrg() << " +/- " << E.Err_Avrg() << endl;	
	cout << "Average Magnetisation	: " << M.Avrg() << " +/- " << M.Err_Avrg() << endl;

	// Printing to datafile dat2.txt for quesion 2.6
	out2 <<  T << "," << E.Avrg() << "," << M_2.Avrg() << "," << E_2.Avrg() << "," <<  E.Err_Avrg() << "," << M_2.Err_Avrg() << "," << E_2.Err_Avrg() << "," << M_4.Avrg() << endl;
	




	

}

int main(){

	
	

	// PARAMETERS for 2.4
	
	// ofstream out_del("null.txt");
    // out_del << "Temp, Energy, Mag_Sq, Energy_Sq, Err_E, Err_E_Sq, Err_Mag_Sq" << endl;


	// int L = 50;
	// double T = 2.26;
	// int NWarmup = 5000;
	// int Nmeas = 5;
	// int NStep = 5;

	// SQIsing(L,T,NWarmup,Nmeas,NStep, out_del);

	// out_del.close();
	

	// PARAMETERS for 2.5

	// ofstream out2("dat2.txt");
    // out2 << "Temp, Energy, Mag_Sq, Energy_Sq, Err_E, Err_E_Sq, Err_Mag_Sq, Mag_Frth" << endl;
	
    // for (float temp = 2; temp <= 2.5; temp+=0.05){

    // 	int L = 10;
	// 	double T = temp;
	// 	int NWarmup = 2000;
	// 	int Nmeas = 5;
	// 	int NStep = 10;

	// 	SQIsing(L,T,NWarmup,Nmeas,NStep, out2);
    // }

	// out2.close();

	// PARAMETERS FOR 2.6
	
    
    for (int L = 8; L<= 16; L+=4){

    	string name = "dat_" + to_string(L) + ".txt";
    	ofstream out3(name);
    	out3 << "Temp, Energy, Mag_Sq, Energy_Sq, Err_E, Err_E_Sq, Err_Mag_Sq, Mag_Frth" << endl;
    
	    for (float temp = 2; temp <= 2.5; temp+=0.05){

			double T = temp;
			int NWarmup = 2000;
			int Nmeas = 5;
			int NStep = 10;

			SQIsing(L,T,NWarmup,Nmeas,NStep, out3);
	    }

	    out3.close();

	}

}