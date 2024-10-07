#include "MCvar.h"
using namespace std;

// USEFUL FUNCTION DEFINITIONS AND TYPE DEFINITIONS

// Periodization correction function to deal with boundary conditions.
int per(int coord, int L){
	if (coord < 0)
		return L-1;
	else if (coord >= L)
		return 0; 

	return coord;
}


// A type for the int array so that I don't call that whole line every time. 
typedef vector<vector<int>> spinvec;  

// Function to calculate the local energy from the neighbours, to be called upon later. 
double local_Energy(int x, int y, spinvec& SpinConf, int L){

	double neighbour_energy = 0;

	// Spin at the site investigated.
	int spin = SpinConf[x][y];

	neighbour_energy += -1 * spin * ( 
		SpinConf[x][per(y+1,L)] + 
		SpinConf[x][per(y-1,L)] + 
		SpinConf[per(x+1,L)][y] + 
		SpinConf[per(x-1,L)][y]
		); 

	return neighbour_energy;
}

// PROBABILITY DISTRIBUTIONS AND RNG ENGINES...

// Making a random engine.
random_device rand; 
seed_seq seed{rand(), rand(), rand(), rand()};
mt19937 engine(seed);
 

void MCSweeps(int sweeps, double J, spinvec& SpinConf, int L, const unordered_map<int, double>& deltaE_exp){


// Random coordinate generator
uniform_int_distribution<int> coordrand(0,L-1);
// Linear dist. random from 0 to 1
uniform_real_distribution<> linrand(0,1);

// Begin the sweeps defined by the number of sweeps demanded.
for (int n = 0; n < sweeps; n++){

	// Each sweep is LxL flip tries.
	int flips = L*L;

	for (int i = 0; i < flips; i++){
		
		// Get random coordinates 
		int x = coordrand(engine);
		int y = coordrand(engine);

		// Find the local energy before a flip.
		double energy_before = local_Energy(x,y,SpinConf,L);
		// Flip the bit. Remember we flipped it. 
		SpinConf[x][y] *= -1;
		// Find local energy after a flip.
		double energy_after = local_Energy(x,y,SpinConf,L);
		// Delta E calc.
		double deltaE = energy_after - energy_before;
		// Lookup the exp with the energy value found.
		double exp = deltaE_exp[deltaE];

		double Prob_of_flip = min(1.0,exp);

		double r = linrand(engine);

		// Remember the spin is already flipped. So if the flip check FAILS, we have to flip it back. Otherwise, move on.
		if  (Prob_of_flip <= r)
			SpinConf[x][y] *= -1;

	}

}

// To calculate the energy of the whole lattice, we find the local energy for each position and divide by 2 to take care of double counting.
double Energy_lattice(const spinvec& SpinConf, int L){
	
	double energy = 0;

	for (int i=0; i<L; i++){
		for (int k=0; k<L; k++){
			energy += local_Energy(i,k,SpinConf, L);
		}
	}

	energy = energy/2;
	return energy;
}

// Same idea for magnetization.
int Magnet_lattice(const spinvec& SpinConf, int L){
	int magnet = 0;

	for (int i=0; i<L; i++){
		for (int k=0; k<L; k++){
			magnet += SpinConf[i][k];
		}
	}
	return magnet;
}

