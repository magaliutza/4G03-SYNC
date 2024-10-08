#ifndef MCSWEEPS_H
#define MCSWEEPS_H

#include <vector>
#include <unordered_map>
#include <random>

using namespace std;

// Type definition for the spin configuration
typedef vector<vector<int>> spinvec;

extern random_device ran; 
extern seed_seq seed;
extern mt19937 engine;

// Function prototypes
void MCSweeps(int sweeps, spinvec& SpinConf, int L, const unordered_map<int, double>& deltaE_exp);
double local_Energy(int x, int y, const spinvec& SpinConf, int L);
double Energy_lattice(const spinvec& SpinConf, int L);
int Magnet_lattice(const spinvec& SpinConf, int L);
void do_measurement(int L, const spinvec& SpinConf, MCvar<double>& E, MCvar<double>& E_2, MCvar<double>& E_4, MCvar<double>& M, MCvar<double>& M_2, MCvar<double>& M_4);
int per(int coord, int L);

#endif // MCSWEEPS_H
