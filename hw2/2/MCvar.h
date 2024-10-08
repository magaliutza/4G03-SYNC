#include <iostream>
#include <algorithm>
#include <math.h>
#include <cmath>
#include <random>
#include <iomanip>
#include <fstream>
#include <array>
#include <vector>
#include <unordered_map>
using namespace std;

template <class Type> class MCvar {

// Data
public:

private:
	int bin_size;
    int count;
    int no_of_bins;
	Type av;
	vector<double> Bins;  

//Methods
public:
        MCvar() : bin_size(1000),count(0),no_of_bins(0) {}
        MCvar(int bsz) : bin_size(bsz),count(0),no_of_bins(0) {}
        void push(const Type v){

                if (count == 0) {
                    av = v;
                    count += 1;
                }

                else if (count < bin_size-1) {
                    av += v;
                    count += 1;
                }

                else {
                    av += v;
                    Bins.push_back(((double) av)/bin_size);
                    count = 0;
                    no_of_bins += 1;
                }
        }
        double Avrg(){

                double res = 0.0;
                for (auto r : Bins) res += r;
                return (res/no_of_bins);
        }

        double SqrAvrg(){
            double res2 = 0.0;  
            for (auto r : Bins) res2 += r*r;
            return (res2/no_of_bins);

        }

        double Err_Avrg(){

                double res = 0.0, res2=0.0;
                
                // use our funcs with avg and :

                double res_avg_sqr = Avrg()*Avrg();
                double res2_avg = SqrAvrg();

                double SDOM = sqrt( (res2_avg - res_avg_sqr) / (no_of_bins - 1) );

                return SDOM;
        }
};
