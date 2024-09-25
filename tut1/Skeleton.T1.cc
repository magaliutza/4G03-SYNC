#include <iostream>
#include <math.h>
#include <random>
#include <cstddef>
#include <cstdlib>
#include <iostream>

using namespace std;

int main() {

//////1.1: test random number generator
  // Seed with a real random value, if available
    #include <random>
    random_device r;
    std::uniform_int_distribution<int> uniform_dist(1, 2);

    std::mt19937 e2;
    std::uniform_real_distribution<> my_r(0,1);
    cout << "HI!" << uniform_dist(r);
    cout << my_r(e2);

//////1.2: define lambda function with mu=5.75, sigsq=0.01
  // Explicitly create a lambda function for the normal distribution
  const double mu = 5.75;
  const double sigsq = 0.01;
  auto normdist = [mu,sigsq] (double x) {
    return /*Define lambda function, Fillng your code here*/;
  };

//////1.3: Lets do a naive Monte Carlo
  //generate random number between interval [-30,30]
  double const xlim_low = -30.0;
  double const xlim_high = 30.0;
  std::uniform_real_distribution<> gen2(xlim_low,xlim_high);

  double x;
  const int N_Meas1 = 100000;

  double x_avrg=0.0;
  double x2_avrg=0.0;

  /*
    Evaluate moments using Simple Monte Carlo
    Filling youre code here
    ...
    ...
    ...
  */

  x_avrg = x_avrg*(xlim_high-xlim_low)/N_Meas1;
  x2_avrg = x2_avrg*(xlim_high-xlim_low)/N_Meas1;
  
  cout << "Naive MC Average: " << x_avrg << endl;
  cout << "Naive MC <x*x>: " << x2_avrg << endl;
  cout << "Naive Variance <(x-mu)*(x-mu)>: " << x2_avrg-x_avrg*x_avrg << endl;
  cout << "mu:" << mu << endl;
  cout << "sigsq:" << sigsq << endl;
  cout << endl;

//////1.4: Random walk

  double xnew;
  const double xstep =0.001;
  const int N_Warmup = 10000000;
  const int N_Step = 1000;
  const int N_Meas2 = 100000;

  //warm up
  x = 5.0;
  /*
    Warm up to reach equilibrium
    Filling your code here
    ...
    ...
    ...
  */
  cout << "Starting measurements from position: " << x << endl;

  x_avrg=0.0;
  x2_avrg=0.0;
  /*
    Perform the random walk
    Filling your code here
    ...
    ...
    ...
  */

  x_avrg /= N_Meas2;
  x2_avrg /= N_Meas2;
  cout << "True MC Average: " << x_avrg << endl;
  cout << "True MC <x*x>: " << x2_avrg << endl;
  cout << "True Variance <(x-mu)*(x-mu)>: " << x2_avrg-x_avrg*x_avrg << endl;
  cout << endl;

  return 0;
}
