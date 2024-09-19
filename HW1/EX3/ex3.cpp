// Exercise 3 Bessel Funcs
// 3.1 



// EQ 5: J_(n+1)(x) = (2n/x)*J_n(x)-J_(n-1)(x)

#include <iostream>
#include <iomanip>
#include <vector>

int main(){

	//initiate the first 2 bessels at x = 0.5, given. 
	//double j0 = 0.938469807240813; 
	//double j1 = 0.2422684577;

	std::vector<double> js;
    js.push_back(0.938469807240813);
    js.push_back(0.2422684577);

	std::cout << std::setprecision(15) << js[0] << std::endl;
    double jnew;
	for (int i=0; i<20;i++){
		jnew = (2*i/0.5)*js[i+1]-js[i];
		js.push_back(jnew);

        }
	std::cout << "20th Bessel using a ascending manner: " <<jnew << std::endl;
//Not a great approximation compared to online sources. 

std::cout << "duckmount worked as intended" << std::endl;
//commenting for git recognition

}