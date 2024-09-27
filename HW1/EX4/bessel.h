double bessel(double x, int o, int N){

	int order = N + o;

	double js[order+1];
	js[order] = 0;
	js[order-1] = 1;
	// We want the iterator in the loop to conform to n in the equation, so we start at n = N - 2 since we want n to be the first descended order,
	// n + 1 is the highest order which is 0, and n - 1 is the next decsended order.
	for (int n = order - 1; n > 0; n--){

		// J_(n-1)(x) = (2n/x)J_n(x) - J_(n+1)(x) --> TRANSLATE TO CODE!
		js[n-1] = (2*(n)/x)*js[n] - js[n+1]; 
	}
	// Let's do the normalization.
	
	int j = 2;
	double sum;
	
	sum += js[0];
	while(j<=order){
		sum+= 2*(js[j]);
		j+=2;
	}
	double besselval = js[o]/sum;
	return besselval;
}
