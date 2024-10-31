#define BIT_SET(a,b) ((a) |= (1U<<(b)))
#define BIT_CLEAR(a,b) ((a) &= ~(1U<<(b)))
#define BIT_FLIP(a,b) ((a) ^= (1U<<(b)))
#define BIT_CHECK(a,b) ((bool)((a) & (1U<<(b))))
// Set on the condition f else clear
//bool f;         // conditional flag
//unsigned int m; // the bit mask
//unsigned int w; // the word to modify:  if (f) w |= m; else w &= ~m; 
#define COND_BIT_SET(a,b,f) ((a) = ((a) & ~(1U<<(b))) | ((-(unsigned int)f) & (1U<<(b))))
#include <vector>
#include <cmath>
#include <iostream>
using namespace std;


void hv(vector<double>& y, const vector<double>& x,int L)
{ 
  for (vector<double>::iterator it = y.begin() ; it != y.end(); ++it)
                    *it=0.0;
  bool b ;
  unsigned int k;
  for (unsigned int i=0;i<x.size();i++){
    if (abs(x[i])>2.2e-16) {
      int jm = L-1;
      double xov2=x[i]/2.0;
      for (int j=0;j<L;j++){
        k=i;
        COND_BIT_SET(k,jm,BIT_CHECK(i,j));
        COND_BIT_SET(k,j,BIT_CHECK(i,jm));
        y[k] += xov2;
        jm = j;
      }
    }
  }
  for (unsigned int i=0;i<x.size();i++)
    y[i]=y[i]-((double) L)/2.0*x[i]/2.0;
}


// copy pasting the jacobi method!!!

#include <iomanip>
#include <algorithm>
using namespace std;
// MOO DENG MOMENT 2 (MDM2)
// setting the types for a matrix
typedef vector<double> Row;
typedef vector<Row> Matrix;



// basic definitions because im absurdly lazy (And traumatized from manually setting equations.)
double alphcalc(double aqq, double app, double apq){
    return (aqq-app)/(2*apq);
  }

double tcalc(double alph){
  return (alph <0 ? -1 : 1)/(abs(alph)+sqrt(alph*alph+1));
}

double ccalc(double t){
  return 1/(sqrt(t*t+1));
}

double scalc(double t, double c){
  return t*c;
}

double taucalc(double s, double c){
  return s / (c+1);
}

void jacdiag(Matrix& A, vector<double>& d){

  double thresh = 2.2e-16;
  int len = A.size();

  for (int sweeps = 0; sweeps <=100; sweeps++){

    // sum checker so we don't sweep unnecessarily (if we've converged to 0 already...)
    double diagsum = 0;
    for (int p = 0; p < len; p++){
      for (int q = p + 1; q < len; q++){
        diagsum += A[p][q];
      }
    } 
    //cout << "Diagonal element sum at sweep " << sweeps << " = " << diagsum << endl;   
    if (diagsum ==0){
      break;
    }


    // set the threshold values depending on what sweep number we're on , based on handout
    if (sweeps>4){
      thresh = 2e-25;
    }
    // If an entire sweep occurs with no rotations, we crunch the loop.
    bool significant = false;

    for (int p = 0; p < len; p++){
      for (int q = p + 1; q < len; q++){
        
        // break check for ignoring low value rotations.
        // if (abs(A[p][q])<thresh){
        //  //cout << "BROKEN!" << endl;
        //  continue;
        // }

        significant = true;

        // pull the values of App, Aqq, and Apq for calcs
        double app = A[p][p];
        double aqq = A[q][q];
        double apq = A[p][q];

      
        // calculate the necessary variables. some might just be to calculate others 
        // and im just lazy
        double alph = alphcalc(aqq, app, apq);
        double t = tcalc(alph);
        double c = ccalc(t);
        double s = scalc(t,c);
        double tau = taucalc(s,c);
        //cout << "tau at " << p << "," << q << " is " << tau << endl;

        // the golden moment; adjust the diagonals and set the chosen point to 0.
        A[p][p] = app - t*apq;
        A[q][q] = aqq + t*apq;
        A[p][q] = 0;
        A[q][p] = 0;
        // deal with the rotation of the row and column.
        for (int i = 0; i < len; ++i) {
                    if (i != p && i != q) {
                        double aip = A[i][p];
                        double aiq = A[i][q];
                        A[i][p] = aip - s * (aiq + tau * aip);
                        A[i][q] = aiq + s * (aip - tau * aiq);
                        A[p][i] = A[i][p];
                        A[q][i] = A[i][q];
                    }
                }
      }
    }
    
    // if (!significant){
    //  break;
    // }

    // square sum printout as demanded by 1.2 
    double squaresum = 0;
    for (int p = 0; p < len; p++){
      for (int q = p + 1; q < len; q++){
        squaresum += A[p][q]*A[p][q];
      }
    }
    //cout << "square sum at sweep " << sweeps << " = " << squaresum << endl; 
  }

  for (int k = 0; k < len; k++) {
          d[k] = A[k][k];
  }
}

void Hnm(Matrix& H){
  int len = H.size();
  for (int n=1; n <= len; n++){
    for (int m=1; m <= len; m++){
      double val = (n == m ? 1 : 0)*(n*n) + ( 0.05*min(n,m) + 5*pow(-1,abs(m-n))*min(n,m));
      H[n-1][m-1] = val;
    }
  }
}
