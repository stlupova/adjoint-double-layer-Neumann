#include <iostream>
#include <vector>

#include "utilities.h" //header file

using namespace std;

//*********************************************************************

double phi(const vector<double>& x) {

  vector<double> x1(3,0), x2(3,0), x3(3,0), x4(3,0);
  x1[0] = sqrt(3.0)/3.0;  x1[1] = 0.0;  x1[2] =-sqrt(6.0)/12.0;
  x2[0] =-sqrt(3.0)/6.0;  x2[1] = 0.5;  x2[2] =-sqrt(6.0)/12.0;
  x3[0] =-sqrt(3.0)/6.0;  x3[1] =-0.5;  x3[2] =-sqrt(6.0)/12.0;
  x4[0] = 0.0;            x4[1] = 0.0;  x4[2] = sqrt(6.0)/4.0;
  
  double r2 = 0.5*0.5;
  double c = 0.6;

  double a1, a2, a3, a4;
  a1 = (x[0]-x1[0])*(x[0]-x1[0]) + (x[1]-x1[1])*(x[1]-x1[1]) + (x[2]-x1[2])*(x[2]-x1[2]); 
  a2 = (x[0]-x2[0])*(x[0]-x2[0]) + (x[1]-x2[1])*(x[1]-x2[1]) + (x[2]-x2[2])*(x[2]-x2[2]); 
  a3 = (x[0]-x3[0])*(x[0]-x3[0]) + (x[1]-x3[1])*(x[1]-x3[1]) + (x[2]-x3[2])*(x[2]-x3[2]); 
  a4 = (x[0]-x4[0])*(x[0]-x4[0]) + (x[1]-x4[1])*(x[1]-x4[1]) + (x[2]-x4[2])*(x[2]-x4[2]); 

  double val = c - exp(-a1/r2) - exp(-a2/r2) - exp(-a3/r2) - exp(-a4/r2);
  return val;
}

//*********************************************************************

// D_phi/D_x(i): i-th derivative of phi

double Dphi(int i, const vector<double>& x) {

  vector<double> x1(3,0), x2(3,0), x3(3,0), x4(3,0);
  x1[0] = sqrt(3.0)/3.0;  x1[1] = 0.0;  x1[2] =-sqrt(6.0)/12.0;
  x2[0] =-sqrt(3.0)/6.0;  x2[1] = 0.5;  x2[2] =-sqrt(6.0)/12.0;
  x3[0] =-sqrt(3.0)/6.0;  x3[1] =-0.5;  x3[2] =-sqrt(6.0)/12.0;
  x4[0] = 0.0;            x4[1] = 0.0;  x4[2] = sqrt(6.0)/4.0;
  
  double r2 = 0.5*0.5;

  double a1, a2, a3, a4;
  a1 = (x[0]-x1[0])*(x[0]-x1[0]) + (x[1]-x1[1])*(x[1]-x1[1]) + (x[2]-x1[2])*(x[2]-x1[2]); 
  a2 = (x[0]-x2[0])*(x[0]-x2[0]) + (x[1]-x2[1])*(x[1]-x2[1]) + (x[2]-x2[2])*(x[2]-x2[2]); 
  a3 = (x[0]-x3[0])*(x[0]-x3[0]) + (x[1]-x3[1])*(x[1]-x3[1]) + (x[2]-x3[2])*(x[2]-x3[2]); 
  a4 = (x[0]-x4[0])*(x[0]-x4[0]) + (x[1]-x4[1])*(x[1]-x4[1]) + (x[2]-x4[2])*(x[2]-x4[2]); 

  double d1, d2, d3, d4;
  d1 = x[i] - x1[i];
  d2 = x[i] - x2[i];
  d3 = x[i] - x3[i];
  d4 = x[i] - x4[i];

  double val = 2.0/r2* (exp(-a1/r2)*d1 + exp(-a2/r2)*d2 + exp(-a3/r2)*d3 + exp(-a4/r2)*d4);
  return val;
}

//*********************************************************************

// Second derivatives of phi
void D2phi(const vector<double>& x, double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23, double& phi31,
	   double& phi32, double& phi33) {
  
  vector<double> x1(3,0), x2(3,0), x3(3,0), x4(3,0);
  x1[0] = sqrt(3.0)/3.0;  x1[1] = 0.0;  x1[2] =-sqrt(6.0)/12.0;
  x2[0] =-sqrt(3.0)/6.0;  x2[1] = 0.5;  x2[2] =-sqrt(6.0)/12.0;
  x3[0] =-sqrt(3.0)/6.0;  x3[1] =-0.5;  x3[2] =-sqrt(6.0)/12.0;
  x4[0] = 0.0;            x4[1] = 0.0;  x4[2] = sqrt(6.0)/4.0;
  
  double r2 = 0.5*0.5;

  vector<double> d1 = x; d1[0]-=x1[0]; d1[1]-=x1[1]; d1[2]-=x1[2];
  vector<double> d2 = x; d2[0]-=x2[0]; d2[1]-=x2[1]; d2[2]-=x2[2];
  vector<double> d3 = x; d3[0]-=x3[0]; d3[1]-=x3[1]; d3[2]-=x3[2];
  vector<double> d4 = x; d4[0]-=x4[0]; d4[1]-=x4[1]; d4[2]-=x4[2];
  double a1, a2, a3, a4;
  a1 = exp(-dot_product(d1,d1)/r2);
  a2 = exp(-dot_product(d2,d2)/r2);
  a3 = exp(-dot_product(d3,d3)/r2); 
  a4 = exp(-dot_product(d4,d4)/r2); 

  phi11 = 2.0/r2*(a1*(1.0-2.0/r2*d1[0]*d1[0]) + a2*(1.0-2.0/r2*d2[0]*d2[0]) +
		  a3*(1.0-2.0/r2*d3[0]*d3[0]) + a4*(1.0-2.0/r2*d4[0]*d4[0]));
  phi22 = 2.0/r2*(a1*(1.0-2.0/r2*d1[1]*d1[1]) + a2*(1.0-2.0/r2*d2[1]*d2[1]) +
		  a3*(1.0-2.0/r2*d3[1]*d3[1]) + a4*(1.0-2.0/r2*d4[1]*d4[1]));
  phi33 = 2.0/r2*(a1*(1.0-2.0/r2*d1[2]*d1[2]) + a2*(1.0-2.0/r2*d2[2]*d2[2]) +
		  a3*(1.0-2.0/r2*d3[2]*d3[2]) + a4*(1.0-2.0/r2*d4[2]*d4[2]));
  
  phi12 = -4.0/(r2*r2)*(a1*d1[0]*d1[1] + a2*d2[0]*d2[1] +
			a3*d3[0]*d3[1] + a4*d4[0]*d4[1]);
  phi13 = -4.0/(r2*r2)*(a1*d1[0]*d1[2] + a2*d2[0]*d2[2] +
			a3*d3[0]*d3[2] + a4*d4[0]*d4[2]);
  phi23 = -4.0/(r2*r2)*(a1*d1[1]*d1[2] + a2*d2[1]*d2[2] +
			a3*d3[1]*d3[2] + a4*d4[1]*d4[2]);

  phi21 = phi12;
  phi31 = phi13;
  phi32 = phi23;
}

//*********************************************************************
//*********************************************************************

double Lapl_SL_5order(double r, double d) {
  if (r < 1e-14) {
    return 16.0/3.0/rootPI/d;
  }
  else {
    return s_Lapl_SL_5order(r/d)/r;
  }
}

//*********************************************************************

double s_Lapl_SL_5order(double r) {
  double r2 = r*r;
  return erf(r) - 2.0*(2.0*r2 - 5.0)*r*exp(-r2)/3.0/rootPI;
}

//*********************************************************************

double Lapl_DL_5order(double r, double d) {
  if (r < 1e-14) {
    return 8.0/3.0/rootPI/(d*d*d);
  }
  else {
    return s_Lapl_DL_5order(r/d)/(r*r*r);
  }
}

//*********************************************************************

double s_Lapl_DL_5order(double r) {
  double r2 = r*r;
  return erf(r) + 2.0*(2.0*r2 - 3.0)*r*exp(-r2)/3.0/rootPI;
}

//*********************************************************************

double s_Lapl_DL_3order(double r) {
  double r2 = r*r;
  return erf(r) - (2.0 / rootPI)*r*exp(-r2);
}

//*********************************************************************
double Lapl_DL_3order(double r, double d) {
  if (r < 1e-14) {
    return 8.0/3.0/rootPI/(d*d*d);
  }
  else {
    return s_Lapl_DL_3order(r/d)/(r*r*r);
  }
}

//*********************************************************************
//*********************************************************************

// Signed distance between two points in 3D

double distance(const vector<double>& x,
		const vector<double>& y,
		const vector<double>& normal) {
  
  return (x[0]-y[0])*normal[0] + (x[1]-y[1])*normal[1] + (x[2]-y[2])*normal[2];
}

//*********************************************************************

double dot_product(const vector<double>& x,
		   const vector<double>& y) {
  
  return x[0]*y[0] + x[1]*y[1] + x[2]*y[2];
}

//*********************************************************************

vector<double> cross_product(const vector<double>& x,
			     const vector<double>& y) {
  
  vector<double> cross_prod(3,0);
  cross_prod[0] = x[1]*y[2] - x[2]*y[1];
  cross_prod[1] = x[2]*y[0] - x[0]*y[2];
  cross_prod[2] = x[0]*y[1] - x[1]*y[0];
  return cross_prod;
}

//*********************************************************************
//*********************************************************************

void error_vector(int n,
		  const vector<double>& x_comp,
		  const vector<double>& x_ex,
		  double &err_max,
		  int &ind_max,
		  double &err_l2) {

  err_max = abs(x_comp[0] - x_ex[0]);
  ind_max = 0;
  err_l2 = err_max * err_max;
  
  for (int i=1; i<n; i++) {
    double err_temp = abs(x_comp[i] - x_ex[i]);
    if ( err_temp > err_max) {
      err_max = err_temp;
      ind_max = i;
    }
    err_l2 += err_temp * err_temp;
  }
  err_l2 = sqrt(err_l2 / n);

}

