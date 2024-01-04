#include <iostream>
#include <vector>

#include "utilities.h" //header file

using namespace std;

//*********************************************************************

double phi(const vector<double>& x) {
  vector<double> a2(3,0);
  a2[0] = ellipse_a*ellipse_a;
  a2[1] = ellipse_b*ellipse_b;
  a2[2] = ellipse_c*ellipse_c;
  
  return x[0]*x[0]/a2[0] + x[1]*x[1]/a2[1] + x[2]*x[2]/a2[2] - 1.0; 
}

//*********************************************************************

// D_phi/D_x(i): i-th derivative of phi

double Dphi(int i, const vector<double>& x) {
  vector<double> a2(3,0);
  a2[0] = ellipse_a*ellipse_a;
  a2[1] = ellipse_b*ellipse_b;
  a2[2] = ellipse_c*ellipse_c;
  
  return 2.0*x[i]/a2[i]; 
}

//*********************************************************************

// Second derivatives of phi
void D2phi(const vector<double>& x, double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23, double& phi31,
	   double& phi32, double& phi33) {
  
  vector<double> a2(3,0);
  a2[0] = ellipse_a*ellipse_a;
  a2[1] = ellipse_b*ellipse_b;
  a2[2] = ellipse_c*ellipse_c;

  phi11 = 2.0/a2[0]; phi12 = 0.0;       phi13 = 0.0;
  phi21 = 0.0;       phi22 = 2.0/a2[1]; phi23 = 0.0;
  phi31 = 0.0;       phi32 = 0.0;       phi33 = 2.0/a2[2];
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

