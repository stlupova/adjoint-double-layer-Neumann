#include <iostream>
#include <vector>
#include <iomanip>

#include "laplace.h" 

using namespace std;

//*********************************************************************
// evaluates laplace single layer potential without forming the matrix
// high order on-surface regularization is used

void eval_Lapl_SL(double h, double DEL,
		  int N_quad, const vector<Surf_point>& Surfc,
		  const vector<double>& g,
		  vector<double>& SL) {
  
  vector<double> pt(3,0), dx(3,0);
  double FourPI = -1.0/(4.0*PI);

  for (int i=0; i<N_quad; i++)  SL[i] = 0.0;
  
  for (int i=0; i<N_quad; i++) {
    pt[0] = Surfc[i].x[0];
    pt[1] = Surfc[i].x[1];
    pt[2] = Surfc[i].x[2];
    
    for (int j=0; j<N_quad; j++) {
      for (int k=0; k<3; k++)  dx[k] = pt[k]-Surfc[j].x[k];
      
      double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      
      SL[i] += g[j] * Lapl_SL_5order(r,DEL) * Surfc[j].Area;      
    }
    SL[i] *= FourPI;
  }
}

//*********************************************************************
// evaluates laplace single layer potential at the origin
// since origin is away from boundary, regularization is not needed

void eval_Lapl_SL_atOrigin(double h, double DEL,
			   int N_quad, const vector<Surf_point>& Surfc,
			   const vector<double>& g,
			   double& SL) {
  
  vector<double> pt(3,0), dx(3,0);
  double FourPI = -1.0/(4.0*PI);

  SL = 0.0;
  
  for (int j=0; j<N_quad; j++) {
    for (int k=0; k<3; k++)  dx[k] = 0.0-Surfc[j].x[k];    

    double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);

    SL += g[j] * Surfc[j].Area / r;      
  }
  SL *= FourPI;
}

//*********************************************************************
// evaluates laplace double layer potential without forming the matrix
// subtraction is used:
// H*phi = int(dG/dn * phi) = int(dG/dn * (phi-phi0)) + 0.5*phi0
// high order on-surface regularization is used

void eval_Lapl_DL_onSurf(double h, double DEL,
			 int N_quad, const vector<Surf_point>& Surfc,
			 const vector<double>& g,
			 vector<double>& DL) {
  
  vector<double> pt(3,0), dx(3,0);
  double FourPI = -1.0/(4.0*PI);
  
  for (int i=0; i<N_quad; i++)  DL[i] = 0.0;
  
  for (int i=0; i<N_quad; i++) {
    
    pt[0] = Surfc[i].x[0];
    pt[1] = Surfc[i].x[1];
    pt[2] = Surfc[i].x[2];
    
    for (int j=0; j<N_quad; j++) {
      
      for (int k=0; k<3; k++)  dx[k] = pt[k]-Surfc[j].x[k];
      
      double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      double n_dot_dx = dot_product(Surfc[j].Nrml,dx);
      
      DL[i] += (g[j]-g[i]) * n_dot_dx * Lapl_DL_5order(r,DEL) * Surfc[j].Area;
      
    }
    DL[i] = DL[i] * FourPI + 0.5*g[i];
  }
}

//*********************************************************************
// evaluates laplace double layer potential without forming the matrix

// H*phi = int(dG/dn * phi) = int(dG/dn * (phi-phi0)) + 0.5*phi0
// high order on-surface regularization is used

void eval_Lapl_DL_onSurf_NoSubtraction(double h, double DEL,
				       int N_quad, const vector<Surf_point>& Surfc,
				       const vector<double>& g,
				       vector<double>& DL) {
  
  vector<double> pt(3,0), dx(3,0);
  double FourPI = -1.0/(4.0*PI);
  
  for (int i=0; i<N_quad; i++)  DL[i] = 0.0;
  
  for (int i=0; i<N_quad; i++) {
    
    pt[0] = Surfc[i].x[0];
    pt[1] = Surfc[i].x[1];
    pt[2] = Surfc[i].x[2];
    
    for (int j=0; j<N_quad; j++) {
      
      for (int k=0; k<3; k++)  dx[k] = pt[k]-Surfc[j].x[k];
      
      double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      double n_dot_dx = dot_product(Surfc[j].Nrml,dx);
      
      DL[i] += g[j] * n_dot_dx * Lapl_DL_5order(r,DEL) * Surfc[j].Area;
      
    }
    DL[i] = DL[i] * FourPI;//  +  0.5 * g[i];
  }
}

//*********************************************************************
// evaluates laplace adjoint double layer potential without forming the matrix
// high order on-surface regularization is used

void eval_Lapl_ADL_onSurf(double h, double DEL,
				       int N_quad, const vector<Surf_point>& Surfc,
				       const vector<double>& g,
				       vector<double>& DL) {
  
  vector<double> pt(3,0), dx(3,0);
  double FourPI = 1.0/(4.0*PI);
  
  for (int i=0; i<N_quad; i++)  DL[i] = 0.0;
  
  for (int i=0; i<N_quad; i++) {
    
    pt[0] = Surfc[i].x[0];
    pt[1] = Surfc[i].x[1];
    pt[2] = Surfc[i].x[2];
    
    for (int j=0; j<N_quad; j++) {
      
      for (int k=0; k<3; k++)  dx[k] = pt[k]-Surfc[j].x[k];
      
      double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      double n_dot_dx = dot_product(Surfc[i].Nrml,dx);
      
      DL[i] += g[j] * n_dot_dx * Lapl_DL_5order(r,DEL) * Surfc[j].Area;
      
    }
    DL[i] = DL[i] * FourPI;// + 0.5 * g[i];
  }
}

//*********************************************************************
// evaluates modified laplace adjoint double layer potential on surface
// 5th order regularization

void eval_Lapl_MADL_onSurf_5order(double h, double DEL,
				  int N_quad, const vector<Surf_point>& Surfc,
				  const vector<double>& g,
				  vector<double>& ADL) {
  
  vector<double> pt(3,0), dx(3,0);
  double FourPI = 1.0/(4.0*PI);
  
  for (int i=0; i<N_quad; i++)  ADL[i] = 0.0;
  
  for (int i=0; i<N_quad; i++) {
    
    pt[0] = Surfc[i].x[0];
    pt[1] = Surfc[i].x[1];
    pt[2] = Surfc[i].x[2];
    
    for (int j=0; j<N_quad; j++) {
      
      for (int k=0; k<3; k++)  dx[k] = pt[k]-Surfc[j].x[k];
      
      double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      double modif = g[j] * dot_product(Surfc[i].Nrml,dx)
	           + g[i] * dot_product(Surfc[j].Nrml,dx);
      
      ADL[i] += modif * Lapl_DL_5order(r,DEL) * Surfc[j].Area;
      
    }
    ADL[i] = ADL[i] * FourPI + 0.5*g[i];
  }
}

//*********************************************************************
// evaluates modified laplace adjoint double layer potential on surface
// 3rd order regularization

void eval_Lapl_MADL_onSurf_3order(double h, double DEL,
				  int N_quad, const vector<Surf_point>& Surfc,
				  const vector<double>& g,
				  vector<double>& ADL) {
  
  vector<double> pt(3,0), dx(3,0);
  double FourPI = 1.0/(4.0*PI);
  
  for (int i=0; i<N_quad; i++)  ADL[i] = 0.0;
  
  for (int i=0; i<N_quad; i++) {
    
    pt[0] = Surfc[i].x[0];
    pt[1] = Surfc[i].x[1];
    pt[2] = Surfc[i].x[2];
    
    for (int j=0; j<N_quad; j++) {
      
      for (int k=0; k<3; k++)  dx[k] = pt[k]-Surfc[j].x[k];
      
      double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
      double modif = g[j] * dot_product(Surfc[i].Nrml,dx)
	           + g[i] * dot_product(Surfc[j].Nrml,dx);
      
      ADL[i] += modif * Lapl_DL_3order(r,DEL) * Surfc[j].Area;
      
    }
    ADL[i] = ADL[i] * FourPI + 0.5*g[i];
  }
}

//*********************************************************************
// evaluates modified laplace adjoint double layer potential on surface
// without regularization
// the point x=y is omitted

void eval_Lapl_MADL_onSurf_no_reg(double h, double DEL,
				  int N_quad, const vector<Surf_point>& Surfc,
				  const vector<double>& g,
				  vector<double>& ADL) {
  
  vector<double> pt(3,0), dx(3,0);
  double FourPI = 1.0/(4.0*PI);
  
  for (int i=0; i<N_quad; i++)  ADL[i] = 0.0;
  
  for (int i=0; i<N_quad; i++) {
    
    pt[0] = Surfc[i].x[0];
    pt[1] = Surfc[i].x[1];
    pt[2] = Surfc[i].x[2];
    
    for (int j=0; j<N_quad; j++) {

      if (i != j) {
	for (int k=0; k<3; k++)  dx[k] = pt[k]-Surfc[j].x[k];
	
	double r = sqrt(dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2]);
	double modif = g[j] * dot_product(Surfc[i].Nrml,dx)
	             + g[i] * dot_product(Surfc[j].Nrml,dx);
	
	ADL[i] += modif * Surfc[j].Area / (r * r * r);
      }
    }
    ADL[i] = ADL[i] * FourPI + 0.5*g[i];
  }
}
