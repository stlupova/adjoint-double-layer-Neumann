#include <iostream>
#include <fstream>
#include <iomanip>
#include <vector>
#include <time.h>

#include "laplace.h"

using namespace std;
    
// direct evaluation of MADL
void eval_MADL_direct(double h, double DEL,
		      int N_quad, vector<Surf_point>& Surfc,
		      int N_target, vector<Target_point>& Target);

// solving integral equation with MADL
void solve_IE_MADL(double h, double DEL,
		   int N_quad, vector<Surf_point>& Surfc,
		   int N_target, vector<Target_point>& Target);


void set_Neumann_exact_soln1(int N, vector<Target_point>& Target,
			     vector<double>& f_ex,
			     vector<double>& u_ex,
			     vector<double>& g_ex,
			     double& u_ex_at0);

void set_Neumann_exact_soln2(int N, vector<Surf_point>& Surfc,
			     vector<Target_point>& Target,
			     vector<double>& u_ex,
			     vector<double>& g_ex,
			     double& u_ex_at0);

//*********************************************************************

int main(int argc, char** argv) {
  int N = N_gridlines;
  double L = 1.0;     // bounding box: [-L,L]^3
  double h = L / N;   // grid spacing
  double DEL = 1.5*pow(h,4.0/5.0); // regularization parameter
  

  vector<Surf_point> Surfc; // an array of surface (quadrature)
                            // points with all their attributes

  Generate_Surface(N, h, &Surfc);
  int N_quad = Surfc.size();
  cout << "Ellipsoid (or sphere) x^2/a^2 + y^2/b^2 + z^2/c^2 = 1; with a=" << ellipse_a
       << ", b=" << ellipse_b << ", c=" << ellipse_c << endl; 
  cout << "Number of quadrature points generated: N = " << N_quad << endl;
  cout << "h = 1/" << N << " = " << h << "     DEL = " << DEL << endl;
    
  clock_t tm = clock();

  // generate target points
  vector<Target_point> Target;
  Generate_Targets_OnSurface(N_quad, Surfc, &Target);

  
  if (direct_or_IE == 1) {
    cout << "Direct evaluation of Modified ADL" << endl;
    eval_MADL_direct(h, DEL,
		     N_quad, Surfc,
		     N_quad, Target);  
  }
  else {
    cout << "Integral equation solution with Modified ADL" << endl;
    solve_IE_MADL(h, DEL,
		  N_quad, Surfc,
		  N_quad, Target);  
  } 

  tm = clock() - tm;
  cout << "CPU time elapsed: " << ((float)tm)/CLOCKS_PER_SEC
       << " seconds" << endl << endl;
  
  return 0;
}

//*********************************************************************

void eval_MADL_direct(double h, double DEL,
		      int N_quad, vector<Surf_point>& Surfc,
		      int N_target, vector<Target_point>& Target) {
  
  vector<double> f_ex(N_quad,0), u_ex(N_quad,0);
  vector<double> g_comp(N_target,0), g_ex(N_target,0);
  double u_ex_at0;

  // set the exact values for boundary data, g_ex
  // and the density function, f_ex
  set_Neumann_exact_soln1(N_target, Target,
			  f_ex,
			  u_ex,
			  g_ex,
			  u_ex_at0);
  
  // compute the modified adjoint double layer integral
  // on the surface, g_comp
  if ( regularization == 0 ) {
    cout << "Modified ADL without regularization" << endl;
    eval_Lapl_MADL_onSurf_no_reg(h, DEL,
				 N_quad, Surfc,
				 f_ex,
				 g_comp);
  }
  else if ( regularization == 3 ) {
    cout << "Modified ADL with 3rd order regularization" << endl;
    eval_Lapl_MADL_onSurf_3order(h, DEL,
				 N_quad, Surfc,
				 f_ex,
				 g_comp);
  }
  else {
    cout << "Modified ADL with 5th order regularization" << endl;
    eval_Lapl_MADL_onSurf_5order(h, DEL,
				 N_quad, Surfc,
				 f_ex,
				 g_comp);
  }
  
  // apply jump condition on the surface
  for (int i=0; i<N_quad; i++)  g_comp[i] -= 0.5*f_ex[i];
  

  // compute and display error
  double err_max, err_l2;
  int ind_max;
  error_vector(N_target, g_comp, g_ex, err_max, ind_max, err_l2);
  cout << "Error between g_comp and g_ex: " << endl;
  cout << "max error = " << err_max << "   l2 error = " << err_l2 << endl;
  
}

//*********************************************************************

void solve_IE_MADL(double h, double DEL,
		   int N_quad, vector<Surf_point>& Surfc,
		   int N_target, vector<Target_point>& Target) {
  
  // set the exact values for solution u_ex, boundary data g_ex
  // and the density function f_ex (in applicable)
  vector<double> f_ex(N_quad,0), g_ex(N_quad,0), u_ex(N_quad,0);
  double u_ex_at0;

  if (example == 1) {
    set_Neumann_exact_soln1(N_target, Target,
			    f_ex,
			    u_ex,
			    g_ex,
			    u_ex_at0);
    
  }
  else {
    set_Neumann_exact_soln2(N_target, Surfc,
			    Target,
			    u_ex,
			    g_ex,
			    u_ex_at0);
  }    

  //compute the surface area so we can adjust f by a constant so that int(fdS)=0
  double surfArea = Compute_surface_area(N_quad, Surfc);

  vector<double> f_old(N_quad,0), f_new(N_quad,0);
  for (int i=0; i<N_quad; i++)  f_old[i] = 0.0;

  
  double f_diff_max = 1.0;
  int n = 0;
  while ( (n < IE_iter_max) && (f_diff_max > IE_tol) ) {
    
    // compute the modified adjoint double layer integral
    vector<double> SOL_comp(N_target,0);
    eval_Lapl_MADL_onSurf_5order(h, DEL,
				 N_quad, Surfc,
				 f_old,
				 SOL_comp);
  
    for (int i=0; i<N_quad; i++) {
      f_new[i] = (1.0 - beta) * f_old[i] + 2.0 * beta * (SOL_comp[i] - g_ex[i]);
    }


    // start modify f
    // find and add constant a to f, ff=f+a, so that
    // int(ff*dS) = int((f+a)*dS) = int(f*dS) + a*int(dS) = int(f*dS) + a*Surface_area = 0
    // the constant then is a = -int(f*dS) / Surface_area
    double int_fdS = 0.0;
    for (int i=0; i<N_quad; i++) {
      int_fdS += f_new[i] * Surfc[i].Area; // first compute int(f*dS)
    }
    double a_const = -int_fdS / surfArea; // second, compute constant a
    
    cout << "Iteration " << n+1 << "   a = " << a_const;
    
    for (int i=0; i<N_quad; i++) {
      f_new[i] += a_const;  // add constant a to f
    }
    // end modify f


    // compute difference between f_new and f_old ---are we meeting the set tolerance?
    double f_diff_l2;
    int f_ind_max = 0;
    error_vector(N_target, f_new, f_old, f_diff_max, f_ind_max, f_diff_l2);

    cout << "  f_diff_max = " << f_diff_max << endl; 

    for (int i=0; i<N_quad; i++) {
      f_old[i] = f_new[i];
    }
    n++; 
  }

  cout << endl;
  cout << "Computing the single layer potential" << endl;
  vector<double> u_comp(N_quad,0), u_comp0(N_quad,0);

  eval_Lapl_SL(h, DEL,
	       N_quad, Surfc,
	       f_new,
	       u_comp);

  
  // start adjust constant in u
  // compute u(0) numerical value
  double u_at0;
  eval_Lapl_SL_atOrigin(h, DEL,
			N_quad, Surfc,
			f_new,
			u_at0);
  
  cout << "u_comp(0) = " << u_at0 << "    u_ex(0) = " << u_ex_at0 << endl;

  // adjust u = u-u(0)+u_ex(0)
  // so that u(0) = u_ex(0) ---same constant in "up to a constant"
  for (int i=0; i<N_quad; i++) u_comp[i] += u_ex_at0 - u_at0;

  // end adjust constant in u
  
  
  double err_max, err_l2;
  int ind_max;

  if (example == 1) {
    error_vector(N_target, f_new, f_ex, err_max, ind_max, err_l2);
    cout << "Error in f:" << endl;
    cout << "max error = " << err_max << "  l2 error = " << err_l2 << endl << endl;
  }
  
  error_vector(N_target, u_comp, u_ex, err_max, ind_max, err_l2);
  cout << "Error in u as SL with f_comp:" << endl;
  cout << "max error = " << err_max << "  l2 error = " << err_l2 << endl << endl;


  // (optional) check accuracy in u if f_ex is used --- only applies to example 1
  if (example == 1) {
    eval_Lapl_SL(h, DEL,
		 N_quad, Surfc,
		 f_ex,
		 u_comp0);
    
    error_vector(N_target, u_comp0, u_ex, err_max, ind_max, err_l2);
    cout << "Error in u as SL with f_ex:" << endl;
    cout << "max error = " << err_max << "  l2 error = " << err_l2 << endl; 
  }
  
}

//*********************************************************************

// example 1
// exact solution similar to example on p.617 of Beale (2004)
// but for the Neumann problem
// f_ex is the density function based on a spherical harmonic
// g_ex is the Neumann boundary condition du/dn = g

void set_Neumann_exact_soln1(int N, vector<Target_point>& Target,
			     vector<double>& f_ex,
			     vector<double>& u_ex,
			     vector<double>& g_ex,
			     double& u_ex_at0) {
  
  double m11 = sqrt(2.0/6.0), m21 = m11, m31 = m11;
  double m12 = 0.0, m22 = sqrt(3.0/6.0), m32 = -m22;
  double m13 = -2.0/sqrt(6.0), m23 = 1.0/sqrt(6.0), m33 = m23;
  
  for (int i=0; i<N; i++) {
    
    double x1 = Target[i].x;
    double x2 = Target[i].y;
    double x3 = Target[i].z;

    double y1 = m11 * x1 + m12 * x2 + m13 * x3;
    double y2 = m21 * x1 + m22 * x2 + m23 * x3;
    double y3 = m31 * x1 + m32 * x2 + m33 * x3;

    f_ex[i] = 1.75  *  (y1 - 2.0 * y2)  *  (7.5 * y3 * y3 - 1.5);

    u_ex[i] = -1.0 / 7.0 * f_ex[i];

    g_ex[i] = -3.0 / 7.0 * f_ex[i];
  }
  u_ex_at0 = 0.0;
}

//*********************************************************************

// example 2
// exact solution similar to example on p.618 of Beale (2004)
// but for the Neumann problem
// u_ex is the harmonic function exp(y1+2y2)*cos(sqrt(5)y3)
// g_ex is the Neumann boundary condition du/dn = g

void set_Neumann_exact_soln2(int N, vector<Surf_point>& Surfc,
			     vector<Target_point>& Target,
			     vector<double>& u_ex,
			     vector<double>& g_ex,
			     double& u_ex_at0)  {

  double m11 = sqrt(2.0/6.0), m21 = m11, m31 = m11;
  double m12 = 0.0, m22 = sqrt(3.0/6.0), m32 = -m22;
  double m13 = -2.0/sqrt(6.0), m23 = 1.0/sqrt(6.0), m33 = m23;

  for (int i=0; i<N; i++) {
    
    double x1 = Target[i].x;
    double x2 = Target[i].y;
    double x3 = Target[i].z;

    double y1 = m11 * x1 + m12 * x2 + m13 * x3;
    double y2 = m21 * x1 + m22 * x2 + m23 * x3;
    double y3 = m31 * x1 + m32 * x2 + m33 * x3;

    double var1 = exp(y1+2.0*y2);
    double var2 = cos(sqrt(5.0)*y3);
    
    u_ex[i] = var1 * var2;
    
    double DuDy1 = var1 * var2;
    double DuDy2 = 2.0 * var1 * var2;
    double DuDy3 = -sqrt(5.0) * var1 * sin(sqrt(5.0) * y3);

    vector<double> DuDx(3,0);
    DuDx[0] = m11 * DuDy1 + m21 * DuDy2 + m31 * DuDy3;
    DuDx[1] = m12 * DuDy1 + m22 * DuDy2 + m32 * DuDy3;
    DuDx[2] = m13 * DuDy1 + m23 * DuDy2 + m33 * DuDy3;

    g_ex[i] = dot_product(DuDx, Surfc[i].Nrml);
  }
  u_ex_at0 = 1.0;
}
