#include <iostream>
#include <vector>
#include <cmath>

#include "utilities.h"

using namespace std;

void Generate_Surface(int N, double h, vector<Surf_point>* Surfc);
int sign_phi(const vector<double>& x);
double Find_hypersurf_pt(const vector<double>& A, double H_seg, int i);
double Newton(const vector<double>& PT, double H_seg, int i, double a, double b);
double bisection(const vector<double>& PT, double H_seg, int i, double a, double b);
double Part_Unity(int i, const vector<double>& Nrml);
double b(double r);
double Compute_surface_area(int N_quad, const vector<Surf_point>& Surfc);
double Integrand(const vector<double>& x);

double mean_curvature(const vector<double>& x);
void select_Monge_Patch(int i, const vector<double>& x0, int& i1, int& i2, int& i3);

void generate_Interp_Stencil(double h, double a1, double a2,
			     vector<double>& p1, vector<double>& p2,
			     vector< vector<double> >& M);
void generate_Interp_Stencil_Cubic(double h, double a1, double a2,
				   vector<double>& p1, vector<double>& p2,
				   vector< vector<double> >& M);

void first_Derivatives(double a1, double a2, vector<double>& soln,
		       double& D1, double& D2);
void first_Derivatives_Cubic(double a1, double a2, vector<double>& soln,
			     double& D1, double& D2);
void second_Derivatives(double a1, double a2, vector<double>& soln,
			double& D11, double& D12, double& D22);
void second_Derivatives_Cubic(double a1, double a2, vector<double>& soln,
			      double& D11, double& D12, double& D22);

void find_Third_Coordinate(int i, int n, int i1, int i2, int i3, double z0,
			   const vector<double>& p1,
			   const vector<double>& p2, vector<double>& p3);
void compute_Dual_Tangents(int i1, int i2, int i3, double dhd1, double dhd2,
			   vector<double>& Dual_T1, vector<double>& Dual_T2);

void Generate_Targets_OnSurface(int N_quad, const vector<Surf_point>& Surfc, vector<Target_point>* Target);

