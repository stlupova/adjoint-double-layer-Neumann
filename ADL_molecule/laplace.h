#include <iostream>
#include <vector>
#include <cmath>

#include "surface.h"

using namespace std;

void eval_Lapl_SL(double h, double DEL,
		  int N_quad, const vector<Surf_point>& Surfc,
		  const vector<double>& g,
		  vector<double>& SL);

void eval_Lapl_SL_atOrigin(double h, double DEL,
			   int N_quad, const vector<Surf_point>& Surfc,
			   const vector<double>& g,
			   double& SL);

void eval_Lapl_DL_onSurf(double h, double DEL,
			 int N_quad, const vector<Surf_point>& Surfc,
			 const vector<double>& g,
			 vector<double>& DL);

void eval_Lapl_DL_onSurf_NoSubtraction(double h, double DEL,
				       int N_quad, const vector<Surf_point>& Surfc,
				       const vector<double>& g,
				       vector<double>& DL);

void eval_Lapl_ADL_onSurf(double h, double DEL,
			  int N_quad, const vector<Surf_point>& Surfc,
			  const vector<double>& g,
			  vector<double>& DL);


void eval_Lapl_MADL_onSurf_5order(double h, double DEL,
				  int N_quad, const vector<Surf_point>& Surfc,
				  const vector<double>& g,
				  vector<double>& DL);

void eval_Lapl_MADL_onSurf_3order(double h, double DEL,
				  int N_quad, const vector<Surf_point>& Surfc,
				  const vector<double>& g,
				  vector<double>& DL);

void eval_Lapl_MADL_onSurf_no_reg(double h, double DEL,
				  int N_quad, const vector<Surf_point>& Surfc,
				  const vector<double>& g,
				  vector<double>& DL);

