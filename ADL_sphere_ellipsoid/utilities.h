#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

static const int N_gridlines = 16;
static const int regularization = 5; //0:no regularization; 3:3rd order; 5:5th order
static const int example = 2; //which example (1 or 2) to run
static const int direct_or_IE = 2; //1:direct evaluation of MADL (example 1), 2:IE (example 1or2)

static const double ellipse_a = 1.0;
static const double ellipse_b = 0.5;
static const double ellipse_c = 0.5;

// constant used in solution of Integral Equation
static const double beta = 0.7;
static const double IE_tol = 1.e-08;
static const int IE_iter_max = 100;


//------------------------------------------------

static const double PI = 3.14159265358979323846;
static const double rootPI = sqrt(PI);
static const double theta = 70.0*PI/180.0; //angle in Beale et al (2.2)
static const double tol = 1e-14; //tolerance in the search of quadrature points (Newton's/bisection method)

struct Surf_point
{
  Surf_point() : x(3,0), f(3,0), g(3,0), vel(3,0), Nrml(3,0),  T1(3,0),  T2(3,0), Area(0) {}
  vector<double> x;
  vector<double> f;
  vector<double> g;
  vector<double> vel;
  vector<double> Nrml;
  vector<double> T1;
  vector<double> T2;
  double Area;
};

struct Target_point
{
Target_point() : flag(0), x(0), y(0), z(0), S1(0), S2(0), S3(0), T1(0), T2(0), T3(0) {}
  int flag;
  double S1, S2, S3;
  double T1, T2, T3;
  double x, y, z;
};

double phi(const vector<double>& x);
double Dphi(int i, const vector<double>& x);
void D2phi(const vector<double>& x, double& phi11, double& phi12, double& phi13,
	   double& phi21, double& phi22, double& phi23, double& phi31,
	   double& phi32, double& phi33);

double Lapl_SL_5order(double r, double d);
double s_Lapl_SL_5order(double r);

double Lapl_DL_5order(double r, double d);
double s_Lapl_DL_5order(double r);

double Lapl_DL_3order(double r, double d);
double s_Lapl_DL_3order(double r);

double distance(const vector<double>& x,
		const vector<double>& y,
		const vector<double>& normal);

double dot_product(const vector<double>& x,
		   const vector<double>& y);

vector<double> cross_product(const vector<double>& x,
			     const vector<double>& y);

void error_vector(int n,
		  const vector<double>& x_comp,
		  const vector<double>& x_ex,
		  double &err_max,
		  int &ind_max,
		  double &err_l2);
