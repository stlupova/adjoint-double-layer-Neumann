#include <iostream>
#include <vector>

#include "surface.h" //header file containing constants, structure and function definitions

using namespace std;

//*********************************************************************

// Find surface (quadrature) points using the method of coordinate patches

void Generate_Coordinate_Patches(int N, double h, vector<Surf_point>* Surfc) {
  double rad = 1.25;
  double ellipse_c = 1.0;
  int N1 = 2.0*rad/h + 1;
  double a1, a2, a_sq, A;
  double phi1, phi2, a1_1, a2_1, a_sq_1, a1_2, a2_2, a_sq_2, part_unity;
  vector<double> pt(3,0); // second argument puts 0 as every value

  
  for (int i = 0; i < N1; i++) {
    for (int j = 0; j < N1; j++) {
      a1 = -rad + i*h;
      a2 = -rad + j*h;
      a_sq = a1*a1 + a2*a2;

      // first coordinate system
      pt[0] = 2.0*a1/(1.0+a_sq);
      pt[1] = 2.0*a2/(1.0+a_sq);
      pt[2] = ellipse_c*(1.0-a_sq)/(1.0+a_sq);
      
      Surf_point temp;
      temp.x = pt;
      temp.Nrml = pt;

      part_unity = 0.0;
      if (a_sq < rad*rad) {
	phi1 = 0.0;
	phi2 = 0.0;
	// the first stereographic projection
        a1_1 = pt[0]/(1.0+pt[2]/ellipse_c);
        a2_1 = pt[1]/(1.0+pt[2]/ellipse_c);
	a_sq_1 = a1_1*a1_1 + a2_1*a2_1;
	if (a_sq_1 < rad*rad) phi1 = exp(-rad*rad/(rad*rad-a_sq_1));
	// the second stereographic projection
        a1_2 = pt[0]/(1.0-pt[2]/ellipse_c);
        a2_2 = pt[1]/(1.0-pt[2]/ellipse_c);
	a_sq_2 = a1_2*a1_2 + a2_2*a2_2;
	if (a_sq_2 < rad*rad) phi2 = exp(-rad*rad/(rad*rad-a_sq_2));

	part_unity = phi1/(phi1+phi2);
      }
      A = 4.0/pow(1.0+a_sq,2)*sqrt(1.0+4.0*(pow(ellipse_c,2)-1)*a_sq/pow(1.0+a_sq,2));
      temp.Area = A*part_unity*h*h;
      Surfc->push_back(temp); //add the point to the quadrature list
      
      // second coordinate system
      pt[2] = -pt[2];
      temp.x = pt;
      temp.Nrml = pt;

      part_unity = 0.0;
      if (a_sq < pow(rad,2)) {
	phi1 = 0.0;
	phi2 = 0.0;
	// the first stereographic projection
        a1_1 = pt[0]/(1.0+pt[2]/ellipse_c);
        a2_1 = pt[1]/(1.0+pt[2]/ellipse_c);
	a_sq_1 = pow(a1_1,2) + pow(a2_1,2);
	if (a_sq_1 < pow(rad,2)) phi1 = exp(-pow(rad,2)/(pow(rad,2)-a_sq_1));
	// the second stereographic projection
        a1_2 = pt[0]/(1.0-pt[2]/ellipse_c);
        a2_2 = pt[1]/(1.0-pt[2]/ellipse_c);
	a_sq_2 = pow(a1_2,2) + pow(a2_2,2);
	if (a_sq_2 < pow(rad,2)) phi2 = exp(-pow(rad,2)/(pow(rad,2)-a_sq_2));

	part_unity = phi2/(phi1+phi2);
      }
      temp.Area = A*part_unity*h*h;
      Surfc->push_back(temp); //add the point to the quadrature list

    }                
  }
}

//*********************************************************************

// Find surface (quadrature) points using the method of Beale et al 2017
//      Uses the simple bracketing algorithm (Wilson's thesis, sec. 2.2.1)

void Generate_Surface(int N, double h, vector<Surf_point>* Surfc) {
  double H_seg = h; // segment spacing
  vector<double> A(3,0), B(3,0); // second argument puts 0 as every value
  int count = 0;
  
  for (int i = 0; i < 3; i++) {
    for (int l1 = -N; l1 <= N; l1++) {
      for (int l2 = -N; l2 <= N; l2++) {
          int ind1,ind2;
	if (i==0) {
        ind1=1; ind2=2;
	}
	if (i==1) {
        ind1=0; ind2=2;
	}
	if (i==2) {
        ind1=0; ind2=1;
	}
    A[ind1] = l1*h; A[ind2] = l2*h;
    B[ind1] = l1*h; B[ind2] = l2*h;
	for (int j = -N; j <= N; j++) { //j-th segment: x_i in [jH,(j+1)H]
	  A[i] = j*H_seg;      //end points of segment
	  B[i] = A[i] + H_seg;

	  if ( sign_phi(A)*sign_phi(B) < 0 ) { //check that S(lambda,i,j) is a bracket
	    
	    double tau = Find_hypersurf_pt(A,H_seg,i); //find a hypersurface point in S
	    Surf_point temp;
	    temp.x = A;
	    temp.x[i] += tau * H_seg; //the point then has i-th coordinate = (j+tau)*H
	    
	    for (int k=0; k<3; k++) temp.Nrml[k] = Dphi(k,temp.x); //compute the normal vector

	    //first tangent
	    temp.T1[ind1]=1.0;
	    temp.T1[ind2]=0.0;
	    temp.T1[i] = -temp.Nrml[ind1]/temp.Nrml[i];

	    //second tangent
	    temp.T2[ind1]=0.0;
	    temp.T2[ind2]=1.0;
	    temp.T2[i] = -temp.Nrml[ind2]/temp.Nrml[i]; 
          
	    double Norm_Dphi = sqrt(temp.Nrml[0]*temp.Nrml[0]
				  + temp.Nrml[1]*temp.Nrml[1]
				  + temp.Nrml[2]*temp.Nrml[2]);
	    
	    //make it a unit normal
	    for (int k=0; k<3; k++) temp.Nrml[k] /= Norm_Dphi; 

	    if ( abs(temp.Nrml[i]) >= cos(theta) ) { //see Beale et al (2.4)
	      temp.Area = h*h*Part_Unity(i,temp.Nrml)/abs(temp.Nrml[i]); //"area element"	      
	      Surfc->push_back(temp); //add the point to the quadrature list
	      count = count + 1;
	    }
	  }                
	}
      }
    }
  }
}

//*********************************************************************

// Sign function as defined in Wilson's thesis p.12

int sign_phi(const vector<double>& x) {
  if ( phi(x) < 0.0 )
      return -1;
  else
    return 1;
}

//*********************************************************************

// Find a quadrature point by a line search algorithm, essentially
// looking for a root of phi on a given interval using Newton's/bisection method

double Find_hypersurf_pt(const vector<double>& A, double H_seg, int i) {
  //return bisection(A,H_seg,i,0.0,1.0);
  return Newton(A,H_seg,i,0,1);
}

//*********************************************************************

double Newton(const vector<double>& PT, double H_seg, int i, double a, double b) {

    int n_max = 100;
    vector<double> x_left = PT;
    vector<double> x_right = PT;
    x_left[i] += a*H_seg;
    x_right[i] += b*H_seg;
    vector<double> x_n = PT;
    double phi_n, Dphi_n;
    int n = 0;
    double p = 0.0;

    if ( phi(x_left)*phi(x_right) > 0.0 ) 
      cout << "Root finding will fail: endpoints do not have opposite sign" << endl;
    else if ( abs(phi(x_left)) <= tol )
        p = a;
    else if ( abs(phi(x_right)) <= tol )
        p = b;
    else {
      // start Newton's iteration using the midpoint as initial guess
      p = (a + b)/2;
      x_n[i] = PT[i] + p*H_seg;
      phi_n = phi(x_n);
      Dphi_n = Dphi(i,x_n)*H_seg;
      while ((abs(phi_n) > tol) && (abs(Dphi_n) > tol) && (n < n_max)) {
	p -= phi_n/Dphi_n; 
	x_n[i] = PT[i] + p*H_seg;
	n = n+1;
	phi_n = phi(x_n);
	Dphi_n = Dphi(i,x_n)*H_seg;
      }
      // if Newton's failed to converge, restart with the bisection method
      if ((n >= n_max) || (abs(Dphi_n) <= tol)) {
	cout << "Newton's didn't converge... restart with bisection..." << endl;
	p = bisection(PT,H_seg,i,a,b);
      }
    }
    return p;
}

//*********************************************************************

double bisection(const vector<double>& PT, double H_seg, int i, double a, double b) {

    vector<double> x_left = PT;
    vector<double> x_right = PT;
    x_left[i] += a*H_seg;
    x_right[i] += b*H_seg;
    vector<double> x_mid = PT;
    double phi_mid;
    double p = 0.0;

    if ( phi(x_left)*phi(x_right) > 0.0 ) 
      cout << "Bisection will fail: endpoints do not have opposite sign" << endl;
    else if ( abs(phi(x_left)) <= tol )
      p = a;
    else if ( abs(phi(x_right)) <= tol )
      p = b;
    else {
      p = (a + b)/2;
      x_mid[i] = PT[i] + p*H_seg;
      phi_mid = phi(x_mid);
      while ( abs(phi_mid) > tol ) {
	if ( phi(x_left)*phi(x_mid) < 0.0 ) 
	  b = p;
	else {
	  a = p;     
	  x_left = x_mid;
	}
	p = (a + b)/2; 
	x_mid[i] = PT[i] + p*H_seg;
	phi_mid = phi(x_mid);
      }
    }
    return p;
}

//*********************************************************************

// Computes the sigma functions defined in Beale et al p.5.
//          These are universal partitions of unity on the unit sphere,
//          but applied to the unit normal of the given surface

double Part_Unity(int i, const vector<double>& Nrml) {
  double n1 = abs(Nrml[0]); // abs(n*e1)
  double n2 = abs(Nrml[1]); // abs(n*e2)
  double n3 = abs(Nrml[2]); // abs(n*e3)

  vector<double> var(3);
  var[0] = b(acos(n1)/theta);
  var[1] = b(acos(n2)/theta);
  var[2] = b(acos(n3)/theta);
  
  double varsum = var[0] + var[1] + var[2];
  return var[i]/varsum;
}

//*********************************************************************

// Smooth bump function (Beale et al p.4)

double b(double r) {
    
    double bump = 0.0;
    if ( abs(r) < 1.0 ) {
      double r2 = r*r;
      //      double temp = exp(r2/(r2-1.0));
      //      if ( temp > 1.e-15 ) bump = temp;
      bump = exp(r2/(r2-1.0));
    }
    return bump;
}

//*********************************************************************

// Compute surface area (surface integral of function f=1) to test the quadrature

double Compute_surface_area(int N_quad, const vector<Surf_point>& Surfc) {
  double Surf_area = 0.0;
  for (int i=0; i<N_quad; i++) {
        Surf_area += Integrand(Surfc[i].x)*Surfc[i].Area;
  }
  return Surf_area;
}

//*********************************************************************

// Integrand is 1 for computing surface area

double Integrand(const vector<double>& x) {
  return 1.0;
}

//*********************************************************************
//*********************************************************************
//*********************************************************************

// Compute mean curvature at a given point on the surface

double mean_curvature(const vector<double>& x) {

  double phi1,phi2,phi3,phi11,phi12,phi13,phi21,phi22,phi23,phi31,phi32,phi33;
  
  phi1 = Dphi(0,x);
  phi2 = Dphi(1,x);
  phi3 = Dphi(2,x);
  D2phi(x,phi11,phi12,phi13,phi21,phi22,phi23,phi31,phi32,phi33);

  double phi1_sq = phi1*phi1;
  double phi2_sq = phi2*phi2;
  double phi3_sq = phi3*phi3;
  double res;
  res  = phi1_sq*phi22 - 2.0*phi1*phi2*phi12 + phi2_sq*phi11;
  res += phi1_sq*phi33 - 2.0*phi1*phi3*phi13 + phi3_sq*phi11;
  res += phi2_sq*phi33 - 2.0*phi2*phi3*phi23 + phi3_sq*phi22;

  double Dphi_sq = phi1_sq + phi2_sq + phi3_sq;
  res = -res/(2.0*Dphi_sq*sqrt(Dphi_sq));

  return res;
}

//*********************************************************************

// decide which Monge parameterization to use based on components of normal vector
void select_Monge_Patch(int i, const vector<double>& x0, int& i1, int& i2, int& i3) {

  double nx = abs(Dphi(0,x0));
  double ny = abs(Dphi(1,x0));
  double nz = abs(Dphi(2,x0));

  i1 = 0;
  i2 = 0;
  i3 = 0;
  
  if ( nx >= ny ) {
    if ( nx >= nz ) { //normal is mostly in the x-direction
      i1 = 1;
      i2 = 2;
      i3 = 0;
    }
    else { //normal is mostly in the z-direction
      i1 = 0;
      i2 = 1;
      i3 = 2;
    }
  }
  else {
    if ( ny >= nz ) { //normal is mostly in the y-direction
      i1 = 0;
      i2 = 2;
      i3 = 1;
    }
    else { //normal is mostly in the z-direction
      i1 = 0;
      i2 = 1;
      i3 = 2;
    }
  }
}

//*********************************************************************

void generate_Interp_Stencil(double h, double a1, double a2,
			     vector<double>& p1, vector<double>& p2,
			     vector< vector<double> >& M) {
  p1[0] = a1;       p2[0] = a2;
  p1[1] = a1;       p2[1] = a2-h;
  p1[2] = a1;       p2[2] = a2+h;
  p1[3] = a1-h;     p2[3] = a2;
  p1[4] = a1+h;     p2[4] = a2;
  p1[5] = a1-h/2.0; p2[5] = a2-h/2.0;
  // p1[3] = a1-h;     p2[3] = a2-h;
  // p1[4] = a1-h;     p2[4] = a2+h;
  // p1[5] = a1-h/2.0; p2[5] = a2;

  for (int j=0; j<6; j++) {
    M[j][0] = p1[j]*p1[j];
    M[j][1] = p2[j]*p2[j];
    M[j][2] = p1[j]*p2[j];
    M[j][3] = p1[j];
    M[j][4] = p2[j];
    M[j][5] = 1.0;
  }
}

//*********************************************************************

void generate_Interp_Stencil_Cubic(double h, double a1, double a2,
			     vector<double>& p1, vector<double>& p2,
			     vector< vector<double> >& M) {
  p1[0] = a1;       p2[0] = a2;
  p1[1] = a1;       p2[1] = a2-h/2.0;
  //p1[1] = a1;       p2[1] = a2+h/2.0;
  p1[2] = a1;       p2[2] = a2+h;

  p1[3] = a1-h;     p2[3] = a2;
  p1[4] = a1-h;     p2[4] = a2-h;
  p1[5] = a1-h;     p2[5] = a2+h;

  p1[6] = a1+h;     p2[6] = a2;
  p1[7] = a1+h;     p2[7] = a2-h;
  p1[8] = a1+h;     p2[8] = a2+h;
  p1[9] = a1-h/2.0; p2[9] = a2;

  for (int j=0; j<10; j++) {
    M[j][0] = p1[j]*p1[j]*p1[j]; // a1^3
    M[j][1] = p2[j]*p2[j]*p2[j]; // a2^3
    M[j][2] = p1[j]*p1[j]*p2[j]; // a1^2*a2
    M[j][3] = p1[j]*p2[j]*p2[j]; // a1*a2^2
    
    M[j][4] = p1[j]*p1[j]; // a1^2
    M[j][5] = p2[j]*p2[j]; // a2^2
    M[j][6] = p1[j]*p2[j]; // a1*a2

    M[j][7] = p1[j]; // a1
    M[j][8] = p2[j]; // a2
    M[j][9] = 1.0;
  }
}

//*********************************************************************

// find first derivatives using coefficients of interpolation
void first_Derivatives(double a1, double a2, vector<double>& soln,
		       double& D1, double& D2) {
  D1 = 2.0*soln[0]*a1 + soln[2]*a2 + soln[3];
  D2 = 2.0*soln[1]*a2 + soln[2]*a1 + soln[4];
}

void first_Derivatives_Cubic(double a1, double a2, vector<double>& soln,
			     double& D1, double& D2) {
  D1 = 3.0*soln[0]*a1*a1 + 2.0*soln[2]*a1*a2 + soln[3]*a2*a2
     + 2.0*soln[4]*a1 + soln[6]*a2 + soln[7];
  D2 = 3.0*soln[1]*a2*a2 + 2.0*soln[3]*a1*a2 + soln[2]*a1*a1
     + 2.0*soln[5]*a2 + soln[6]*a1 + soln[8];
}

//*********************************************************************

// find second derivatives using coefficients of interpolation
void second_Derivatives(double a1, double a2, vector<double>& soln,
			double& D11, double& D12, double& D22) {
  D11 = 2.0*soln[0];
  D12 = soln[2];
  D22 = 2.0*soln[1];  
}

void second_Derivatives_Cubic(double a1, double a2, vector<double>& soln,
			      double& D11, double& D12, double& D22) {
  D11 = 6.0*soln[0]*a1 + 2.0*soln[2]*a2 + 2.0*soln[4];
  D22 = 6.0*soln[1]*a2 + 2.0*soln[3]*a1 + 2.0*soln[5];  
  D12 = 2.0*soln[2]*a1 + 2.0*soln[3]*a2 + soln[6];
}

//*********************************************************************

// find the third coordinate at interpolation points (solve phi=0 by Newton's method)
void find_Third_Coordinate(int i, int n, int i1, int i2, int i3, double z0, const vector<double>& p1,
			   const vector<double>& p2, vector<double>& p3) {
  vector<double> a2(3,0);
  a2[0] = ellipse_a*ellipse_a;
  a2[1] = ellipse_b*ellipse_b;
  a2[2] = ellipse_c*ellipse_c;

  for (int j=0; j<n; j++) {
    double temp = 1.0 - p1[j]*p1[j]/a2[i1] - p2[j]*p2[j]/a2[i2];
    p3[j] = sqrt(temp*a2[i3]);
    if (z0 < 0.0) p3[j] = -p3[j];
  }
  // p3[0] = z0;
  // for (int j=1; j<n; j++) {
  //   p3[j] = sqrt(1.0 - p1[j]*p1[j] - p2[j]*p2[j]);
  //   if (z0 < 0.0) p3[j] = -p3[j];
  // }
}

//*********************************************************************

void compute_Dual_Tangents(int i1, int i2, int i3, double dhd1, double dhd2,
			   vector<double>& Dual_T1, vector<double>& Dual_T2) {

  // find the inverse metric
  vector<double> ln(2,0);
  vector< vector<double> > inv_g(2,ln);
  double dhd1_sq = dhd1*dhd1;
  double dhd2_sq = dhd2*dhd2;
  double det = 1.0 + dhd1_sq + dhd2_sq;
  inv_g[0][0] = (1.0 + dhd2_sq)/det;
  inv_g[0][1] = -dhd1*dhd2/det;
  inv_g[1][0] = inv_g[0][1];
  inv_g[1][1] = (1.0 + dhd1_sq)/det;
  
  // find the tangents
  vector<double> T1(3,0), T2(3,0);
  T1[i1] = 1.0; T1[i2] = 0.0; T1[i3] = dhd1;
  T2[i1] = 0.0; T2[i2] = 1.0; T2[i3] = dhd2;

  // find the dual tangents
  Dual_T1[i1] = inv_g[0][0]*T1[i1] + inv_g[0][1]*T2[i1];
  Dual_T1[i2] = inv_g[0][0]*T1[i2] + inv_g[0][1]*T2[i2];
  Dual_T1[i3] = inv_g[0][0]*T1[i3] + inv_g[0][1]*T2[i3];
  
  Dual_T2[i1] = inv_g[1][0]*T1[i1] + inv_g[1][1]*T2[i1];
  Dual_T2[i2] = inv_g[1][0]*T1[i2] + inv_g[1][1]*T2[i2];
  Dual_T2[i3] = inv_g[1][0]*T1[i3] + inv_g[1][1]*T2[i3];
      
}

//*********************************************************************
//*********************************************************************

void Generate_Targets_OnSurface(int N_quad, const vector<Surf_point>& Surfc, vector<Target_point>* Target) { 

  Target_point temp;

  for (int i=0; i<N_quad; i++) {
    temp.x = Surfc[i].x[0];
    temp.y = Surfc[i].x[1];
    temp.z = Surfc[i].x[2];
    Target->push_back(temp); //add the point to the target list
  }
}
