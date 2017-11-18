#ifndef EIGENSOLVERSYM3X3_H
#define EIGENSOLVERSYM3X3_H

#include "CLHEP/Matrix/Matrix.h"
#include <algorithm> // for sorting eigenvalues
#include <limits> // For NaN

class EigenSolverSym3x3 {
public:
  EigenSolverSym3x3(CLHEP::HepMatrix A): _A(A) {}

  bool operator()(double root[3]) {
    double max = computeMaximumMagnitudeEntry();
    if(max > 1) _A /= max;
    bool status = computeRoots(root);
    if(max > 1) for(int i = 0; i != 3; i++) root[i] *= max;
    return status;
  }

private:
  double computeMaximumMagnitudeEntry() {
    double max = 0.;
    for(int i = 0; i != _A.num_row(); i++) 
      for(int j = 0; j != _A.num_col(); j++) 
	if(fabs(_A[i][j]) > max) max = fabs(_A[i][j]);
    return max;
  }

  bool computeRoots(double root[3]) {
    const double inv3 = 1./3.;
    const double root3 = sqrt(3.);
    double a00 = _A[0][0];
    double a01 = _A[0][1];
    double a02 = _A[0][2];
    double a11 = _A[1][1];
    double a12 = _A[1][2];
    double a22 = _A[2][2];
    double c0 = a00*a11*a22 + 2.0*a01*a02*a12 - a00*a12*a12 - a11*a02*a02 - a22*a01*a01;
    double c1 = a00*a11 - a01*a01 + a00*a22 - a02*a02 + a11*a22 - a12*a12;
    double c2 = a00 + a11 + a22;
    double c2Div3 = c2*inv3;
    double aDiv3 = (c1 - c2*c2Div3)*inv3;
    if(aDiv3 > 0.0) aDiv3 = 0.0;
    double mbDiv2 = 0.5*(c0 + c2Div3*(2.0*c2Div3*c2Div3 - c1));
    double q = mbDiv2*mbDiv2 + aDiv3*aDiv3*aDiv3;
    if(q > 0.0) q = 0.0;
    double magnitude = -aDiv3 >= 0.? sqrt(-aDiv3): std::numeric_limits<double>::quiet_NaN();
    double angle = -q >= 0.? atan2(sqrt(-q), mbDiv2)*inv3: std::numeric_limits<double>::quiet_NaN();
    double cs = cos(angle);
    double sn = sin(angle);
    root[0] = c2Div3 + 2.0*magnitude*cs;
    root[1] = c2Div3 - magnitude*(cs + root3*sn);
    root[2] = c2Div3 - magnitude*(cs - root3*sn);

    // Sort the roots here to obtain root[0] > root[1] > root[2].
    std::sort(root, root+3, comp);
    return -aDiv3 >= 0. && -q >= 0.;
  }
  CLHEP::HepMatrix _A;

  struct Compare {
    bool operator()(double i, double j) { return i > j;}
  } comp;
};

#endif
