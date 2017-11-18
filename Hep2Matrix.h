#ifndef HEP2MATRIX_H
#define HEP2MATRIX_H

#include "CLHEP/Matrix/Matrix.h"
#include "CLHEP/Vector/TwoVector.h"
#include <algorithm> // for sorting eigenvalues
#include <limits> // for NaNs

class Hep2Matrix: public CLHEP::HepMatrix {
public:
  Hep2Matrix(): CLHEP::HepMatrix(2, 2, 0) {}
  Hep2Matrix(int i): CLHEP::HepMatrix(2, 2, i) {}
  Hep2Matrix(double i): CLHEP::HepMatrix(2, 2, i) {}
  ~Hep2Matrix() {}

  const std::pair<double, double> eigenvalues() const {

    const double NaN = std::numeric_limits<double>::quiet_NaN();
    std::pair<double, double> noResult = std::pair<double, double>(NaN, NaN);

    if((*this)[1][0] == 0 && (*this)[0][1] == 0) {
      return std::pair<double, double>((*this)[0][0], (*this)[1][1]);
    }

    double T = (*this)[0][0] - (*this)[1][1];
    double D = T * T + 4.0 * (*this)[0][1] * (*this)[1][0];
    double e[2];
    if(D >= 0) {
      double m = (*this)[0][0] + (*this)[1][1];
      double n = sqrt(D);
      e[0] = 0.5 * (m + n);
      e[1] = 0.5 * (m - n);
    }
    else {
      return noResult;
    }
    
    std::sort(e, e+2, comp); // e1 > e2

    return std::pair<double, double>(e[0], e[1]);
  }
  
  const std::pair<CLHEP::Hep2Vector, CLHEP::Hep2Vector> eigenvectors() const {

    const double NaN = std::numeric_limits<double>::quiet_NaN();
    CLHEP::Hep2Vector vectorNaN = Hep2Vector(NaN, NaN);
    std::pair<CLHEP::Hep2Vector, CLHEP::Hep2Vector> noResult = std::pair<CLHEP::Hep2Vector, CLHEP::Hep2Vector>(vectorNaN, vectorNaN);

    if((*this)[1][0] == 0 && (*this)[0][1] == 0) {
      CLHEP::Hep2Vector evec1 = Hep2Vector(1, 0);
      CLHEP::Hep2Vector evec2 = Hep2Vector(0, 1);
      return std::pair<CLHEP::Hep2Vector, CLHEP::Hep2Vector>(evec1, evec2);
    }

    double e1 = eigenvalues().first;
    double e2 = eigenvalues().second;
    if(e1 != e1 || e2 != e2) return noResult;

    if((*this)[1][0] != 0) {
      CLHEP::Hep2Vector evec1 = Hep2Vector((*this)[1][1] - e1, -(*this)[1][0]);
      CLHEP::Hep2Vector evec2 = Hep2Vector((*this)[1][1] - e2, -(*this)[1][0]);
      return std::pair<CLHEP::Hep2Vector, CLHEP::Hep2Vector>(evec1, evec2);
    }
    
    if((*this)[0][1] != 0) {
      CLHEP::Hep2Vector evec1 = Hep2Vector(-(*this)[0][1], (*this)[0][0] - e1);
      CLHEP::Hep2Vector evec2 = Hep2Vector(-(*this)[0][1], (*this)[0][0] - e2);
      return std::pair<CLHEP::Hep2Vector, CLHEP::Hep2Vector>(evec1, evec2);
    }
  }

private:
  struct Compare {
    bool operator()(double i, double j) { return i > j;}
  } comp;
};

#endif
