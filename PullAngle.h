#ifndef PULLANGLE_H
#define PULLANGLE_H

#include "JetEvent/Jet.h"
#include "CLHEP/Vector/TwoVector.h"
#include "UserAnalysis/Numeric.h"

using namespace Numeric;

class PullAngle { // Helper function
public:
  PullAngle(const CLHEP::Hep2Vector* v): _pull(v) {}

  double operator()(const Jet*& j) const {
    CLHEP::Hep2Vector J = Hep2Vector(j->rapidity(), phiCorr(j->phi()));
    double theta = phiCorr(_pull->angle(J));
    // AxB = |A||B|sin(theta)n^; is n^ pointing up or down?
    // {x1, y1, 0}x{x2, y2, 0} = {0, 0, x1*y2-x2*y1}
    double x1 = j->rapidity();
    double y1 = phiCorr(j->phi());
    double x2 = _pull->x();
    double y2 = phiCorr(_pull->y());
    double Z = (x1 * y2) - (x2 * y1); // out of eta/phi plane, or into?
    return Z >= 0.? theta: -theta;
  }

  double operator()(const Jet* const& j) const {
    CLHEP::Hep2Vector J = Hep2Vector(j->rapidity(), phiCorr(j->phi()));
    double theta = phiCorr(_pull->angle(J));
    // AxB = |A||B|sin(theta)n^; is n^ pointing up or down?
    // {x1, y1, 0}x{x2, y2, 0} = {0, 0, x1*y2-x2*y1}
    double x1 = j->rapidity();
    double y1 = phiCorr(j->phi());
    double x2 = _pull->x();
    double y2 = phiCorr(_pull->y());
    double Z = (x1 * y2) - (x2 * y1); // out of eta/phi plane, or into?
    return Z >= 0.? theta: -theta;
  }

  double operator()(double z = 1.) const { // positive z by default
    CLHEP::Hep2Vector beam = Hep2Vector(z, 0.); // In 2D plane (eta,phi), z is parallel to eta axis
    double theta = phiCorr(_pull->angle(beam));
    // AxB = |A||B|sin(theta)n^; is n^ pointing up or down?
    // {z, 0, 0}x{x, y, 0} = {0, 0, y*z};
    double y = phiCorr(_pull->y());
    double Z = (y * z); // out of eta/phi plane, or into?
    return Z >= 0.? theta: -theta;
  }

private:
  PhiCorr phiCorr;
  const CLHEP::Hep2Vector* _pull;
};

#endif
