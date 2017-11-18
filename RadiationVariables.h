#ifndef RADIATIONVARIABLES_H
#define RADIATIONVARIABLES_H

#include "JetEvent/Jet.h"
#include "JetEvent/JetConstituentIterator.h"
#include "JetUtils/JetCollectionHelper.h" // For sorting (sub)jets
#include "JetUtils/JetSorters.h" // For sorting (sub)jets

#include "UserAnalysis/Numeric.h" // Common numeric stuff
#include "UserAnalysis/JetShapes.h" // For pull, BCF, etc.
#include "UserAnalysis/PullAngle.h"

#include <vector>

// Note from
// https://cdsweb.cern.ch/record/1418210/files/ATLAS-COM-CONF-2012-002.pdf,
// event selection for ttbar, one of the cuts reads: "a minimum
// distance between the two b-tagged jets ∆R(b, b) > 1.2, to remove bb
// pairs originating from gluon splitting which contribute to the QCD
// multi-jet background at low ∆R whereas the tt signal is
// characterised by the two b-jets with a ∆R around π;"

// debugging macros so we can pin down message provenance at a glance
#include <iostream>
#define DEBUG(x)							\
  std::cout << "DEBUG:\t" << "\t" << #x << "\t" << x << "\t"		\
  << __FUNCTION__ << "\t" << __FILE__ << "\t" << __LINE__ << std::endl;
#define DEBUG(x)

class RadiationVariables {
public:
  RadiationVariables(const JetCollection* theJetCollection):
  _theJetCollection(theJetCollection),
  _Njets(numJets()),
  _dijets(_Njets >= 2) {
    if(_theJetCollection) { DEBUG(_Njets);
      JetCollectionHelper::jetcollection_t theJets(_theJetCollection->begin(), _theJetCollection->end()); DEBUG(0);
      JetCollectionHelper::sort_jets(theJets, JetSorters::sortJetByEtDown()); DEBUG(0);
      JetCollectionHelper::jetcollection_t::iterator firstJet = theJets.begin(); DEBUG(0);
      JetCollectionHelper::jetcollection_t::iterator lastJet  = theJets.end(); DEBUG(0);
      _jets.clear(); DEBUG(0);
      for(unsigned int i = 0; firstJet != lastJet && i != 4; ++firstJet, i++) { DEBUG(i);
    	if((*firstJet) == 0) continue; DEBUG(*firstJet);
    	const Jet* aJet = (*firstJet)->clone(true, true); DEBUG(aJet);
    	if(aJet == 0) continue;; DEBUG(aJet);
    	_jets.push_back(aJet);; DEBUG(_jets.size());
      }
      _firstJet = _Njets >= 1? _jets[0]: 0; DEBUG(_firstJet);
      _secondJet = _Njets >= 2? _jets[1]: 0; DEBUG(_secondJet);
      _thirdJet = _Njets >= 3? _jets[2]: 0; DEBUG(_thirdJet);
      _fourthJet = _Njets >= 4? _jets[3]: 0; DEBUG(_fourthJet);
    }
    DEBUG(0);
  }
  
  ~RadiationVariables() {
    if(_theJetCollection) {
      for(unsigned int i = 0; i != _jets.size(); i++) {
  	delete _jets[i];
  	_jets[i] = 0;
      }
      _jets.clear();
    }
  }

  int numJets() const { DEBUG(_theJetCollection);
    return _theJetCollection? _theJetCollection->size(): 0;
  }

  bool isDijets() const { DEBUG(_dijets); return _dijets; } 

  // Original bcf
  double/* jStr_radEi_*/bcfJ1(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets > 0) {
      BCF bcf(_firstJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcf();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfJ2(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_dijets) {
      BCF bcf(_secondJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcf();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfJ3(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets >= 3) {
      BCF bcf(_thirdJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_thirdJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_thirdJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcf();
    }
    return result;
  }

  // "transverse" bcf
  double/* jStr_radEi_*/bcfTJ1(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets > 0) {
      BCF bcf(_firstJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfT();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfTJ2(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_dijets) {
      BCF bcf(_secondJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfT();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfTJ3(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets >= 3) {
      BCF bcf(_thirdJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_thirdJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_thirdJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfT();
    }
    return result;
  }

  // bcfAsymY
  double/* jStr_radEi_*/bcfAsymYJ1(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets > 0) {
      BCF bcf(_firstJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymY();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfAsymYJ2(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_dijets) {
      BCF bcf(_secondJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymY();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfAsymYJ3(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets >= 3) {
      BCF bcf(_thirdJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_thirdJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_thirdJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymY();
    }
    return result;
  }

  // bcfAsymPhi
  double/* jStr_radEi_*/bcfAsymPhiJ1(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets > 0) {
      BCF bcf(_firstJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymPhi();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfAsymPhiJ2(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_dijets) {
      BCF bcf(_secondJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymPhi();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfAsymPhiJ3(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets >= 3) {
      BCF bcf(_thirdJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_thirdJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_thirdJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymPhi();
    }
    return result;
  }

  // bcfAsymYPhi 
  double/* jStr_radEi_*/bcfAsymYPhiJ1(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets > 0) {
      BCF bcf(_firstJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymYPhi();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfAsymYPhiJ2(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_dijets) {
      BCF bcf(_secondJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymYPhi();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfAsymYPhiJ3(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets >= 3) {
      BCF bcf(_thirdJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_thirdJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_thirdJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymYPhi();
    }
    return result;
  }

  // bcfAsymYPhi2
  double/* jStr_radEi_*/bcfAsymYPhi2J1(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets > 0) {
      BCF bcf(_firstJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymYPhi2();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfAsymYPhi2J2(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_dijets) {
      BCF bcf(_secondJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymYPhi2();
    }
    return result;
  }

  double/* jStr_radEi_*/bcfAsymYPhi2J3(int version = 0, unsigned int a = 1) {
    double result = NaN;
    if(_Njets >= 3) {
      BCF bcf(_thirdJet, version, a);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_thirdJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_thirdJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, bcf(firstConstituent));
      result = bcf.bcfAsymYPhi2();
    }
    return result;
  }

  double/* jStr_radEi_*/BJ1() {
    double result = NaN;
    if(_Njets > 0) {
      Pull pull(_firstJet);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, pull(firstConstituent));
      result = pull.B();
    }
    return result;
  }

  double/* jStr_radEi_*/BJ2() {
    double result = NaN;
    if(_dijets) {
      Pull pull(_secondJet);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, pull(firstConstituent));
      result = pull.B();
    }
    return result;
  }

  double/* jStr_radEi_*/BJ3() {
    double result = NaN;
    if(_Njets >= 3) {
      Pull pull(_thirdJet);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_thirdJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_thirdJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, pull(firstConstituent));
      result = pull.B();
    }
    return result;
  }

  // arXiv: 1010.3698v1 [hep-ph] ATLAS-COM-CONF-2011-056.pdf
  double/* jStr_radEi_*/girthJ1() {
    double result = NaN;
    if(_Njets > 0) {
      Pull pull(_firstJet);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, pull(firstConstituent));
      result = pull.girth();
    }
    return result;
  }
  // ATLAS-COM-CONF-2011-056.pdf
  double/* jStr_radEi_*/girthJ2() {
    double result = NaN;
    if(_dijets) {
      Pull pull(_secondJet);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, pull(firstConstituent));
      result = pull.girth();
    }
    return result;
  }

  double/* jStr_radEi_*/girthJ3() {
    double result = NaN;
    if(_Njets >= 3) {
      Pull pull(_thirdJet);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_thirdJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_thirdJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, pull(firstConstituent));
      result = pull.girth();
    }
    return result;
  }

  double/* jStr_radEi_*/girth32() {
    return Divide(girthJ3(), girthJ2());
  }

  double/* jStr_radEi_*/girth21() {
    return Divide(girthJ2(), girthJ1());
  }

  double/* jStr_radEi_*/girthAsymJ1J2() {
    return Asymmetry(girthJ1(), girthJ2());
  }

  // arXiv:1010.3698v1 [hep-ph]
  double/* jStr_radEi_*/alpha1() {
    double result = NaN;
    if(_dijets) {
      Pull pull(_firstJet);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, pull(firstConstituent));
      PullAngle p(pull.t());
      result = p(_secondJet);
    }
    return result;
  }

  double/* jStr_radEi_*/alpha2() {
    double result = NaN;
    if(_dijets) {
      Pull pull(_secondJet);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, pull(firstConstituent));
      PullAngle p(pull.t());
      result =  p(_firstJet);
    }
    return result;
  }

  double/* jStr_radEi_*/alpha() {
    double result = NaN;
    if(_dijets) {
      double a1 = alpha1();
      double a2 = alpha2();
      result = Sqrt((a1 * a1) + (a2 * a2));
    }
    return result;
  }

  double/* jStr_radEi_*/beta1() { 
    double result = NaN;
    if(_dijets) {
      Pull pull(_firstJet);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, pull(firstConstituent));
      PullAngle p(pull.t());
      double angle1 = p(1.);
      double angle2 = p(-1.);
      result = _firstJet->rapidity() >= 0.? angle1: angle2;
    }
    return result;
  }

  double/* jStr_radEi_*/beta2() {
    double result = NaN;
    if(_dijets) {
      Pull pull(_secondJet);
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent != lastConstituent; ++firstConstituent, pull(firstConstituent));
      PullAngle p(pull.t());
      double angle1 = p(1.);
      double angle2 = p(-1.);
      result = _secondJet->rapidity() >= 0.? angle1: angle2;
    }
    return result;
  }
  
  double/* jStr_radEi_*/beta() {
    double result = NaN;
    if(_dijets) {
      double b1 = beta1();
      double b2 = beta2();
      result = Sqrt((b1 * b1) + (b2 * b2));
    }
    return result;
  }

  // // Should provide more discrimination than alphaX and betaX alone
  // double gamma1() const {
  //   return _dijets? Divide(fabs(alpha1()), fabs(beta1())): NaN;
  // }

  // double gamma2() const {
  //   return _dijets? Divide(fabs(alpha2()), fabs(beta2())): NaN;
  // }

  // double gamma() const {
  //   double result = NaN;
  //   if(_dijets) {
  //     double g1 = gamma1();
  //     double g2 = gamma2();
  //     result = Sqrt((g1 * g1) + (g2 * g2));
  //   }
  //   return result;
  // }

  // My idea: delta between jet ellipse orientations
  double/* jStr_radEi_*/thetaJ1J2() {
    double result = NaN;
    if(_dijets) {
      DeltaPhi deltaPhi;
      Pull pull1(_firstJet);
      JetConstituentIterator firstConstituent1 = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent1 = JetConstituentIterator::last(_firstJet);
      for(; firstConstituent1 != lastConstituent1; ++firstConstituent1, pull1(firstConstituent1));
      Pull pull2(_secondJet);
      JetConstituentIterator firstConstituent2 = JetConstituentIterator::first(_secondJet);
      JetConstituentIterator lastConstituent2 = JetConstituentIterator::last(_secondJet);
      for(; firstConstituent2 != lastConstituent2; ++firstConstituent2, pull2(firstConstituent2));
      double theta1 = pull1.orientation();
      double theta2 = pull2.orientation();
      double thetaJ1J2 = deltaPhi(theta1, theta2);
      result = thetaJ1J2;
    }
    return result;
  }

  // arXiv: 1102.1012v1 [hep-ph], promoted to be an event-shape variable
  double/* jStr_radEi_*/dipolarityInfLine() const { // dipolarityPt
    double result = NaN;
    if(_dijets) {
      DeltaPhi deltaPhi;
      if(_firstJet == 0) return NaN;
      Jet* fatJet = _firstJet->clone(true, true);
      if(fatJet == 0) return NaN;

      double rapidity1 = _firstJet->rapidity(); // x1
      double rapidity2 = _secondJet->rapidity(); // x2
      double phi1 = _firstJet->phi(); // y1
      double phi2 = _secondJet->phi(); // y2
      double dRapidity = rapidity2 - rapidity1; // ux = x2 - x1
      double dPhi = deltaPhi(phi2, phi1); // uy = y2 - y1
      double jetSep2 = (dRapidity * dRapidity) + (dPhi * dPhi); // length = ux2 + uy2
      double jetSep = sqrt(jetSep2); // sqrt(length)

      double sum = 0;
      bool status = fatJet->addJet(_secondJet);
      if(status) {
	double et_J = fatJet->et();

	JetConstituentIterator firstConstituent = JetConstituentIterator::first(fatJet);
	JetConstituentIterator lastConstituent = JetConstituentIterator::last(fatJet);
      
	for(; firstConstituent != lastConstituent; ++firstConstituent) {
	  double et_i = firstConstituent.et();

	  double rapidity_i = firstConstituent.rapidity();
	  double dRapidity_i = rapidity1 - rapidity_i; // vx = x1 - X
	  double phi_i = firstConstituent.phi();
	  double dPhi_i = deltaPhi(phi1, phi_i); // vy = y1 - Y

	  double Ri;
	  //double det = (dRapidity_i * dRapidity) + (dPhi_i * dPhi); // det = (vx * ux) + (vy * uy) = v.u = |v||u|cos(theta)

	  if(false) {
	    //if(det < 0 || det > jetSep) {
	    dRapidity = rapidity2 - rapidity_i; // ux = x2 - X
	    dPhi = deltaPhi(phi2, phi_i); // uy = y2 - Y
	    Ri = sqrt(std::min(
			       (dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i), // vx2 + vy2
			       (dRapidity * dRapidity) + (dPhi * dPhi) // ux2 + uy2
			       ));
	  }
	  else {
	    Ri = Divide(fabs(dRapidity_i * dPhi - dRapidity * dPhi_i), jetSep); // |vx * uy - ux * vy| / length
	  }

	  sum += Divide(et_i, et_J) * (Ri * Ri);
	}
	result = Divide(sum, jetSep2);
      }
      delete fatJet;
    }
    return result;
  }

  double/* jStr_radEi_*/dipolarityLineSeg() const { // dipolarityEt
    double result = NaN;
    if(_dijets) {
      DeltaPhi deltaPhi;
      if(_firstJet == 0) return NaN;
      Jet* fatJet = _firstJet->clone(true, true);
      if(fatJet == 0) return NaN;

      double rapidity1 = _firstJet->rapidity(); // x1
      double rapidity2 = _secondJet->rapidity(); // x2
      double phi1 = _firstJet->phi(); // y1
      double phi2 = _secondJet->phi(); // y2
      double dRapidity = rapidity2 - rapidity1; // ux = x2 - x1
      double dPhi = deltaPhi(phi2, phi1); // uy = y2 - y1
      double jetSep2 = (dRapidity * dRapidity) + (dPhi * dPhi); // length = ux2 + uy2
      double jetSep = sqrt(jetSep2); // sqrt(length)

      double sum = 0;
      bool status = fatJet->addJet(_secondJet);
      if(status) {
	double et_J = fatJet->et();

	JetConstituentIterator firstConstituent = JetConstituentIterator::first(fatJet);
	JetConstituentIterator lastConstituent = JetConstituentIterator::last(fatJet);
      
	for(; firstConstituent != lastConstituent; ++firstConstituent) {
	  double et_i = firstConstituent.et();

	  double rapidity_i = firstConstituent.rapidity();
	  double dRapidity_i = rapidity1 - rapidity_i; // vx = x1 - X
	  double phi_i = firstConstituent.phi();
	  double dPhi_i = deltaPhi(phi1, phi_i); // vy = y1 - Y

	  double Ri;
	  double det = (dRapidity_i * dRapidity) + (dPhi_i * dPhi); // det = (vx * ux) + (vy * uy) = v.u = |v||u|cos(theta)

	  if(det < 0 || det > jetSep) {
	    dRapidity = rapidity2 - rapidity_i; // ux = x2 - X
	    dPhi = deltaPhi(phi2, phi_i); // uy = y2 - Y
	    // Ri = sqrt(std::min(
	    // 		       (dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i), // vx2 + vy2
	    // 		       (dRapidity * dRapidity) + (dPhi * dPhi) // ux2 + uy2
	    // 		       ));
	    Ri = 0.;
	  }
	  else {
	    Ri = Divide(fabs(dRapidity_i * dPhi - dRapidity * dPhi_i), jetSep); // |vx * uy - ux * vy| / length
	  }

	  sum += Divide(et_i, et_J) * (Ri * Ri);
	}
	result = Divide(sum, jetSep2);
      }
      delete fatJet;
    }
    return result;
  }

  double/* jStr_radEi_*/BCF1() const { // dipolarityPtY // infinite horizontal lines (-inf, inf)
    double result = NaN;
    if(_dijets) {
      double maxRapidity = 1000.;
      DeltaPhi deltaPhi;
      double sum = 0;

      // first jet
      double rapidity1 = _firstJet->rapidity(); // x1
      double rapidity2 = rapidity1 >= 0.? maxRapidity: -maxRapidity; // x2
      double phi1 = _firstJet->phi(); // y1 (phi2 = phi1 => dPhi = uy = 0)
      double dRapidity = rapidity2 - rapidity1; // ux = x2 - x1

      double et_J = _firstJet->et();
      
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      
      for(; firstConstituent != lastConstituent; ++firstConstituent) {
	double et_i = firstConstituent.et();

	//double rapidity_i = firstConstituent.rapidity();
	//double dRapidity_i = rapidity1 - rapidity_i; // vx = x1 - X
	double phi_i = firstConstituent.phi();
	double dPhi_i = deltaPhi(phi1, phi_i); // vy = y1 - Y
	
	double Ri;
	//double det = (dRapidity_i * dRapidity); // det = (vx * ux) + (vy * uy) = v.u = |v||u| cos(theta)	
	
	if(false) { // will never do this!
	  //Ri = sqrt((dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i));
	  Ri = 0.;
	}
	else {
	  Ri = Divide(fabs(dRapidity * dPhi_i), fabs(dRapidity)); // |vx * uy - ux * vy| / length
	}
	sum += Divide(et_i, et_J) * (Ri * Ri);
      }

      // second jet
      rapidity1 = _secondJet->rapidity(); // x1
      rapidity2 = rapidity1 >= 0.? maxRapidity: -maxRapidity; // x2
      phi1 = _secondJet->phi(); // y1 (phi2 = phi1 => dPhi = uy = 0)
      dRapidity = rapidity2 - rapidity1; // ux = x2 - x1

      et_J = _secondJet->et();
      
      firstConstituent = JetConstituentIterator::first(_secondJet);
      lastConstituent = JetConstituentIterator::last(_secondJet);
      
      for(; firstConstituent != lastConstituent; ++firstConstituent) {
	double et_i = firstConstituent.et();

	//double rapidity_i = firstConstituent.rapidity();
	//double dRapidity_i = rapidity1 - rapidity_i; // vx = x1 - X
	double phi_i = firstConstituent.phi();
	double dPhi_i = deltaPhi(phi1, phi_i); // vy = y1 - Y
	
	double Ri;
	//double det = (dRapidity_i * dRapidity); // det = (vx * ux) + (vy * uy) = v.u = |v||u| cos(theta)	
	
	if(false) { // will never do this!
	  //Ri = sqrt((dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i)); 
	  Ri = 0.;
	}
	else {
	  Ri = Divide(fabs(dRapidity * dPhi_i), fabs(dRapidity)); // |vx * uy - ux * vy| / length
	}
	sum += Divide(et_i, et_J) * (Ri * Ri);
      }

      rapidity1 = _firstJet->rapidity(); // x1
      rapidity2 = _secondJet->rapidity(); // x2
      phi1 = _firstJet->phi(); // y1
      double phi2 = _secondJet->phi(); // y2
      dRapidity = rapidity2 - rapidity1; // ux = x2 - x1
      double dPhi = deltaPhi(phi2, phi1); // uy = y2 - y1
      double jetSep2 = (dRapidity * dRapidity) + (dPhi * dPhi); // length = ux2 + uy2
      //double jetSep = sqrt(jetSep2); // sqrt(length)

      result = Divide(sum, jetSep2);
    }
    return result;
  }

  double/* jStr_radEi_*/BCF2() const { // dipolarityPtY // infinite horizontal lines [eta, inf); (-inf, -eta]
    double result = NaN;
    if(_dijets) {
      double maxRapidity = 1000.;
      DeltaPhi deltaPhi;
      double sum = 0;

      // first jet
      double rapidity1 = _firstJet->rapidity(); // x1
      double rapidity2 = rapidity1 >= 0.? maxRapidity: -maxRapidity; // x2
      double phi1 = _firstJet->phi(); // y1 (phi2 = phi1 => dPhi = uy = 0)
      double dRapidity = rapidity2 - rapidity1; // ux = x2 - x1

      double et_J = _firstJet->et();
      
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      
      for(; firstConstituent != lastConstituent; ++firstConstituent) {
	double et_i = firstConstituent.et();

	double rapidity_i = firstConstituent.rapidity();
	double dRapidity_i = rapidity1 - rapidity_i; // vx = x1 - X
	double phi_i = firstConstituent.phi();
	double dPhi_i = deltaPhi(phi1, phi_i); // vy = y1 - Y
	
	double Ri;
	double det = (dRapidity_i * dRapidity); // det = (vx * ux) + (vy * uy) = v.u = |v||u| cos(theta)	
	
	if(det < 0) { // we're off the end of the line
	  //Ri = sqrt((dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i));
	  Ri = 0.;
	}
	else {
	  Ri = Divide(fabs(dRapidity * dPhi_i), fabs(dRapidity)); // |vx * uy - ux * vy| / length
	}
	sum += Divide(et_i, et_J) * (Ri * Ri);
      }

      // second jet
      rapidity1 = _secondJet->rapidity(); // x1
      rapidity2 = rapidity1 >= 0.? maxRapidity: -maxRapidity; // x2
      phi1 = _secondJet->phi(); // y1 (phi2 = phi1 => dPhi = uy = 0)
      dRapidity = rapidity2 - rapidity1; // ux = x2 - x1

      et_J = _secondJet->et();
      
      firstConstituent = JetConstituentIterator::first(_secondJet);
      lastConstituent = JetConstituentIterator::last(_secondJet);
      
      for(; firstConstituent != lastConstituent; ++firstConstituent) {
	double et_i = firstConstituent.et();

	double rapidity_i = firstConstituent.rapidity();
	double dRapidity_i = rapidity1 - rapidity_i; // vx = x1 - X
	double phi_i = firstConstituent.phi();
	double dPhi_i = deltaPhi(phi1, phi_i); // vy = y1 - Y
	
	double Ri;
	double det = (dRapidity_i * dRapidity); // det = (vx * ux) + (vy * uy) = v.u = |v||u| cos(theta)	
	
	if(det < 0) { // we're off the end of the line
	  //Ri = sqrt((dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i));
	  Ri = 0.;
	}
	else {
	  Ri = Divide(fabs(dRapidity * dPhi_i), fabs(dRapidity)); // |vx * uy - ux * vy| / length
	}
	sum += Divide(et_i, et_J) * (Ri * Ri);
      }

      rapidity1 = _firstJet->rapidity(); // x1
      rapidity2 = _secondJet->rapidity(); // x2
      phi1 = _firstJet->phi(); // y1
      double phi2 = _secondJet->phi(); // y2
      dRapidity = rapidity2 - rapidity1; // ux = x2 - x1
      double dPhi = deltaPhi(phi2, phi1); // uy = y2 - y1
      double jetSep2 = (dRapidity * dRapidity) + (dPhi * dPhi); // length = ux2 + uy2
      //double jetSep = sqrt(jetSep2); // sqrt(length)

      result = Divide(sum, jetSep2);
    }
    return result;
  }

  double/* jStr_radEi_*/BCF3() const { // dipolarityPtY // infinite horizontal lines [eta, inf); (-inf, -eta]
    double result = NaN;
    if(_dijets) {
      double maxRapidity = 1000.;
      DeltaPhi deltaPhi;
      double sum = 0;

      // first jet
      double rapidity1 = _firstJet->rapidity(); // x1
      double rapidity2 = rapidity1 >= 0.? maxRapidity: -maxRapidity; // x2
      double phi1 = _firstJet->phi(); // y1 (phi2 = phi1 => dPhi = uy = 0)
      double dRapidity = rapidity2 - rapidity1; // ux = x2 - x1

      double et_J = _firstJet->et();
      
      JetConstituentIterator firstConstituent = JetConstituentIterator::first(_firstJet);
      JetConstituentIterator lastConstituent = JetConstituentIterator::last(_firstJet);
      
      for(; firstConstituent != lastConstituent; ++firstConstituent) {
	double et_i = firstConstituent.et();

	double rapidity_i = firstConstituent.rapidity();
	double dRapidity_i = rapidity1 - rapidity_i; // vx = x1 - X
	double phi_i = firstConstituent.phi();
	double dPhi_i = deltaPhi(phi1, phi_i); // vy = y1 - Y
	
	double Ri;
	double det = (dRapidity_i * dRapidity); // det = (vx * ux) + (vy * uy) = v.u = |v||u| cos(theta)	
	
	if(det < 0) { // we're off the end of the line
	  Ri = sqrt((dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i));
	  //Ri = 0.;
	}
	else {
	  Ri = Divide(fabs(dRapidity * dPhi_i), fabs(dRapidity)); // |vx * uy - ux * vy| / length
	}
	sum += Divide(et_i, et_J) * (Ri * Ri);
      }

      // second jet
      rapidity1 = _secondJet->rapidity(); // x1
      rapidity2 = rapidity1 >= 0.? maxRapidity: -maxRapidity; // x2
      phi1 = _secondJet->phi(); // y1 (phi2 = phi1 => dPhi = uy = 0)
      dRapidity = rapidity2 - rapidity1; // ux = x2 - x1

      et_J = _secondJet->et();
      
      firstConstituent = JetConstituentIterator::first(_secondJet);
      lastConstituent = JetConstituentIterator::last(_secondJet);
      
      for(; firstConstituent != lastConstituent; ++firstConstituent) {
	double et_i = firstConstituent.et();

	double rapidity_i = firstConstituent.rapidity();
	double dRapidity_i = rapidity1 - rapidity_i; // vx = x1 - X
	double phi_i = firstConstituent.phi();
	double dPhi_i = deltaPhi(phi1, phi_i); // vy = y1 - Y
	
	double Ri;
	double det = (dRapidity_i * dRapidity); // det = (vx * ux) + (vy * uy) = v.u = |v||u| cos(theta)	
	
	if(det < 0) { // we're off the end of the line
	  Ri = sqrt((dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i));
	  //Ri = 0.;
	}
	else {
	  Ri = Divide(fabs(dRapidity * dPhi_i), fabs(dRapidity)); // |vx * uy - ux * vy| / length
	}
	sum += Divide(et_i, et_J) * (Ri * Ri);
      }

      rapidity1 = _firstJet->rapidity(); // x1
      rapidity2 = _secondJet->rapidity(); // x2
      phi1 = _firstJet->phi(); // y1
      double phi2 = _secondJet->phi(); // y2
      dRapidity = rapidity2 - rapidity1; // ux = x2 - x1
      double dPhi = deltaPhi(phi2, phi1); // uy = y2 - y1
      double jetSep2 = (dRapidity * dRapidity) + (dPhi * dPhi); // length = ux2 + uy2
      //double jetSep = sqrt(jetSep2); // sqrt(length)

      result = Divide(sum, jetSep2);
    }
    return result;
  }

  double/* jStr_radEi_*/dipolarity() const { // dipolarityEtY, i.e. the original meaning
    // http://www.codeguru.com/forum/showthread.php?t=194400&page=2
    // http://mathworld.wolfram.com/Point-LineDistance2-Dimensional.html
    // http://stackoverflow.com/questions/849211/shortest-distance-between-a-point-and-a-line-segment
    
    double result = NaN;
    if(_dijets) {
      DeltaPhi deltaPhi;
      if(_firstJet == 0) return NaN;
      Jet* fatJet = _firstJet->clone(true, true);
      if(fatJet == 0) return NaN;

      double rapidity1 = _firstJet->rapidity(); // x1
      double rapidity2 = _secondJet->rapidity(); // x2
      double phi1 = _firstJet->phi(); // y1
      double phi2 = _secondJet->phi(); // y2
      double dRapidity = rapidity2 - rapidity1; // ux = x2 - x1
      double dPhi = deltaPhi(phi2, phi1); // uy = y2 - y1
      double jetSep2 = (dRapidity * dRapidity) + (dPhi * dPhi); // length = ux2 + uy2
      double jetSep = sqrt(jetSep2); // sqrt(length)

      double sum = 0;
      bool status = fatJet->addJet(_secondJet);
      if(status) {
	double et_J = fatJet->et();

	JetConstituentIterator firstConstituent = JetConstituentIterator::first(fatJet);
	JetConstituentIterator lastConstituent = JetConstituentIterator::last(fatJet);
      
	for(; firstConstituent != lastConstituent; ++firstConstituent) {
	  double et_i = firstConstituent.et();

	  double rapidity_i = firstConstituent.rapidity();
	  double dRapidity_i = rapidity1 - rapidity_i; // vx = x1 - X
	  double phi_i = firstConstituent.phi();
	  double dPhi_i = deltaPhi(phi1, phi_i); // vy = y1 - Y

	  double Ri;
	  double det = (dRapidity_i * dRapidity) + (dPhi_i * dPhi); // det = (vx * ux) + (vy * uy) = v.u = |v||u|cos(theta)

	  if(det < 0 || det > jetSep) {
	    dRapidity = rapidity2 - rapidity_i; // ux = x2 - X
	    dPhi = deltaPhi(phi2, phi_i); // uy = y2 - Y
	    Ri = sqrt(std::min(
			       (dRapidity_i * dRapidity_i) + (dPhi_i * dPhi_i), // vx2 + vy2
			       (dRapidity * dRapidity) + (dPhi * dPhi) // ux2 + uy2
			       ));
	  }
	  else {
	    Ri = Divide(fabs(dRapidity_i * dPhi - dRapidity * dPhi_i), jetSep); // |vx * uy - ux * vy| / length
	  }

	  sum += Divide(et_i, et_J) * (Ri * Ri);
	}
	result = Divide(sum, jetSep2);
      }
      delete fatJet;
    }
    return result;
  }
  
private:
  const JetCollection* _theJetCollection;
  std::vector<const Jet*> _jets;
  const unsigned int _Njets;
  const bool _dijets;
  const Jet* _firstJet;
  const Jet* _secondJet;
  const Jet* _thirdJet;
  const Jet* _fourthJet;
};

#endif
